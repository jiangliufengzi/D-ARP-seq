#!/usr/bin/env python3
"""
================================================================================
Script: annotate_tRNA_position_sprinzl.py
Purpose: Add Sprinzl numbering (standard tRNA site names) to chr/position tables
================================================================================

Background:
  Map raw coordinates from a tRNA mutation/modification pipeline
  (tRNAscan-SE ID + raw position) onto Sprinzl numbering so sites can be
  compared across studies and with the literature.

  Typical inputs:
    - TSV tables from tRNA mutation calling (chr/trnascan_id and position/raw_pos)
    - Site summaries that need Sprinzl annotation

  Outputs:
    - Insert sprinzl_position after the position column
    - Count annotated vs unmatched rows
    - Unmatched sites keep the original position as a fallback

  Sprinzl numbering:
    Community-standard tRNA site numbers, anchored on the anticodon loop
    (positions 34/35/36), so homologous sites can be compared across tRNAs.

Workflow:
  parse args (in-place and --output are mutually exclusive)
    -> load_mapping() from hg38_trna_raw_to_sprinzl.tsv
    -> resolve_jobs() (single file, directory, or --in-place)
    -> annotate_file() / write_annotated_file() (optional ProcessPool)
    -> write {stem}.sprinzl.tsv, sprinzl_annotated/, or overwrite in place

Mapping lookup in write_annotated_file():
  1. Locate chr/position columns (names are configurable)
  2. Drop an existing output_col to avoid duplicate inserts
  3. Insert sprinzl_position after position
  4. key = (normalize_id(chr), normalize_position(position))
  5. Hit -> write sprinzl_pos; miss -> keep the raw position

In-place writes:
  write .{name}.sprinzl.tmp.{pid}, then os.replace() onto the original file

CLI (summary):
  -i/--input, -o/--output, --in-place, --suffix
  --map-file, --map-id-col, --map-position-col, --map-sprinzl-col
  --chr-col, --position-col, --output-col, --sep, --map-sep
  --pattern, --recursive, --workers, --continue-on-error

Notes:
  1. Mapping keys must match input chr/position exactly
  2. Unmatched sites are not errors; they keep the raw position
  3. --in-place uses a temp file plus atomic replace
  4. Directory mode skips the mapping file itself
  5. Conflicting Sprinzl values for the same key abort the run

================================================================================
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_MAP_FILE = SCRIPT_DIR / "hg38_trna_raw_to_sprinzl.tsv"

_WORKER_MAPPING = None
_WORKER_CONFIG = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Add Sprinzl position annotation to tables with chr and position columns.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    input_group = parser.add_argument_group("Input / output")
    input_group.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input file or directory.",
    )
    input_group.add_argument(
        "-o",
        "--output",
        default=None,
        help=(
            "Output file or directory. Optional for a single input file "
            "(a suffix is inserted into the original name); optional for a "
            "directory (defaults to sprinzl_annotated/ under the input directory)."
        ),
    )
    input_group.add_argument(
        "--in-place",
        action="store_true",
        help="Overwrite input files in place. Cannot be combined with --output.",
    )
    input_group.add_argument(
        "--suffix",
        default=".sprinzl",
        help="Suffix inserted into output names when not overwriting in place (default: .sprinzl).",
    )

    map_group = parser.add_argument_group("Sprinzl mapping table")
    map_group.add_argument(
        "--map-file",
        default=str(DEFAULT_MAP_FILE),
        help=f"Path to the Sprinzl mapping table (default: {DEFAULT_MAP_FILE}).",
    )
    map_group.add_argument(
        "--map-id-col",
        default="trnascan_id",
        help="Mapping-table column matching the input chr column (default: trnascan_id).",
    )
    map_group.add_argument(
        "--map-position-col",
        default="raw_pos",
        help="Mapping-table column matching the input position column (default: raw_pos).",
    )
    map_group.add_argument(
        "--map-sprinzl-col",
        default="sprinzl_pos",
        help="Mapping-table column written to the output (default: sprinzl_pos).",
    )

    column_group = parser.add_argument_group("Input column names")
    column_group.add_argument(
        "--chr-col",
        default="chr",
        help="ID/chr column in the input file (default: chr).",
    )
    column_group.add_argument(
        "--position-col",
        default="position",
        help="Position column in the input file (default: position).",
    )
    column_group.add_argument(
        "--output-col",
        default="sprinzl_position",
        help="Name of the added/replaced output column (default: sprinzl_position).",
    )

    format_group = parser.add_argument_group("File format")
    format_group.add_argument(
        "--sep",
        default="tab",
        help="Input/output delimiter: tab, comma, semicolon, pipe, or a single character (default: tab).",
    )
    format_group.add_argument(
        "--map-sep",
        default="tab",
        help="Mapping-table delimiter: tab, comma, semicolon, pipe, or a single character (default: tab).",
    )

    batch_group = parser.add_argument_group("Directory batch mode")
    batch_group.add_argument(
        "--pattern",
        action="append",
        default=None,
        help="Glob pattern(s) to process in directory mode; may be repeated (default: *.tsv).",
    )
    batch_group.add_argument(
        "--recursive",
        action="store_true",
        help="Search input files recursively in directory mode.",
    )
    batch_group.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Process-pool size in directory mode (default: 1).",
    )
    batch_group.add_argument(
        "--continue-on-error",
        action="store_true",
        help="In directory mode, continue after a per-file failure and report all failures at the end.",
    )

    args = parser.parse_args()

    if args.in_place and args.output:
        parser.error("--in-place cannot be combined with --output.")
    if args.workers < 1:
        parser.error("--workers must be >= 1.")

    args.sep = parse_delimiter(args.sep, "--sep")
    args.map_sep = parse_delimiter(args.map_sep, "--map-sep")
    args.pattern = args.pattern or ["*.tsv"]

    return args


def parse_delimiter(value: str, arg_name: str) -> str:
    aliases = {
        "tab": "\t",
        "\\t": "\t",
        "t": "\t",
        "comma": ",",
        ",": ",",
        "semicolon": ";",
        ";": ";",
        "pipe": "|",
        "|": "|",
        "space": " ",
    }
    delimiter = aliases.get(value, value)
    if len(delimiter) != 1:
        raise SystemExit(f"[ERROR] {arg_name} must be a single character or tab/comma/semicolon/pipe.")
    return delimiter


def open_text(path: Path, mode: str):
    if "b" in mode:
        raise ValueError("open_text only supports text mode")
    encoding = "utf-8" if any(flag in mode for flag in ("w", "a", "x")) else "utf-8-sig"
    if path.name.endswith(".gz"):
        return gzip.open(path, mode, encoding=encoding, newline="")
    return open(path, mode, encoding=encoding, newline="")


def normalize_position(value: str) -> str:
    text = str(value).strip()
    if re.fullmatch(r"[+-]?\d+", text):
        return str(int(text))
    if re.fullmatch(r"[+-]?\d+\.0+", text):
        return str(int(text.split(".", 1)[0]))
    return text


def normalize_id(value: str) -> str:
    return str(value).strip()


def load_mapping(
    map_file: Path,
    map_sep: str,
    map_id_col: str,
    map_position_col: str,
    map_sprinzl_col: str,
) -> dict[tuple[str, str], str]:
    if not map_file.is_file():
        raise FileNotFoundError(f"Mapping table not found: {map_file}")

    mapping: dict[tuple[str, str], str] = {}
    row_count = 0
    duplicate_count = 0

    with open_text(map_file, "rt") as handle:
        reader = csv.DictReader(handle, delimiter=map_sep)
        if reader.fieldnames is None:
            raise ValueError(f"Mapping table is empty: {map_file}")

        required_cols = [map_id_col, map_position_col, map_sprinzl_col]
        missing_cols = [col for col in required_cols if col not in reader.fieldnames]
        if missing_cols:
            raise ValueError(
                f"Mapping table is missing columns: {', '.join(missing_cols)}; "
                f"found: {', '.join(reader.fieldnames)}"
            )

        for row in reader:
            row_count += 1
            key = (
                normalize_id(row[map_id_col]),
                normalize_position(row[map_position_col]),
            )
            value = str(row[map_sprinzl_col]).strip()
            if key in mapping:
                if mapping[key] != value:
                    raise ValueError(
                        "Mapping table has conflicting Sprinzl values for the same ID/position: "
                        f"{key[0]} {key[1]} -> {mapping[key]} / {value}"
                    )
                duplicate_count += 1
                continue
            mapping[key] = value

    if not mapping:
        raise ValueError(f"Mapping table has no usable records: {map_file}")

    print(
        f"[INFO] Loaded Sprinzl mapping: {len(mapping)} keys "
        f"({row_count} rows, {duplicate_count} duplicate consistent keys)"
    )
    return mapping


def insert_suffix(path: Path, suffix: str) -> Path:
    if path.name.endswith(".gz"):
        without_gz = path.with_suffix("")
        new_name = f"{without_gz.stem}{suffix}{without_gz.suffix}.gz"
        return path.with_name(new_name)
    return path.with_name(f"{path.stem}{suffix}{path.suffix}")


def collect_input_files(
    input_dir: Path,
    patterns: Iterable[str],
    recursive: bool,
    map_file: Path,
) -> list[Path]:
    files: list[Path] = []
    seen: set[Path] = set()
    resolved_map = map_file.resolve()

    for pattern in patterns:
        iterator = input_dir.rglob(pattern) if recursive else input_dir.glob(pattern)
        for path in sorted(iterator):
            if not path.is_file():
                continue
            resolved = path.resolve()
            if resolved == resolved_map:
                continue
            if resolved in seen:
                continue
            files.append(path)
            seen.add(resolved)

    return files


def resolve_jobs(args: argparse.Namespace) -> tuple[list[tuple[Path, Path]], bool]:
    input_path = Path(args.input).resolve()

    if not input_path.exists():
        raise FileNotFoundError(f"Input path not found: {input_path}")

    if input_path.is_file():
        if args.in_place:
            return [(input_path, input_path)], False

        if args.output:
            output_path = Path(args.output).resolve()
            if output_path.exists() and output_path.is_dir():
                output_path = output_path / insert_suffix(input_path, args.suffix).name
        else:
            output_path = insert_suffix(input_path, args.suffix)
        return [(input_path, output_path)], False

    if not input_path.is_dir():
        raise ValueError(f"Input path is neither a file nor a directory: {input_path}")

    input_files = collect_input_files(
        input_path,
        args.pattern,
        args.recursive,
        Path(args.map_file).resolve(),
    )
    if not input_files:
        raise ValueError(
            f"No files matching {', '.join(args.pattern)} under {input_path}"
        )

    if args.in_place:
        return [(path, path) for path in input_files], True

    output_dir = (
        Path(args.output).resolve()
        if args.output
        else input_path / "sprinzl_annotated"
    )
    if output_dir.exists() and not output_dir.is_dir():
        raise ValueError(f"For a directory input, --output must be a directory: {output_dir}")

    jobs: list[tuple[Path, Path]] = []
    for input_file in input_files:
        rel_path = input_file.relative_to(input_path)
        output_file = output_dir / rel_path.parent / insert_suffix(rel_path, args.suffix).name
        jobs.append((input_file, output_file))

    return jobs, True


def build_config(args: argparse.Namespace) -> dict[str, str]:
    return {
        "chr_col": args.chr_col,
        "position_col": args.position_col,
        "output_col": args.output_col,
        "sep": args.sep,
    }


def _init_worker(mapping: dict[tuple[str, str], str], config: dict[str, str]) -> None:
    global _WORKER_MAPPING, _WORKER_CONFIG
    _WORKER_MAPPING = mapping
    _WORKER_CONFIG = config


def _worker_annotate(job: tuple[Path, Path]) -> dict[str, object]:
    if _WORKER_MAPPING is None or _WORKER_CONFIG is None:
        raise RuntimeError("worker is not initialized")
    input_file, output_file = job
    return annotate_file(input_file, output_file, _WORKER_MAPPING, **_WORKER_CONFIG)


def annotate_file(
    input_file: Path,
    output_file: Path,
    mapping: dict[tuple[str, str], str],
    chr_col: str,
    position_col: str,
    output_col: str,
    sep: str,
) -> dict[str, object]:
    output_file.parent.mkdir(parents=True, exist_ok=True)

    same_file = input_file.resolve() == output_file.resolve()
    write_target = output_file
    if same_file:
        write_target = input_file.with_name(
            f".{input_file.name}.sprinzl.tmp.{os.getpid()}"
        )

    try:
        result = write_annotated_file(
            input_file=input_file,
            output_file=write_target,
            mapping=mapping,
            chr_col=chr_col,
            position_col=position_col,
            output_col=output_col,
            sep=sep,
        )
        if same_file:
            os.replace(write_target, input_file)
            result["output"] = str(input_file)
        return result
    except Exception:
        if same_file and write_target.exists():
            write_target.unlink()
        raise


def write_annotated_file(
    input_file: Path,
    output_file: Path,
    mapping: dict[tuple[str, str], str],
    chr_col: str,
    position_col: str,
    output_col: str,
    sep: str,
) -> dict[str, object]:
    total_rows = 0
    annotated_rows = 0

    with open_text(input_file, "rt") as in_handle, open_text(output_file, "wt") as out_handle:
        reader = csv.reader(in_handle, delimiter=sep)
        writer = csv.writer(out_handle, delimiter=sep, lineterminator="\n")

        try:
            header = next(reader)
        except StopIteration:
            raise ValueError(f"Input file is empty: {input_file}")

        if chr_col not in header:
            raise ValueError(f"{input_file} is missing column {chr_col}; found: {', '.join(header)}")
        if position_col not in header:
            raise ValueError(
                f"{input_file} is missing column {position_col}; found: {', '.join(header)}"
            )

        chr_idx = header.index(chr_col)
        position_idx = header.index(position_col)

        existing_output_indices = {
            idx for idx, col_name in enumerate(header) if col_name == output_col
        }
        base_header = [
            col_name
            for idx, col_name in enumerate(header)
            if idx not in existing_output_indices
        ]
        insert_idx = base_header.index(position_col) + 1
        new_header = base_header[:insert_idx] + [output_col] + base_header[insert_idx:]
        writer.writerow(new_header)

        base_header_len = len(base_header)
        for row in reader:
            total_rows += 1

            if len(row) < len(header):
                row = row + [""] * (len(header) - len(row))

            chr_value = row[chr_idx] if chr_idx < len(row) else ""
            position_value = row[position_idx] if position_idx < len(row) else ""
            key = (normalize_id(chr_value), normalize_position(position_value))
            sprinzl_position = mapping.get(key)

            if sprinzl_position is None:
                sprinzl_position = position_value
            else:
                annotated_rows += 1

            base_row = [
                value
                for idx, value in enumerate(row)
                if idx not in existing_output_indices
            ]
            if len(base_row) < base_header_len:
                base_row.extend([""] * (base_header_len - len(base_row)))

            new_row = (
                base_row[:insert_idx]
                + [sprinzl_position]
                + base_row[insert_idx:]
            )
            writer.writerow(new_row)

    return {
        "input": str(input_file),
        "output": str(output_file),
        "rows": total_rows,
        "annotated": annotated_rows,
        "unannotated": total_rows - annotated_rows,
    }


def print_result(result: dict[str, object]) -> None:
    print(
        "[DONE] {input} -> {output} | rows={rows}, "
        "annotated={annotated}, unannotated={unannotated}".format(**result)
    )


def main() -> int:
    args = parse_args()

    try:
        mapping = load_mapping(
            map_file=Path(args.map_file).resolve(),
            map_sep=args.map_sep,
            map_id_col=args.map_id_col,
            map_position_col=args.map_position_col,
            map_sprinzl_col=args.map_sprinzl_col,
        )
        jobs, is_batch = resolve_jobs(args)
        config = build_config(args)

        print(f"[INFO] Files to process: {len(jobs)}")

        if not is_batch or args.workers == 1:
            for job in jobs:
                result = annotate_file(job[0], job[1], mapping, **config)
                print_result(result)
            return 0

        failed: list[tuple[str, str]] = []
        with ProcessPoolExecutor(
            max_workers=args.workers,
            initializer=_init_worker,
            initargs=(mapping, config),
        ) as executor:
            future_to_input = {
                executor.submit(_worker_annotate, job): str(job[0])
                for job in jobs
            }
            for future in as_completed(future_to_input):
                input_name = future_to_input[future]
                try:
                    result = future.result()
                except Exception as exc:
                    failed.append((input_name, str(exc)))
                    print(f"[ERROR] {input_name}: {exc}", file=sys.stderr)
                    if not args.continue_on_error:
                        for pending in future_to_input:
                            pending.cancel()
                        break
                else:
                    print_result(result)

        if failed:
            print("[ERROR] The following files failed:", file=sys.stderr)
            for input_name, message in failed:
                print(f"  - {input_name}: {message}", file=sys.stderr)
            return 1

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())

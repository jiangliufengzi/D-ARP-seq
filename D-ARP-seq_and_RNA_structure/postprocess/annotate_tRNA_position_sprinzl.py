#!/usr/bin/env python3
"""
================================================================================
脚本名称: annotate_tRNA_position_sprinzl.py
功能: 为 chr/position 位点表添加 Sprinzl 编号注释 (tRNA 标准位点命名)
================================================================================

背景/用途:
  本脚本用于将 tRNA 突变/修饰检测流水线中的原始坐标 (tRNAscan-SE ID +
  raw position)，映射为 Sprinzl 标准编号体系，便于跨研究比较和文献对照。

  典型输入来源:
    - tRNA 突变检测输出的 TSV 文件 (含 chr/trnascan_id 和 position/raw_pos 列)
    - 需要 Sprinzl 注释的位点汇总表

  输出用于下游分析:
    - 在 position 列后插入 sprinzl_position 列
    - 统计注释成功/未匹配位点数
    - 未匹配位点保留原始 position 值作为回退

  Sprinzl 编号体系:
    tRNA 学界通用的位点编号标准，以反密码子环为基准 (位置 34/35/36)，
    便于不同 tRNA 之间的位点对应与比较。

================================================================================

工作流程 (按执行顺序):

    ┌─────────────────────────────────────┐
    │  输入: -i 位点文件或目录            │
    │  映射表: hg38_trna_raw_to_sprinzl.tsv│
    └──────────────┬──────────────────────┘
                   │
                   ↓
    ┌─────────────────────────────────────┐
    │  main() / parse_args()              │
    │  ───────────────────────────────   │
    │  • 解析命令行参数                   │
    │  • 验证 in-place / output 互斥      │
    └──────────────┬──────────────────────┘
                   │
                   ↓
    ┌─────────────────────────────────────┐
    │  load_mapping()                     │
    │  加载 Sprinzl 映射表                │
    │  ───────────────────────────────   │
    │  • 读取 (trnascan_id, raw_pos) 键   │
    │  • 构建 → sprinzl_pos 字典          │
    │  • 校验重复键一致性                 │
    └──────────────┬──────────────────────┘
                   │
                   ↓
    ┌─────────────────────────────────────┐
    │  resolve_jobs()                     │
    │  确定输入/输出文件对                │
    │  ───────────────────────────────   │
    │  • 单文件 → 单 job                  │
    │  • 目录 → collect_input_files()     │
    │  • 支持 --in-place 原地覆盖         │
    └──────────────┬──────────────────────┘
                   │
                   ↓
    ┌─────────────────────────────────────┐
    │  annotate_file() × N                │
    │  [逐文件注释，可 ProcessPool 并行]  │
    │  ───────────────────────────────   │
    │  │  ┌────────────────────────────┐  │
    │  └→ │ write_annotated_file()    │  │
    │     │ [核心算法]                │  │
    │     │ • 读取表头，定位列索引    │  │
    │     │ • 在 position 后插入新列  │  │
    │     │ • 查映射表获取 Sprinzl 号 │  │
    │     │ • 未匹配则保留原 position │  │
    │     └───────────┬────────────────┘  │
    └─────────────────┼───────────────────┘
                      │
                      ↓
    ┌─────────────────────────────────────┐
    │  输出注释文件                        │
    │  ───────────────────────────────   │
    │  单文件: {stem}.sprinzl.tsv         │
    │  目录:   sprinzl_annotated/         │
    │  原地:   --in-place 覆盖原文件      │
    └─────────────────────────────────────┘

================================================================================
核心函数说明 (按执行顺序):

  主流程函数:
  1. main()                   - 主入口，加载映射表并分发任务
  2. parse_args()             - 解析命令行参数
  3. load_mapping()           - 加载 Sprinzl 映射表为字典
  4. resolve_jobs()           - 解析单文件/目录模式的输入输出对

  文件发现:
  5. collect_input_files()    - 目录模式按 glob pattern 收集文件
  6. insert_suffix()          - 生成带 .sprinzl 后缀的输出文件名

  注释核心:
  7. annotate_file()          - 单文件注释 (含 in-place 临时文件安全写入)
     └→ write_annotated_file() - [核心] 逐行查表并插入 sprinzl_position 列

  并行:
  8. _init_worker()           - ProcessPool worker 初始化 (共享映射表)
  9. _worker_annotate()       - 进程池 worker 封装

  辅助函数:
  - normalize_id()            - 标准化 ID/chr 值 (strip)
  - normalize_position()      - 标准化 position 值 (整数化)
  - open_text()               - 自动处理 .gz 压缩文件

================================================================================
核心算法说明 - write_annotated_file():

  映射查找逻辑:
  ─────────────────────────────────────────────────────────────────────
  步骤    操作                              说明
  ─────────────────────────────────────────────────────────────────────
  1       读取表头，定位 chr/position 列    列名可通过 --chr-col 等指定
  2       移除已有的 output_col 列          避免重复插入
  3       在 position 列后插入新列          sprinzl_position
  4       key = (normalize_id(chr),         标准化后查映射字典
             normalize_position(position))
  5       命中 → 写入 sprinzl_pos           annotated_rows += 1
  6       未命中 → 保留原 position 值       作为回退，不计入 annotated
  ─────────────────────────────────────────────────────────────────────

  映射表结构 (默认 hg38_trna_raw_to_sprinzl.tsv):
    - trnascan_id:  与输入 chr 列对应 (tRNAscan-SE ID)
    - raw_pos:      与输入 position 列对应 (原始坐标)
    - sprinzl_pos:  Sprinzl 标准编号 (输出值)

  in-place 安全写入:
    - 先写入临时文件 .{name}.sprinzl.tmp.{pid}
    - 成功后 os.replace() 原子替换原文件

================================================================================
输入格式:

  位点文件:
    - 带表头的 TSV/CSV (默认 tab 分隔)
    - 必需列: chr (或 --chr-col 指定) 和 position (或 --position-col 指定)
    - 支持 .gz 压缩
    - 单文件或目录批量 (*.tsv 默认)

  Sprinzl 映射表 (--map-file):
    - 默认: 脚本同目录下 hg38_trna_raw_to_sprinzl.tsv
    - 列: trnascan_id, raw_pos, sprinzl_pos (列名可自定义)
    - 键 (trnascan_id, raw_pos) 必须唯一且一致

================================================================================
输出格式:

  注释结果 (原文件 + sprinzl_position 列):
  ─────────────────────────────────────────────────────────────────────
  列名              说明
  ─────────────────────────────────────────────────────────────────────
  (原始列...)       输入文件所有原始列保留
  sprinzl_position  插入在 position 列之后
                    命中映射表 → Sprinzl 编号
                    未命中     → 保留原始 position 值
  ─────────────────────────────────────────────────────────────────────

  终端统计 (每文件):
    rows=总数据行数, annotated=成功映射行数, unannotated=未匹配行数

  输出路径规则:
    - 单文件无 -o: {stem}.sprinzl.tsv
    - 目录无 -o:   {input_dir}/sprinzl_annotated/
    - --in-place:   覆盖原文件

================================================================================
命令行参数:

  输入输出:
    -i, --input             输入文件或文件夹 (必需)
    -o, --output            输出文件或文件夹
    --in-place              原地覆盖 (不可与 -o 同用)
    --suffix                输出文件名后缀 (默认: .sprinzl)

  Sprinzl 映射表:
    --map-file              映射表路径
    --map-id-col            映射表 ID 列 (默认: trnascan_id)
    --map-position-col      映射表 position 列 (默认: raw_pos)
    --map-sprinzl-col       映射表 Sprinzl 列 (默认: sprinzl_pos)

  输入列名:
    --chr-col               输入 ID/chr 列 (默认: chr)
    --position-col          输入 position 列 (默认: position)
    --output-col            输出列名 (默认: sprinzl_position)

  文件格式:
    --sep                   输入/输出分隔符 (默认: tab)
    --map-sep               映射表分隔符 (默认: tab)

  批量模式:
    --pattern               glob 模式 (默认: *.tsv, 可重复)
    --recursive             递归搜索
    --workers               并行进程数 (默认: 1)
    --continue-on-error     单文件失败时继续处理

================================================================================
使用示例:

  1. 注释单个 TSV 文件:
      python /home/pf/14T/scripts/annotate_tRNA_position_sprinzl.py \\
         -i sample.tsv \\
         -o sample.sprinzl.tsv

  2. 目录批量 + 8 并行:
      python /home/pf/14T/scripts/annotate_tRNA_position_sprinzl.py \\
         -i input_dir \\
         -o output_dir \\
         --workers 8

  3. 自定义列名映射:
      python /home/pf/14T/scripts/annotate_tRNA_position_sprinzl.py \\
         -i sample.tsv \\
         --chr-col trnascan_id \\
         --position-col raw_pos

  4. 原地覆盖 + 递归搜索:
      python /home/pf/14T/scripts/annotate_tRNA_position_sprinzl.py \\
         -i mutation_results/ \\
         --in-place --recursive --pattern '*.tsv'

  5. 完整参数示例:
      python /home/pf/14T/scripts/annotate_tRNA_position_sprinzl.py \\
         -i /data/tRNA_mutations \\
         -o /data/tRNA_sprinzl \\
         --map-file /ref/hg38_trna_raw_to_sprinzl.tsv \\
         --chr-col chr --position-col position \\
         --workers 4 --continue-on-error

================================================================================
依赖:
  - Python 3.10+ (使用 type hints 与 union 语法)
  - 标准库: argparse, csv, gzip, re, concurrent.futures, pathlib

================================================================================
注意事项:
  1. 映射表键 (trnascan_id, raw_pos) 必须与输入文件的 chr/position 精确匹配
  2. 未匹配位点不会报错，而是保留原始 position 值 (unannotated 计数增加)
  3. --in-place 使用临时文件 + 原子替换，避免写入中断损坏原文件
  4. 目录模式默认排除映射表自身，避免误处理
  5. 映射表中同一键对应不同 Sprinzl 值时会报错终止

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
        description="给含 chr 和 position 列的文件添加 Sprinzl position 注释。",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    input_group = parser.add_argument_group("输入输出")
    input_group.add_argument(
        "-i",
        "--input",
        required=True,
        help="输入文件或输入文件夹。",
    )
    input_group.add_argument(
        "-o",
        "--output",
        default=None,
        help=(
            "输出文件或输出文件夹。单文件输入时可省略，默认在原文件名中插入后缀；"
            "文件夹输入时可省略，默认写入输入文件夹下的 sprinzl_annotated/。"
        ),
    )
    input_group.add_argument(
        "--in-place",
        action="store_true",
        help="原地覆盖输入文件。启用后不能同时指定 --output。",
    )
    input_group.add_argument(
        "--suffix",
        default=".sprinzl",
        help="未原地覆盖时，插入到输出文件名中的后缀（默认: .sprinzl）。",
    )

    map_group = parser.add_argument_group("Sprinzl 映射表")
    map_group.add_argument(
        "--map-file",
        default=str(DEFAULT_MAP_FILE),
        help=f"Sprinzl 映射表路径（默认: {DEFAULT_MAP_FILE}）。",
    )
    map_group.add_argument(
        "--map-id-col",
        default="trnascan_id",
        help="映射表中与输入 chr 列对应的列名（默认: trnascan_id）。",
    )
    map_group.add_argument(
        "--map-position-col",
        default="raw_pos",
        help="映射表中与输入 position 列对应的列名（默认: raw_pos）。",
    )
    map_group.add_argument(
        "--map-sprinzl-col",
        default="sprinzl_pos",
        help="映射表中要写入输出的 Sprinzl 列名（默认: sprinzl_pos）。",
    )

    column_group = parser.add_argument_group("输入列名")
    column_group.add_argument(
        "--chr-col",
        default="chr",
        help="输入文件中的 ID/chr 列名（默认: chr）。",
    )
    column_group.add_argument(
        "--position-col",
        default="position",
        help="输入文件中的 position 列名（默认: position）。",
    )
    column_group.add_argument(
        "--output-col",
        default="sprinzl_position",
        help="新增/替换的输出列名（默认: sprinzl_position）。",
    )

    format_group = parser.add_argument_group("文件格式")
    format_group.add_argument(
        "--sep",
        default="tab",
        help="输入/输出文件分隔符: tab, comma, semicolon, pipe 或单个字符（默认: tab）。",
    )
    format_group.add_argument(
        "--map-sep",
        default="tab",
        help="映射表分隔符: tab, comma, semicolon, pipe 或单个字符（默认: tab）。",
    )

    batch_group = parser.add_argument_group("文件夹批量模式")
    batch_group.add_argument(
        "--pattern",
        action="append",
        default=None,
        help="文件夹模式下要处理的 glob pattern，可重复指定（默认: *.tsv）。",
    )
    batch_group.add_argument(
        "--recursive",
        action="store_true",
        help="文件夹模式下递归搜索输入文件。",
    )
    batch_group.add_argument(
        "--workers",
        type=int,
        default=1,
        help="文件夹模式并行进程数（默认: 1）。",
    )
    batch_group.add_argument(
        "--continue-on-error",
        action="store_true",
        help="文件夹模式下单个文件失败时继续处理其它文件，最后汇总失败信息。",
    )

    args = parser.parse_args()

    if args.in_place and args.output:
        parser.error("--in-place 不能和 --output 同时使用。")
    if args.workers < 1:
        parser.error("--workers 必须 >= 1。")

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
        raise SystemExit(f"[ERROR] {arg_name} 只能是一个字符或 tab/comma/semicolon/pipe。")
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
        raise FileNotFoundError(f"映射表不存在: {map_file}")

    mapping: dict[tuple[str, str], str] = {}
    row_count = 0
    duplicate_count = 0

    with open_text(map_file, "rt") as handle:
        reader = csv.DictReader(handle, delimiter=map_sep)
        if reader.fieldnames is None:
            raise ValueError(f"映射表为空: {map_file}")

        required_cols = [map_id_col, map_position_col, map_sprinzl_col]
        missing_cols = [col for col in required_cols if col not in reader.fieldnames]
        if missing_cols:
            raise ValueError(
                f"映射表缺少列: {', '.join(missing_cols)}; "
                f"已有列: {', '.join(reader.fieldnames)}"
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
                        "映射表中同一个 ID/position 对应多个 Sprinzl 值: "
                        f"{key[0]} {key[1]} -> {mapping[key]} / {value}"
                    )
                duplicate_count += 1
                continue
            mapping[key] = value

    if not mapping:
        raise ValueError(f"映射表没有可用记录: {map_file}")

    print(
        f"[INFO] 读取 Sprinzl 映射: {len(mapping)} 个键 "
        f"({row_count} 行, {duplicate_count} 个重复一致键)"
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
        raise FileNotFoundError(f"输入路径不存在: {input_path}")

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
        raise ValueError(f"输入路径既不是文件也不是文件夹: {input_path}")

    input_files = collect_input_files(
        input_path,
        args.pattern,
        args.recursive,
        Path(args.map_file).resolve(),
    )
    if not input_files:
        raise ValueError(
            f"未在 {input_path} 中找到匹配文件: {', '.join(args.pattern)}"
        )

    if args.in_place:
        return [(path, path) for path in input_files], True

    output_dir = (
        Path(args.output).resolve()
        if args.output
        else input_path / "sprinzl_annotated"
    )
    if output_dir.exists() and not output_dir.is_dir():
        raise ValueError(f"文件夹输入时 --output 必须是输出文件夹: {output_dir}")

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
        raise RuntimeError("worker 未初始化")
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
            raise ValueError(f"输入文件为空: {input_file}")

        if chr_col not in header:
            raise ValueError(f"{input_file} 缺少列: {chr_col}; 已有列: {', '.join(header)}")
        if position_col not in header:
            raise ValueError(
                f"{input_file} 缺少列: {position_col}; 已有列: {', '.join(header)}"
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

        print(f"[INFO] 待处理文件数: {len(jobs)}")

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
            print("[ERROR] 以下文件处理失败:", file=sys.stderr)
            for input_name, message in failed:
                print(f"  - {input_name}: {message}", file=sys.stderr)
            return 1

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Compare reproduced folder-14 CSV tables with reference tables."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Dict, List

import pandas as pd


def collect_csvs(roots: List[Path]) -> Dict[str, Path]:
    root_label = ", ".join(str(root) for root in roots)
    files = []
    for root in roots:
        files.extend(sorted(root.glob("*.csv")))
    result: Dict[str, Path] = {}
    for path in files:
        name = path.name
        if name in result:
            raise ValueError(
                f"Duplicate CSV basename below {root_label}: {name}\n"
                f"  {result[name]}\n  {path}"
            )
        result[name] = path
    return result


def normalize_rows(frame: pd.DataFrame) -> pd.DataFrame:
    keys = [column for column in ("chr", "position") if column in frame.columns]
    if keys:
        return frame.sort_values(keys, kind="stable").reset_index(drop=True)
    return frame.reset_index(drop=True)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare reference and reproduced 14_merge_DHU_site_rep CSVs."
    )
    parser.add_argument(
        "--reference", type=Path, action="append", required=True,
        help="Reference folder 14; repeat this option for multiple folders.",
    )
    parser.add_argument("--reproduced", type=Path, required=True)
    parser.add_argument("--atol", type=float, default=1e-12)
    parser.add_argument("--rtol", type=float, default=1e-10)
    args = parser.parse_args()

    reference = collect_csvs(args.reference)
    reproduced = collect_csvs([args.reproduced])
    missing = sorted(set(reference) - set(reproduced))
    extra = sorted(set(reproduced) - set(reference))
    failures: List[str] = []

    if missing:
        failures.append("Missing files: " + ", ".join(missing))
    if extra:
        failures.append("Extra files: " + ", ".join(extra))

    for name in sorted(set(reference) & set(reproduced)):
        expected = normalize_rows(pd.read_csv(reference[name]))
        observed = normalize_rows(pd.read_csv(reproduced[name]))
        try:
            pd.testing.assert_frame_equal(
                observed,
                expected,
                check_dtype=False,
                check_exact=False,
                atol=args.atol,
                rtol=args.rtol,
                check_like=False,
            )
            print(f"OK   {name} ({len(expected)} rows)")
        except AssertionError as error:
            failures.append(f"DIFF {name}: {str(error).splitlines()[0]}")

    if failures:
        print("\nComparison failed:", file=sys.stderr)
        for failure in failures:
            print(f"- {failure}", file=sys.stderr)
        return 1

    print(f"\nAll {len(reference)} CSV files match.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

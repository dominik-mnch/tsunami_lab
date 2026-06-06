#!/usr/bin/env python3
"""Clean generated run directories below solutions/."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Remove generated solutions/run_* directories."
    )
    parser.add_argument(
        "--solutions-dir",
        default="solutions",
        type=Path,
        help="Path to the solutions directory. Defaults to ./solutions.",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Also remove solutions/solution.nc and solutions/checkpoint.nc.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be removed without deleting anything.",
    )
    parser.add_argument(
        "--yes",
        "-y",
        action="store_true",
        help="Skip the confirmation prompt.",
    )
    return parser.parse_args()


def collect_targets(solutions_dir: Path, include_main_files: bool) -> list[Path]:
    if not solutions_dir.exists():
        return []

    targets = sorted(path for path in solutions_dir.glob("run_*") if path.is_dir())

    if include_main_files:
        for file_name in ("solution.nc", "checkpoint.nc"):
            path = solutions_dir / file_name
            if path.exists():
                targets.append(path)

    return targets


def remove_target(path: Path) -> None:
    if path.is_dir():
        shutil.rmtree(path)
    else:
        path.unlink()


def main() -> int:
    args = parse_args()
    solutions_dir = args.solutions_dir.resolve()
    targets = collect_targets(solutions_dir, args.all)

    if not targets:
        print(f"No generated solution runs found in {solutions_dir}.")
        return 0

    action = "Would remove" if args.dry_run else "Will remove"
    print(f"{action} {len(targets)} item(s):")
    for target in targets:
        print(f"  {target}")

    if args.dry_run:
        return 0

    if not args.yes:
        answer = input("Continue? [y/N] ").strip().lower()
        if answer not in {"y", "yes"}:
            print("Aborted.")
            return 1

    for target in targets:
        remove_target(target)

    print("Cleanup complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
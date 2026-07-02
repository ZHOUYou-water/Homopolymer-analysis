#!/usr/bin/env python3
"""Filter ONT FASTQ by read length and mean Phred quality."""

from __future__ import annotations

import argparse
from pathlib import Path
from statistics import mean


def mean_phred(qual: str) -> float:
    return mean(ord(c) - 33 for c in qual)


def filter_fastq(
    in_fastq: Path,
    out_fastq: Path,
    min_length: int,
    min_mean_q: float,
) -> tuple[int, int]:
    total = kept = 0
    with in_fastq.open() as src, out_fastq.open("w") as dst:
        while True:
            header = src.readline()
            if not header:
                break
            seq = src.readline().strip()
            plus = src.readline()
            qual = src.readline().strip()
            total += 1

            if len(seq) < min_length:
                continue
            if mean_phred(qual) < min_mean_q:
                continue

            kept += 1
            dst.write(header)
            dst.write(seq + "\n")
            dst.write(plus)
            dst.write(qual + "\n")

            if total % 100000 == 0:
                print(f"scanned {total:,} kept {kept:,}", flush=True)

    return total, kept


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-i", "--input", required=True, type=Path, help="Input FASTQ")
    ap.add_argument("-o", "--output", required=True, type=Path, help="Output FASTQ")
    ap.add_argument("-l", "--min-length", type=int, default=1000, help="Min read length (default: 1000)")
    ap.add_argument("-q", "--min-mean-q", type=float, default=7.0, help="Min mean Q (default: 7)")
    args = ap.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    total, kept = filter_fastq(args.input, args.output, args.min_length, args.min_mean_q)

    stats = args.output.with_suffix(args.output.suffix + ".stats.txt")
    stats.write_text(
        "\n".join(
            [
                f"input={args.input}",
                f"output={args.output}",
                f"min_length={args.min_length}",
                f"min_mean_q={args.min_mean_q}",
                f"reads_total={total}",
                f"reads_kept={kept}",
                f"reads_removed={total - kept}",
            ]
        )
        + "\n"
    )
    print(f"Done: kept {kept:,}/{total:,} reads -> {args.output}")


if __name__ == "__main__":
    main()

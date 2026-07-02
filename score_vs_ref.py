#!/usr/bin/env python3
"""Score assemblies against a reference with MUMmer (indel count)."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


def read_len(path: Path) -> int:
    return len("".join(l.strip() for l in path.read_text().splitlines() if not l.startswith(">")))


def count_variants(mummer: str, ref: str, asm: str, prefix: Path) -> dict:
    subprocess.run(
        [f"{mummer}/nucmer", "--maxmatch", "-p", str(prefix), ref, asm],
        check=True,
        capture_output=True,
    )
    out = subprocess.check_output(
        [f"{mummer}/show-snps", "-T", f"{prefix}.delta"], text=True
    )
    indels = snps = ins = dels = 0
    for line in out.splitlines():
        if line.startswith("NUCMER") or line.startswith("["):
            continue
        p = line.split("\t")
        if len(p) < 4:
            continue
        if p[1] == "." or p[2] == ".":
            indels += 1
            if p[1] == ".":
                ins += 1
            else:
                dels += 1
        else:
            snps += 1
    return {
        "indels": indels,
        "insertions": ins,
        "deletions": dels,
        "snps": snps,
        "total": indels + snps,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ref", required=True, help="Reference FASTA")
    ap.add_argument("--assembly", required=True, help="Assembly FASTA")
    ap.add_argument("--mummer", default="mummer", help="MUMmer install prefix or PATH")
    ap.add_argument("--prefix", default="asm_vs_ref", help="Output prefix for delta files")
    args = ap.parse_args()

    ref = Path(args.ref)
    asm = Path(args.assembly)
    prefix = Path(args.prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    ref_len = read_len(ref)
    asm_len = read_len(asm)
    v = count_variants(args.mummer, str(ref), str(asm), prefix)

    tsv = prefix.with_suffix(".scores.tsv")
    tsv.write_text(
        "assembly\tlength\tlength_diff_vs_ref\tindels\tinsertions\tdeletions\tsnps\ttotal_variants\n"
        f"{asm}\t{asm_len}\t{asm_len - ref_len}\t{v['indels']}\t{v['insertions']}\t"
        f"{v['deletions']}\t{v['snps']}\t{v['total']}\n"
    )
    print(tsv.read_text())


if __name__ == "__main__":
    main()

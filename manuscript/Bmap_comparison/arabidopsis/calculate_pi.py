#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import subprocess
from pathlib import Path

import numpy as np

from Bvalcalc.utils.load_chr_sizes import load_chr_sizes


SCALES = (1_000, 10_000, 100_000, 1_000_000)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bcftools", required=True, type=Path)
    parser.add_argument("--vcf", required=True, type=Path)
    parser.add_argument("--chr-sizes", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    sizes = load_chr_sizes(args.chr_sizes)
    contigs = list(sizes)
    accumulators = {
        scale: {
            chrom: (np.zeros(size // scale, dtype=float), np.zeros(size // scale, dtype=np.int64))
            for chrom, size in sizes.items()
        }
        for scale in SCALES
    }

    view = subprocess.Popen(
        [str(args.bcftools), "view", "--threads", str(args.threads), "-r", ",".join(contigs),
         "-m2", "-M2", "-Ou", str(args.vcf)],
        stdout=subprocess.PIPE,
    )
    fill = subprocess.Popen(
        [str(args.bcftools), "+fill-tags", "-Ou", "--", "-t", "AC,AN"],
        stdin=view.stdout,
        stdout=subprocess.PIPE,
    )
    if view.stdout is not None:
        view.stdout.close()
    query = subprocess.Popen(
        [str(args.bcftools), "query", "-f", "%CHROM\\t%POS\\t%INFO/AC\\t%INFO/AN\\n"],
        stdin=fill.stdout,
        stdout=subprocess.PIPE,
        text=True,
    )
    if fill.stdout is not None:
        fill.stdout.close()
    if query.stdout is None:
        raise RuntimeError("Failed to open bcftools query output")

    records = polymorphic = 0
    for line in query.stdout:
        records += 1
        chrom, pos_text, ac_text, an_text = line.rstrip().split("\t")
        if chrom not in sizes or "," in ac_text:
            continue
        pos, ac, an = int(pos_text), int(ac_text), int(an_text)
        if an < 2 or ac <= 0 or ac >= an:
            continue
        polymorphic += 1
        mpd = (2.0 * ac * (an - ac)) / (an * (an - 1))
        for scale in SCALES:
            index = (pos - 1) // scale
            pi_num, n_vars = accumulators[scale][chrom]
            if index < len(pi_num):
                pi_num[index] += mpd
                n_vars[index] += 1

    query_rc = query.wait()
    fill_rc = fill.wait()
    view_rc = view.wait()
    if (query_rc, fill_rc, view_rc) != (0, 0, 0):
        raise SystemExit(f"bcftools pipeline failed: query={query_rc}, fill={fill_rc}, view={view_rc}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    for scale in SCALES:
        path = args.out_dir / f"pi_south_{scale // 1000}kb.csv.gz"
        with gzip.open(path, "wt", encoding="ascii", newline="") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(("chrom", "start", "end", "center", "pi_num", "n_vars_used",
                             "n_neutral_sites", "pi", "pi_window_bp", "pi_step_bp"))
            for chrom in contigs:
                pi_num, n_vars = accumulators[scale][chrom]
                for index in range(len(pi_num)):
                    start = index * scale + 1
                    end = (index + 1) * scale
                    neutral_sites = scale - 1  # Matches the paper plotting script's [start,end) denominator.
                    pi = pi_num[index] / neutral_sites if n_vars[index] >= 1 else np.nan
                    writer.writerow((f"chr{chrom}", start, end, start + scale / 2 - 0.5,
                                     f"{pi_num[index]:.17g}", int(n_vars[index]), neutral_sites,
                                     f"{pi:.17g}" if np.isfinite(pi) else "", scale, scale))
        print(f"wrote\t{path}")
    print(f"records_biallelic\t{records}")
    print(f"records_polymorphic\t{polymorphic}")


if __name__ == "__main__":
    main()


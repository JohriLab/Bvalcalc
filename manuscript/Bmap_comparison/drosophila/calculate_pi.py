#!/usr/bin/env python3
"""Calculate neutral-site D. melanogaster pi windows used for B-map tests.

Neutral sites are unannotated bases plus fourfold-degenerate CDS bases. UTR,
phastCons, low-confidence, and inaccessible bases are excluded. The VCF INFO
field ``Degeneracy`` supplies fourfold-degenerate positions.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path

import allel
import numpy as np
import pandas as pd


COLUMNS = [
    "chrom", "start", "end", "center", "pi_num", "n_vars_used",
    "n_neutral_sites", "n_syn_sites", "n_intergenic_sites",
    "n_syn_snps_used", "n_intergenic_snps_used", "pi",
    "pi_window_bp", "pi_step_bp",
]


def chrom_name(value: str) -> str:
    value = str(value)
    return value if value.startswith("chr") else f"chr{value}"


def read_sizes(path: Path) -> dict[str, int]:
    frame = pd.read_csv(path, sep=r"[\s,]+", engine="python", header=None,
                        usecols=[0, 1], names=["chrom", "length"])
    return dict(zip(frame["chrom"].map(chrom_name), frame["length"].astype(int)))


def read_bed(path: Path | None) -> pd.DataFrame:
    if path is None or str(path) == "/dev/null":
        return pd.DataFrame(columns=["chrom", "start", "end"])
    frame = pd.read_csv(path, sep=r"\s+", header=None, usecols=[0, 1, 2],
                        names=["chrom", "start", "end"])
    frame["chrom"] = frame["chrom"].map(chrom_name)
    frame[["start", "end"]] = frame[["start", "end"]].astype(int)
    return frame[frame["end"] > frame["start"]].sort_values(
        ["chrom", "start", "end"]
    )


def merge(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    result = [sorted(intervals)[0]]
    for start, end in sorted(intervals)[1:]:
        old_start, old_end = result[-1]
        if start <= old_end:
            result[-1] = (old_start, max(old_end, end))
        else:
            result.append((start, end))
    return result


def within(frame: pd.DataFrame, start: int, end: int) -> list[tuple[int, int]]:
    overlap = frame[(frame["end"] > start) & (frame["start"] < end)]
    return [(max(int(row.start), start), min(int(row.end), end))
            for row in overlap.itertuples()]


def intersect(left, right) -> list[tuple[int, int]]:
    result = []
    i = j = 0
    while i < len(left) and j < len(right):
        start, end = max(left[i][0], right[j][0]), min(left[i][1], right[j][1])
        if end > start:
            result.append((start, end))
        if left[i][1] < right[j][1]:
            i += 1
        else:
            j += 1
    return result


def total_length(intervals) -> int:
    return sum(end - start for start, end in intervals)


def positions_in(positions: np.ndarray, intervals) -> np.ndarray:
    result = np.zeros(len(positions), dtype=bool)
    for start, end in intervals:
        result |= (positions >= start) & (positions < end)
    return result


def fourfold_positions(vcf: Path, chrom: str) -> np.ndarray:
    if shutil.which("bcftools") is None:
        raise SystemExit("bcftools is required to query INFO/Degeneracy")
    command = ["bcftools", "query", "-r", chrom.removeprefix("chr"),
               "-f", "%POS\\t%INFO/Degeneracy\\n", str(vcf)]
    output = subprocess.run(command, check=True, capture_output=True, text=True).stdout
    positions = []
    for line in output.splitlines():
        position, value = line.split("\t", 1)
        tokens = value.replace(" ", "").split(",")
        is_fourfold = False
        for token in tokens:
            try:
                is_fourfold |= float(token) == 4.0
            except ValueError:
                pass
        if is_fourfold:
            positions.append(int(position))
    return np.asarray(sorted(positions), dtype=np.int64)


def count_positions(positions: np.ndarray, intervals) -> int:
    return int(sum(
        np.searchsorted(positions, end, side="left")
        - np.searchsorted(positions, start, side="left")
        for start, end in intervals
    ))


def members(reference: np.ndarray, query: np.ndarray) -> np.ndarray:
    indices = np.searchsorted(reference, query, side="left")
    result = np.zeros(len(query), dtype=bool)
    valid = indices < len(reference)
    result[valid] = reference[indices[valid]] == query[valid]
    return result


def calculate_chromosome(vcf, chrom, length, beds, fourfold, window_bp, samples):
    rows = []
    beds = {name: frame[frame["chrom"] == chrom] for name, frame in beds.items()}
    for zero_start in range(0, max(0, length - window_bp + 1), window_bp):
        # Preserve the coordinates and denominator used in the original workflow.
        start, end = zero_start + 1, zero_start + window_bp
        intervals = {name: merge(within(frame, start, end))
                     for name, frame in beds.items()}
        selected = merge(intervals["utr"] + intervals["phastcons"]
                         + intervals["lowconf"] + intervals["cds"]
                         + intervals["inaccessible"])
        other_masks = merge(intervals["utr"] + intervals["phastcons"]
                            + intervals["lowconf"] + intervals["inaccessible"])
        n_syn = count_positions(fourfold, intervals["cds"])
        n_syn -= count_positions(fourfold, intersect(intervals["cds"], other_masks))
        n_neutral = end - start - total_length(selected) + n_syn

        callset = allel.read_vcf(
            str(vcf), region=f"{chrom.removeprefix('chr')}:{start}-{end}",
            samples=samples, fields=["variants/POS", "calldata/GT"],
        )
        pi_num = np.nan
        n_vars = n_syn_snps = n_intergenic_snps = 0
        if n_neutral > 0 and callset is not None and callset.get("variants/POS") is not None:
            positions = np.asarray(callset["variants/POS"], dtype=np.int64)
            counts = allel.GenotypeArray(callset["calldata/GT"]).count_alleles()
            if counts.shape[1] >= 2:
                ac0, ac1 = counts[:, 0].astype(int), counts[:, 1].astype(int)
                n_hap = ac0 + ac1
                multiallelic = (np.any(counts[:, 2:] > 0, axis=1)
                                if counts.shape[1] > 2 else np.zeros(len(positions), bool))
                masks = {name: positions_in(positions, value)
                         for name, value in intervals.items()}
                is_fourfold = members(fourfold, positions)
                with np.errstate(invalid="ignore", divide="ignore"):
                    mpd = 2.0 * ac1 * (n_hap - ac1) / (n_hap * (n_hap - 1))
                valid = (~(masks["utr"] | masks["phastcons"] | masks["lowconf"]
                           | masks["inaccessible"])
                         & (~masks["cds"] | is_fourfold)
                         & (ac0 > 0) & (ac1 > 0) & ~multiallelic
                         & (n_hap >= 2) & np.isfinite(mpd))
                n_vars = int(valid.sum())
                if n_vars:
                    pi_num = float(mpd[valid].sum())
                    n_syn_snps = int(np.sum(valid & masks["cds"] & is_fourfold))
                    n_intergenic_snps = int(np.sum(valid & ~masks["cds"]))

        pi = pi_num / n_neutral if n_neutral > 0 and np.isfinite(pi_num) else np.nan
        rows.append((chrom, start, end, start + window_bp / 2 - 0.5,
                     pi_num, n_vars, n_neutral, n_syn, n_neutral - n_syn,
                     n_syn_snps, n_intergenic_snps, pi, window_bp, window_bp))
    return pd.DataFrame(rows, columns=COLUMNS)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True, type=Path)
    parser.add_argument("--chrom-sizes", required=True, type=Path)
    parser.add_argument("--cds-bed", required=True, type=Path)
    parser.add_argument("--utr-bed", required=True, type=Path)
    parser.add_argument("--phastcons-bed", required=True, type=Path)
    parser.add_argument("--lowconf-bed", type=Path)
    parser.add_argument("--inaccessible-bed", type=Path)
    parser.add_argument("--contigs", default="2L,2R,3L,3R")
    parser.add_argument("--samples", type=Path, help="Optional one-sample-per-line file")
    parser.add_argument("--window-bp", required=True, type=int)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    sizes = read_sizes(args.chrom_sizes)
    contigs = [chrom_name(value.strip()) for value in args.contigs.split(",")]
    samples = ([line.strip() for line in args.samples.read_text().splitlines()
                if line.strip()] if args.samples else None)
    beds = {"cds": read_bed(args.cds_bed), "utr": read_bed(args.utr_bed),
            "phastcons": read_bed(args.phastcons_bed),
            "lowconf": read_bed(args.lowconf_bed),
            "inaccessible": read_bed(args.inaccessible_bed)}
    results = []
    for chrom in contigs:
        print(f"Processing {chrom}", flush=True)
        results.append(calculate_chromosome(
            args.vcf, chrom, sizes[chrom], beds, fourfold_positions(args.vcf, chrom),
            args.window_bp, samples,
        ))
    output = pd.concat(results, ignore_index=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.out, index=False,
                  compression="gzip" if args.out.suffix == ".gz" else None)
    print(f"Wrote {args.out} ({len(output):,} windows)")


if __name__ == "__main__":
    main()

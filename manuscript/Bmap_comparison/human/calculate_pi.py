#!/usr/bin/env python3
"""Build YRI pi windows using an explicit accessible-interval denominator."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import allel
import numpy as np
import pandas as pd


OUTPUT_COLUMNS = [
    "chrom",
    "start",
    "end",
    "center",
    "pi_num",
    "n_vars_used",
    "n_neutral_sites",
    "n_syn_sites",
    "n_intergenic_sites",
    "n_syn_snps_used",
    "n_intergenic_snps_used",
    "pi",
    "pi_window_bp",
    "pi_step_bp",
]


def parse_contigs(raw: str) -> list[str]:
    return [item.strip().removeprefix("chr") for item in raw.split(",") if item.strip()]


def read_chrom_sizes(path: str) -> dict[str, int]:
    frame = pd.read_csv(path, header=None, names=["chrom", "length"])
    if frame.shape[1] == 1:
        frame = pd.read_csv(path, header=None, names=["record"])
        split = frame["record"].str.split(",", n=1, expand=True)
        frame = pd.DataFrame({"chrom": split[0], "length": split[1]})
    frame["chrom"] = frame["chrom"].astype(str).str.removeprefix("chr")
    frame["length"] = pd.to_numeric(frame["length"], errors="raise").astype(np.int64)
    return dict(zip(frame["chrom"], frame["length"]))


def merge_intervals(starts: np.ndarray, ends: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    order = np.argsort(starts, kind="stable")
    starts = starts[order].astype(np.int64, copy=False)
    ends = ends[order].astype(np.int64, copy=False)
    valid = ends > starts
    starts = starts[valid]
    ends = ends[valid]
    if len(starts) < 2 or np.all(starts[1:] >= ends[:-1]):
        return starts, ends

    merged_starts = [int(starts[0])]
    merged_ends = [int(ends[0])]
    for start, end in zip(starts[1:], ends[1:]):
        if int(start) <= merged_ends[-1]:
            merged_ends[-1] = max(merged_ends[-1], int(end))
        else:
            merged_starts.append(int(start))
            merged_ends.append(int(end))
    return np.asarray(merged_starts), np.asarray(merged_ends)


def read_accessible_intervals(
    path: str, contigs: list[str]
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    frame = pd.read_csv(
        path,
        sep="\t",
        header=None,
        usecols=[0, 1, 2],
        names=["chrom", "start", "end"],
        dtype={"chrom": str, "start": np.int64, "end": np.int64},
    )
    frame["chrom"] = frame["chrom"].str.removeprefix("chr")
    frame = frame[frame["chrom"].isin(contigs)]
    result = {}
    for chrom, group in frame.groupby("chrom", sort=False):
        result[chrom] = merge_intervals(
            group["start"].to_numpy(), group["end"].to_numpy()
        )
    missing = sorted(set(contigs) - set(result))
    if missing:
        raise ValueError(f"Accessible BED lacks requested contigs: {missing}")
    return result


def accessible_bp_per_window(
    starts: np.ndarray,
    ends: np.ndarray,
    chrom_length: int,
    window_bp: int,
) -> np.ndarray:
    n_windows = chrom_length // window_bp
    covered_end = n_windows * window_bp
    starts = np.clip(starts, 0, covered_end)
    ends = np.clip(ends, 0, covered_end)
    valid = ends > starts
    starts = starts[valid]
    ends = ends[valid]

    denominator = np.zeros(n_windows, dtype=np.int64)
    first = starts // window_bp
    last = (ends - 1) // window_bp
    same = first == last
    np.add.at(denominator, first[same], ends[same] - starts[same])

    for start, end, first_bin, last_bin in zip(
        starts[~same], ends[~same], first[~same], last[~same]
    ):
        denominator[first_bin] += (first_bin + 1) * window_bp - start
        denominator[last_bin] += end - last_bin * window_bp
        if last_bin > first_bin + 1:
            denominator[first_bin + 1 : last_bin] += window_bp
    return denominator


def positions_in_intervals(
    positions0: np.ndarray, starts: np.ndarray, ends: np.ndarray
) -> np.ndarray:
    interval_index = np.searchsorted(starts, positions0, side="right") - 1
    valid = interval_index >= 0
    valid_indices = interval_index[valid]
    valid[valid] = positions0[valid] < ends[valid_indices]
    return valid


def process_chromosome(task: tuple) -> pd.DataFrame:
    (
        vcf_path,
        chrom,
        chrom_length,
        starts,
        ends,
        window_bp,
    ) = task
    n_windows = chrom_length // window_bp
    denominator = accessible_bp_per_window(
        starts, ends, chrom_length=chrom_length, window_bp=window_bp
    )
    callset = allel.read_vcf(
        vcf_path,
        region=chrom,
        fields=["variants/POS", "calldata/GT"],
    )
    pi_num = np.zeros(n_windows, dtype=float)
    n_vars = np.zeros(n_windows, dtype=np.int64)

    if callset is not None and callset.get("variants/POS") is not None:
        positions = np.asarray(callset["variants/POS"], dtype=np.int64)
        genotype = allel.GenotypeArray(callset["calldata/GT"])
        allele_counts = genotype.count_alleles()
        accessible = positions_in_intervals(positions - 1, starts, ends)

        if allele_counts.shape[1] > 2:
            multiallelic = np.any(allele_counts[:, 2:] > 0, axis=1)
        else:
            multiallelic = np.zeros(len(positions), dtype=bool)
        ac0 = allele_counts[:, 0].astype(np.int64, copy=False)
        ac1 = allele_counts[:, 1].astype(np.int64, copy=False)
        n_haplotypes = ac0 + ac1
        polymorphic = (ac0 > 0) & (ac1 > 0)
        with np.errstate(invalid="ignore", divide="ignore"):
            mpd = (
                2.0
                * ac1
                * (n_haplotypes - ac1)
                / (n_haplotypes * (n_haplotypes - 1))
            )
        valid = (
            accessible
            & ~multiallelic
            & polymorphic
            & (n_haplotypes >= 2)
            & np.isfinite(mpd)
        )
        window_index = (positions[valid] - 1) // window_bp
        in_complete_windows = window_index < n_windows
        window_index = window_index[in_complete_windows]
        mpd = mpd[valid][in_complete_windows]
        pi_num = np.bincount(
            window_index, weights=mpd, minlength=n_windows
        ).astype(float)
        n_vars = np.bincount(window_index, minlength=n_windows).astype(np.int64)

    start1 = np.arange(n_windows, dtype=np.int64) * window_bp + 1
    end = start1 + window_bp - 1
    with np.errstate(invalid="ignore", divide="ignore"):
        pi = pi_num / denominator
    pi[denominator == 0] = np.nan

    return pd.DataFrame(
        {
            "chrom": f"chr{chrom}",
            "start": start1,
            "end": end,
            "center": start1 + window_bp / 2 - 0.5,
            "pi_num": pi_num,
            "n_vars_used": n_vars,
            "n_neutral_sites": denominator,
            "n_syn_sites": 0,
            "n_intergenic_sites": denominator,
            "n_syn_snps_used": 0,
            "n_intergenic_snps_used": n_vars,
            "pi": pi,
            "pi_window_bp": window_bp,
            "pi_step_bp": window_bp,
        }
    )[OUTPUT_COLUMNS]


def validate(frame: pd.DataFrame, window_bp: int) -> dict[str, float | int]:
    denominator = frame["n_neutral_sites"].to_numpy(dtype=np.int64)
    positive = denominator[denominator > 0]
    minimum_positive = max(20, int(len(frame) * 0.6))
    if len(positive) < minimum_positive:
        raise ValueError(f"Only {len(positive)} windows have accessible sequence")
    minimum_unique = min(100, max(10, int(len(positive) * 0.5)))
    if np.unique(positive).size < minimum_unique:
        raise ValueError("Accessible denominators do not vary enough; mask likely failed")
    if np.all(positive == window_bp - 1) or np.all(positive == window_bp):
        raise ValueError("Every denominator equals window size; mask was not applied")
    median = float(np.median(positive))
    if not 0.10 * window_bp <= median <= 0.95 * window_bp:
        raise ValueError(f"Implausible median accessible denominator: {median}")
    return {
        "n_windows": len(frame),
        "n_positive_accessible_windows": len(positive),
        "n_finite_pi_windows": int(np.isfinite(frame["pi"]).sum()),
        "n_positive_pi_windows": int((frame["pi"] > 0).sum()),
        "median_accessible_bp": median,
        "mean_accessible_bp": float(np.mean(positive)),
        "min_accessible_bp": int(np.min(positive)),
        "max_accessible_bp": int(np.max(positive)),
        "unique_accessible_denominators": int(np.unique(positive).size),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--accessible-bed", required=True)
    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--contigs", default="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22")
    parser.add_argument("--window-bp", type=int, default=1_000_000)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--out", required=True)
    parser.add_argument("--validation-out", required=True)
    args = parser.parse_args()

    contigs = parse_contigs(args.contigs)
    sizes = read_chrom_sizes(args.chrom_sizes)
    accessible = read_accessible_intervals(args.accessible_bed, contigs)
    tasks = [
        (
            args.vcf,
            chrom,
            sizes[chrom],
            accessible[chrom][0],
            accessible[chrom][1],
            args.window_bp,
        )
        for chrom in contigs
    ]

    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(process_chromosome, task): task[1] for task in tasks}
        for future in as_completed(futures):
            chrom = futures[future]
            result = future.result()
            results.append(result)
            print(f"Completed chr{chrom}: {len(result)} windows", flush=True)

    order = {f"chr{chrom}": index for index, chrom in enumerate(contigs)}
    frame = pd.concat(results, ignore_index=True)
    frame["_order"] = frame["chrom"].map(order)
    frame = frame.sort_values(["_order", "start"]).drop(columns="_order")
    diagnostics = validate(frame, args.window_bp)

    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output, index=False, compression="gzip" if output.suffix == ".gz" else None)
    pd.DataFrame([diagnostics]).to_csv(args.validation_out, sep="\t", index=False)
    print(pd.Series(diagnostics).to_string(), flush=True)
    print(f"Wrote {output}", flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Generate all analysis figures for the IBD comparison project.

Produces:
  figures/seg_length_hist.png      - segment length distribution
  figures/pihat_vs_ibd.png         - PLINK PI_HAT vs IBD sharing
  figures/karyogram_<pair>.png     - per-pair chromosome karyograms
"""
import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# ── helpers ──────────────────────────────────────────────────────────────────

def norm_pair(a, b):
    return (a, b) if a <= b else (b, a)


def merge_intervals(ivs):
    if not ivs:
        return []
    ivs = sorted(ivs)
    m = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= m[-1][1]:
            m[-1][1] = max(m[-1][1], e)
        else:
            m.append([s, e])
    return [(a, b) for a, b in m]


def total_bp(ivs):
    return sum(e - s for s, e in merge_intervals(ivs))


def _iid(hap_id):
    if "_" in hap_id:
        parts = hap_id.rsplit("_", 1)
        if parts[1].isdigit():
            return parts[0]
    return hap_id


# ── parsers ───────────────────────────────────────────────────────────────────

def parse_germline_phased(path):
    """Return dict pair -> segments AND list of (pair, length_bp) for histogram."""
    pairs = defaultdict(list)
    seg_lengths = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            p1 = cols[0].split()
            p2 = cols[1].split()
            if len(p1) < 2 or len(p2) < 2:
                continue
            iid1 = _iid(p1[1])
            iid2 = _iid(p2[1])
            start, end = int(cols[3].split()[0]), int(cols[3].split()[1])
            pairs[norm_pair(iid1, iid2)].append((start, end))
            seg_lengths.append(end - start)
    return pairs, seg_lengths


def parse_beagle(path):
    pairs = defaultdict(list)
    seg_lengths = []
    with gzip.open(path, "rb") as f:
        for line in f:
            line = line.decode().strip()
            if not line:
                continue
            sid1, _h1, sid2, _h2, _chrom, start, end, _lod = line.split("\t")
            pairs[norm_pair(sid1, sid2)].append((int(start), int(end)))
            seg_lengths.append(int(end) - int(start))
    return pairs, seg_lengths


def parse_plink_genome(path):
    out = {}
    with open(path) as f:
        header = f.readline().strip().split()
        if not header:
            return out
        idx = {h: i for i, h in enumerate(header)}
        required = {"IID1", "IID2", "PI_HAT", "Z0", "Z1", "Z2"}
        if not required.issubset(idx):
            missing = sorted(required - set(idx))
            raise ValueError(f"PLINK genome file missing required columns: {missing}")
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            if len(parts) < len(header):
                continue
            key = norm_pair(parts[idx["IID1"]], parts[idx["IID2"]])
            out[key] = {
                "pi_hat": float(parts[idx["PI_HAT"]]),
                "z0": float(parts[idx["Z0"]]),
                "z1": float(parts[idx["Z1"]]),
                "z2": float(parts[idx["Z2"]]),
            }
    return out


# ── figure 1: segment length histogram ───────────────────────────────────────

def plot_seg_length_hist(germ_lengths, beagle_lengths, outpath):
    fig, ax = plt.subplots(figsize=(8, 4))
    bins = list(range(0, 12_000_000, 500_000))
    ax.hist(
        [l / 1e6 for l in germ_lengths],
        bins=[b / 1e6 for b in bins],
        alpha=0.65,
        color="#b08968",
        label=f"GERMLINE phased (n={len(germ_lengths)})",
        edgecolor="white",
        linewidth=0.4,
    )
    ax.hist(
        [l / 1e6 for l in beagle_lengths],
        bins=[b / 1e6 for b in bins],
        alpha=0.65,
        color="#0b7285",
        label=f"Beagle (n={len(beagle_lengths)})",
        edgecolor="white",
        linewidth=0.4,
    )
    ax.set_xlabel("IBD Segment Length (Mb)")
    ax.set_ylabel("Count")
    ax.set_title("Distribution of IBD Segment Lengths: GERMLINE vs Beagle")
    ax.legend()
    ax.grid(axis="y", alpha=0.25, linewidth=0.6)
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)
    print(f"Saved: {outpath}")


# ── figure 2: PLINK PI_HAT vs total IBD bp ───────────────────────────────────

def plot_pihat_vs_ibd(germ_pairs, beagle_pairs, plink, chr_len, outpath):
    all_pairs = sorted(set(plink.keys()) & (set(germ_pairs) | set(beagle_pairs)))

    germ_props = []
    beagle_props = []
    pi_hats = []

    for pair in all_pairs:
        if pair not in plink:
            continue
        pi_hats.append(plink[pair]["pi_hat"])
        germ_props.append(total_bp(germ_pairs.get(pair, [])) / chr_len)
        beagle_props.append(total_bp(beagle_pairs.get(pair, [])) / chr_len)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=False)

    for ax, props, label, color in [
        (axes[0], germ_props, "GERMLINE phased", "#b08968"),
        (axes[1], beagle_props, "Beagle", "#0b7285"),
    ]:
        ax.scatter(pi_hats, props, alpha=0.5, s=18, color=color, edgecolors="none")
        ax.set_xlabel("PLINK PI_HAT")
        ax.set_ylabel("Fraction of Chr22 in IBD")
        ax.set_title(f"{label}: PI_HAT vs Chr22 IBD sharing")
        ax.grid(alpha=0.2, linewidth=0.5)
        if not pi_hats:
            ax.text(0.5, 0.5, "No overlapping pairs", ha="center", va="center",
                    transform=ax.transAxes, fontsize=10)

    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)
    print(f"Saved: {outpath}")


# ── figure 3: karyogram for a pair ───────────────────────────────────────────

def plot_karyogram(germ_pairs, beagle_pairs, id1, id2, chr_len, outpath):
    pair = norm_pair(id1, id2)
    germ_segs = merge_intervals(germ_pairs.get(pair, []))
    beagle_segs = merge_intervals(beagle_pairs.get(pair, []))

    fig, ax = plt.subplots(figsize=(11, 3))

    # Chromosome backbone
    ax.barh(0.7, chr_len, left=0, height=0.08, color="#dddddd", zorder=0)
    ax.barh(0.3, chr_len, left=0, height=0.08, color="#dddddd", zorder=0)

    for s, e in germ_segs:
        ax.barh(0.7, e - s, left=s, height=0.16,
                color="#b08968", alpha=0.85, edgecolor="none", zorder=1)
    for s, e in beagle_segs:
        ax.barh(0.3, e - s, left=s, height=0.16,
                color="#0b7285", alpha=0.85, edgecolor="none", zorder=1)

    ax.set_xlim(0, chr_len)
    ax.set_ylim(0, 1)
    ax.set_yticks([0.3, 0.7])
    ax.set_yticklabels([
        f"Beagle (n={len(beagle_segs)}, {total_bp(beagle_segs)/1e6:.1f} Mb)",
        f"GERMLINE (n={len(germ_segs)}, {total_bp(germ_segs)/1e6:.1f} Mb)",
    ])
    ax.set_xlabel("Chr22 Position (bp)")
    ax.set_title(f"IBD Segments: {id1} vs {id2}")
    ax.grid(axis="x", alpha=0.2, linewidth=0.5)

    germ_patch = mpatches.Patch(color="#b08968", label="GERMLINE phased")
    beagle_patch = mpatches.Patch(color="#0b7285", label="Beagle")
    ax.legend(handles=[germ_patch, beagle_patch], loc="upper right", fontsize=8)

    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)
    print(f"Saved: {outpath}")


# ── figure 4: overlap visualisation ──────────────────────────────────────────

def plot_method_overlap_bar(overall_metrics_path, outpath):
    with open(overall_metrics_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        row = next(reader)

    b_bp = int(row["beagle_total_bp"])
    g_bp = int(row["germline_total_bp"])
    ov_bp = int(row["pairwise_overlap_bp"])
    b_only = b_bp - ov_bp
    g_only = g_bp - ov_bp

    categories = ["Beagle only", "Shared", "GERMLINE only"]
    values = [b_only / 1e6, ov_bp / 1e6, g_only / 1e6]
    colors = ["#0b7285", "#4c956c", "#b08968"]

    fig, ax = plt.subplots(figsize=(6, 4))
    bars = ax.bar(categories, values, color=colors, edgecolor="white", linewidth=0.5)
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.5,
                f"{val:.1f} Mb", ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("Total IBD (Mb)")
    ax.set_title("IBD Coverage by Method (Chr22, ASW Population)")
    ax.grid(axis="y", alpha=0.25, linewidth=0.6)
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)
    print(f"Saved: {outpath}")


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--beagle-ibd", required=True)
    p.add_argument("--germline-match", required=True)
    p.add_argument("--plink-genome", required=True)
    p.add_argument("--overall-metrics", required=True)
    p.add_argument("--chr-length", type=int, default=51_244_237)
    p.add_argument("--pairs", nargs="+", default=["NA20317:NA20318", "NA20359:NA20362"],
                   help="Sample pairs for karyograms (format: ID1:ID2)")
    p.add_argument("--outdir", default="results/figures")
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("Parsing GERMLINE...")
    germ_pairs, germ_lengths = parse_germline_phased(args.germline_match)
    print("Parsing Beagle...")
    beagle_pairs, beagle_lengths = parse_beagle(args.beagle_ibd)
    print("Parsing PLINK...")
    plink = parse_plink_genome(args.plink_genome)

    plot_seg_length_hist(germ_lengths, beagle_lengths, outdir / "seg_length_hist.png")
    plot_pihat_vs_ibd(germ_pairs, beagle_pairs, plink, args.chr_length, outdir / "pihat_vs_ibd.png")
    plot_method_overlap_bar(args.overall_metrics, outdir / "method_overlap.png")

    for pair_str in args.pairs:
        id1, id2 = pair_str.split(":")
        fname = outdir / f"karyogram_{id1}_{id2}.png"
        plot_karyogram(germ_pairs, beagle_pairs, id1, id2, args.chr_length, fname)


if __name__ == "__main__":
    main()

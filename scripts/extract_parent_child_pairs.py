#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from collections import Counter


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ped", required=True, help="Path to integrated_call_samples*.ped")
    parser.add_argument(
        "--out",
        required=True,
        help="Output TSV with columns: population,parent,child,relationship",
    )
    args = parser.parse_args()

    rows = []
    with open(args.ped, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            if r.get("phase 3 genotypes") == "1":
                rows.append(r)

    present = {r["Individual ID"] for r in rows}
    pairs: list[tuple[str, str, str, str]] = []
    for r in rows:
        pop = r["Population"]
        child = r["Individual ID"]
        dad = r["Paternal ID"]
        mom = r["Maternal ID"]
        if dad != "0" and dad in present:
            pairs.append((pop, dad, child, "father-child"))
        if mom != "0" and mom in present:
            pairs.append((pop, mom, child, "mother-child"))

    with open(args.out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["population", "parent", "child", "relationship"])
        w.writerows(pairs)

    counts = Counter(pop for pop, _, _, _ in pairs)
    print(f"pairs={len(pairs)} populations={len(counts)}")
    for pop, n in counts.most_common(10):
        print(f"{pop}\t{n}")


if __name__ == "__main__":
    main()

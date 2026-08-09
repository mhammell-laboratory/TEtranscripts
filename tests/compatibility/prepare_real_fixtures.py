#!/usr/bin/env python3
"""Reproduce the compact real RNA-seq count fixtures used in compatibility CI."""

import argparse
import csv
import gzip
import hashlib
import io
import math
import os
import urllib.request


FIXTURE_DIRECTORY = os.path.join(os.path.dirname(__file__), "fixtures")
FIXTURE_SIZE = 120
PER_STRATUM = 10

DATASETS = (
    {
        "name": "pasilla",
        "url": "https://raw.githubusercontent.com/thelovelab/DESeq2/devel/inst/extdata/pasilla_gene_counts.tsv.gz",
        "sha256": "4b5570008936786da365d110df1c83c8cffc646d6c667d6b15ecbfb3d9ee735d",
        "delimiter": "\t",
        "gzip": True,
        "treatment": ("treated1", "treated2", "treated3"),
        "control": ("untreated1", "untreated2", "untreated3"),
    },
    {
        "name": "airway",
        "url": "https://bdsr.stephenturner.us/data/airway_scaledcounts.csv",
        "sha256": "4d34df15affa22fa5d5539874e148378ea18db1c9d731d175ef1b9f73cbc2927",
        "delimiter": ",",
        "gzip": False,
        "treatment": ("SRR1039509", "SRR1039513", "SRR1039517"),
        "control": ("SRR1039508", "SRR1039512", "SRR1039516"),
    },
)


def _stable_key(dataset_name, feature_name):
    value = (dataset_name + "\0" + feature_name).encode("utf-8")
    return hashlib.sha256(value).hexdigest()


def _download(dataset):
    with urllib.request.urlopen(dataset["url"]) as response:
        content = response.read()
    digest = hashlib.sha256(content).hexdigest()
    if digest != dataset["sha256"]:
        raise ValueError(
            "%s source checksum changed: expected %s, got %s"
            % (dataset["name"], dataset["sha256"], digest)
        )
    return gzip.decompress(content) if dataset["gzip"] else content


def _read_rows(dataset, content):
    text = io.StringIO(content.decode("utf-8"))
    reader = csv.DictReader(text, delimiter=dataset["delimiter"])
    sample_names = dataset["treatment"] + dataset["control"]
    missing = set(sample_names) - set(reader.fieldnames or ())
    if missing:
        raise ValueError("%s source is missing samples: %s" % (dataset["name"], sorted(missing)))

    feature_column = reader.fieldnames[0]
    rows = []
    for source_row in reader:
        counts = []
        for sample_name in sample_names:
            value = float(source_row[sample_name])
            rounded = int(round(value))
            if value < 0 or abs(value - rounded) > 1e-6:
                raise ValueError("count sources must contain non-negative integers")
            counts.append(rounded)
        if max(counts) <= 1:
            continue
        treatment_mean = sum(counts[:3]) / 3.0
        control_mean = sum(counts[3:]) / 3.0
        log2_fold_change = math.log(
            (treatment_mean + 0.5) / (control_mean + 0.5), 2
        )
        rows.append(
            {
                "feature": source_row[feature_column],
                "counts": counts,
                "total": sum(counts),
                "effect": (
                    "down" if log2_fold_change < -1.5
                    else "up" if log2_fold_change > 1.5
                    else "stable"
                ),
            }
        )
    return rows


def _subsample(dataset, rows):
    by_depth = sorted(rows, key=lambda row: (row["total"], row["feature"]))
    for rank, row in enumerate(by_depth):
        row["depth"] = min(3, rank * 4 // len(by_depth))

    selected = []
    selected_names = set()
    for depth in range(4):
        for effect in ("down", "stable", "up"):
            candidates = [
                row for row in rows
                if row["depth"] == depth and row["effect"] == effect
            ]
            candidates.sort(key=lambda row: _stable_key(dataset["name"], row["feature"]))
            for row in candidates[:PER_STRATUM]:
                selected.append(row)
                selected_names.add(row["feature"])

    # Some high-depth directional strata contain fewer than ten features. Fill
    # those slots deterministically from all remaining eligible source rows.
    if len(selected) < FIXTURE_SIZE:
        remaining = [row for row in rows if row["feature"] not in selected_names]
        remaining.sort(key=lambda row: _stable_key(dataset["name"], row["feature"]))
        selected.extend(remaining[:FIXTURE_SIZE - len(selected)])
    if len(selected) != FIXTURE_SIZE:
        raise ValueError("%s source did not provide %d eligible rows" % (dataset["name"], FIXTURE_SIZE))
    return sorted(selected, key=lambda row: row["feature"])


def _render(dataset, rows):
    samples = [name + ".T" for name in dataset["treatment"]]
    samples.extend(name + ".C" for name in dataset["control"])
    lines = ["\t".join(["gene/TE"] + samples)]
    for row in rows:
        lines.append("\t".join([row["feature"]] + [str(value) for value in row["counts"]]))
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--check",
        action="store_true",
        help="verify that regenerated content matches the committed fixtures",
    )
    args = parser.parse_args()

    os.makedirs(FIXTURE_DIRECTORY, exist_ok=True)
    for dataset in DATASETS:
        rows = _read_rows(dataset, _download(dataset))
        rendered = _render(dataset, _subsample(dataset, rows))
        destination = os.path.join(
            FIXTURE_DIRECTORY, dataset["name"] + "_subsample.cntTable"
        )
        if args.check:
            with open(destination) as handle:
                committed = handle.read()
            if committed != rendered:
                raise SystemExit("fixture differs from regenerated source: " + destination)
        else:
            with open(destination, "w") as handle:
                handle.write(rendered)
        print("verified" if args.check else "wrote", destination)


if __name__ == "__main__":
    main()

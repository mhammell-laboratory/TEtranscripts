#!/usr/bin/env python3
"""Run the legacy DESeq2 path and native replacement on one count table."""

from __future__ import division

import argparse
import csv
import math
import os
import statistics
import subprocess

from TEToolkit.DifferentialAnalysis import run_differential_analysis


RESULT_COLUMNS = (
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
)


def _number(value):
    if value in (None, "", "NA", "NaN", "Inf", "-Inf"):
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def _read_results(filename):
    with open(filename, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        # R's historical write.table() output has an unnamed row-name field in
        # each record but no corresponding leading header cell.  The native
        # writer emits that empty header cell explicitly.  Normalize both forms
        # before mapping the numeric result columns.
        columns = tuple(header[1:] if header and not header[0] else header)
        rows = {}
        for fields in reader:
            if not fields:
                continue
            name = fields[0]
            values = fields[1:]
            row = dict(zip(columns, values))
            rows[name] = {column: _number(row.get(column)) for column in RESULT_COLUMNS}
    return columns, rows


def _ranks(values):
    ordered = sorted(enumerate(values), key=lambda item: item[1])
    result = [0.0] * len(values)
    index = 0
    while index < len(ordered):
        stop = index + 1
        while stop < len(ordered) and ordered[stop][1] == ordered[index][1]:
            stop += 1
        rank = (index + 1 + stop) / 2.0
        for original_index, _ in ordered[index:stop]:
            result[original_index] = rank
        index = stop
    return result


def _pearson(left, right):
    left_mean = statistics.mean(left)
    right_mean = statistics.mean(right)
    numerator = sum(
        (left_value - left_mean) * (right_value - right_mean)
        for left_value, right_value in zip(left, right)
    )
    left_sum = sum((value - left_mean) ** 2 for value in left)
    right_sum = sum((value - right_mean) ** 2 for value in right)
    denominator = math.sqrt(left_sum * right_sum)
    return numerator / denominator if denominator else 1.0


def _spearman(left, right):
    return _pearson(_ranks(left), _ranks(right))


def _fmt(value, digits=3):
    if value is None:
        return "NA"
    return ("%%.%df" % digits) % value


def _sign(value):
    return 1 if value > 0 else -1 if value < 0 else 0


def _significant(rows):
    return {
        name
        for name, row in rows.items()
        if row["padj"] is not None
        and row["padj"] < 0.05
        and row["log2FoldChange"] is not None
        and abs(row["log2FoldChange"]) > 1.0
    }


def _render_report(legacy_columns, legacy, native_columns, native, legacy_commit, r_version):
    names = sorted(set(legacy) & set(native))
    lfc_pairs = [
        (legacy[name]["log2FoldChange"], native[name]["log2FoldChange"])
        for name in names
        if legacy[name]["log2FoldChange"] is not None
        and native[name]["log2FoldChange"] is not None
    ]
    old_lfc = [pair[0] for pair in lfc_pairs]
    new_lfc = [pair[1] for pair in lfc_pairs]
    direction_pairs = [pair for pair in lfc_pairs if abs(pair[0]) >= 0.5]
    direction_matches = sum(_sign(old) == _sign(new) for old, new in direction_pairs)
    direction_rate = direction_matches / len(direction_pairs) if direction_pairs else 1.0

    padj_pairs = [
        (legacy[name]["padj"], native[name]["padj"])
        for name in names
        if legacy[name]["padj"] not in (None, 0.0)
        and native[name]["padj"] not in (None, 0.0)
    ]
    old_padj_rank = [-math.log10(pair[0]) for pair in padj_pairs]
    new_padj_rank = [-math.log10(pair[1]) for pair in padj_pairs]
    legacy_sig = _significant(legacy)
    native_sig = _significant(native)
    intersection = legacy_sig & native_sig
    union = legacy_sig | native_sig
    sig_jaccard = len(intersection) / len(union) if union else 1.0
    expected = {
        name: 1 if name.startswith("up_") else -1
        for name in names
        if name.startswith("up_") or name.startswith("down_")
    }
    expected_matches = sum(
        _sign(legacy[name]["log2FoldChange"]) == direction
        and _sign(native[name]["log2FoldChange"]) == direction
        for name, direction in expected.items()
    )

    contract_match = legacy_columns == native_columns == RESULT_COLUMNS
    feature_match = set(legacy) == set(native)
    lfc_spearman = _spearman(old_lfc, new_lfc)
    padj_spearman = _spearman(old_padj_rank, new_padj_rank) if padj_pairs else 1.0
    median_lfc_difference = statistics.median(
        abs(old - new) for old, new in lfc_pairs
    )

    lines = [
        "# Legacy DESeq2 versus Python-native comparison",
        "",
        "The legacy result was generated by the exact `write_R_code()` function ",
        "from commit `%s`; its emitted script was run with %s. Both paths used " % (legacy_commit, r_version),
        "`small_counts.cntTable` with three treatment and three control replicates.",
        "",
        "This is a compatibility test, not a claim of numerical identity with the full DESeq2 package.",
        "",
        "| Measure | Result |",
        "| --- | ---: |",
        "| Result-table contract | %s |" % ("match" if contract_match else "DIFFERENT"),
        "| Feature rows | %d legacy / %d native |" % (len(legacy), len(native)),
        "| Feature-set match | %s |" % ("yes" if feature_match else "no"),
        "| Log2 fold-change Spearman correlation | %s |" % _fmt(lfc_spearman),
        "| Median absolute log2 fold-change difference | %s |" % _fmt(median_lfc_difference),
        "| Direction agreement for DESeq2 abs(LFC) >= 0.5 | %d/%d (%s%%) |"
        % (direction_matches, len(direction_pairs), _fmt(direction_rate * 100, 1)),
        "| Engineered up/down direction recovered by both | %d/%d |"
        % (expected_matches, len(expected)),
        "| Adjusted-p-value rank correlation | %s |" % _fmt(padj_spearman),
        "| Significant calls (padj < 0.05 and abs(LFC) > 1) | %d legacy / %d native |"
        % (len(legacy_sig), len(native_sig)),
        "| Significant-call overlap | %d; Jaccard %s |"
        % (len(intersection), _fmt(sig_jaccard)),
        "",
        "## Per-feature results",
        "",
        "| Feature | DESeq2 LFC | Native LFC | DESeq2 padj | Native padj |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for name in names:
        lines.append(
            "| %s | %s | %s | %s | %s |"
            % (
                name,
                _fmt(legacy[name]["log2FoldChange"]),
                _fmt(native[name]["log2FoldChange"]),
                _fmt(legacy[name]["padj"]),
                _fmt(native[name]["padj"]),
            )
        )

    failures = []
    if not contract_match:
        failures.append("result-table columns differ")
    if not feature_match:
        failures.append("feature sets differ")
    if lfc_spearman < 0.90:
        failures.append("fold-change rank correlation is below 0.90")
    if direction_rate < 0.90:
        failures.append("fold-change direction agreement is below 90%")
    if expected_matches != len(expected):
        failures.append("an engineered up/down effect has the wrong direction")
    lines.extend(("", "## Acceptance checks", ""))
    if failures:
        lines.extend("- FAIL: " + failure for failure in failures)
    else:
        lines.extend(
            (
                "- PASS: identical result-table and feature contracts",
                "- PASS: fold-change rank correlation >= 0.90",
                "- PASS: fold-change direction agreement >= 90%",
                "- PASS: both implementations recover every engineered up/down direction",
            )
        )
    return "\n".join(lines) + "\n", failures


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--legacy-revision", required=True)
    parser.add_argument("--counts", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    counts = os.path.abspath(args.counts)
    legacy_project = os.path.join(os.path.abspath(args.output_dir), "legacy")
    native_project = os.path.join(os.path.abspath(args.output_dir), "native")
    legacy_path = args.legacy_revision + ":bin/TEtranscripts"
    legacy_source = subprocess.check_output(
        ["git", "show", legacy_path], text=True
    )
    legacy_namespace = {
        "__file__": legacy_path,
        "__name__": "legacy_TEtranscripts",
    }
    exec(compile(legacy_source, legacy_path, "exec"), legacy_namespace)
    r_code = legacy_namespace["write_R_code"](
        counts,
        ["t1", "t2", "t3"],
        ["c1", "c2", "c3"],
        legacy_project,
        False,
        "DESeq_default",
        0.05,
        2.0,
        [],
        1,
    )
    legacy_r = os.path.join(args.output_dir, "legacy_DESeq2.R")
    with open(legacy_r, "w") as handle:
        handle.write(r_code)
    subprocess.check_call(["Rscript", legacy_r])

    run_differential_analysis(
        counts,
        3,
        3,
        native_project,
        padj_cutoff=0.05,
        fold_change_cutoff=2.0,
        min_read=1,
    )
    legacy_columns, legacy = _read_results(legacy_project + "_gene_TE_analysis.txt")
    native_columns, native = _read_results(native_project + "_gene_TE_analysis.txt")
    r_version = subprocess.check_output(
        ["Rscript", "--version"], stderr=subprocess.STDOUT, text=True
    ).strip()
    report, failures = _render_report(
        legacy_columns,
        legacy,
        native_columns,
        native,
        args.legacy_revision,
        r_version,
    )
    with open(args.report, "w") as handle:
        handle.write(report)
    print(report)
    if failures:
        raise SystemExit("compatibility acceptance checks failed")


if __name__ == "__main__":
    main()

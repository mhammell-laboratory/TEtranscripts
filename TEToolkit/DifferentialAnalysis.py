"""Python-native differential expression helpers for TEtranscripts.

The public functions in this module mirror the small subset of the DESeq and
DESeq2 APIs that TEtranscripts historically used.  They intentionally avoid R
and Bioconductor while preserving TEtranscripts' command-line inputs and output
file contracts.

This is a compact, dependency-free implementation of the DESeq median-ratio
normalization, mean-dispersion fitting, negative-binomial-inspired exact tests,
Wald tests, Benjamini-Hochberg adjustment, and result-table generation.  It is
not intended to be a general replacement for the full Bioconductor packages.
"""

from __future__ import division

import csv
import math


_EPSILON = 1e-12
_LOG_2 = math.log(2.0)

try:
    _STRING_TYPES = (basestring,)
except NameError:
    _STRING_TYPES = (str,)


def _isfinite(value):
    return not math.isnan(value) and not math.isinf(value)


class CountDataSet(object):
    """Minimal Python equivalent of the DESeq ``CountDataSet`` container."""

    def __init__(self, count_data, conditions, row_names=None, column_names=None):
        self.counts = [list(map(float, row)) for row in count_data]
        if not self.counts:
            raise ValueError("count_data must contain at least one feature")

        width = len(self.counts[0])
        if width == 0 or any(len(row) != width for row in self.counts):
            raise ValueError("count_data must be a non-empty rectangular matrix")
        if len(conditions) != width:
            raise ValueError("one condition is required for each sample")
        if any((not _isfinite(value)) or value < 0 for row in self.counts for value in row):
            raise ValueError("counts must be finite and non-negative")

        self.conditions = list(conditions)
        self.row_names = list(row_names or ("row%d" % (i + 1) for i in range(len(self.counts))))
        self.column_names = list(column_names or ("sample%d" % (i + 1) for i in range(width)))
        if len(self.row_names) != len(self.counts):
            raise ValueError("one row name is required for each feature")
        if len(self.column_names) != width:
            raise ValueError("one column name is required for each sample")

        self.size_factors = [1.0] * width
        self.normalized_counts = [row[:] for row in self.counts]
        self.base_means = [_mean(row) for row in self.counts]
        self.raw_dispersions = [0.0] * len(self.counts)
        self.fitted_dispersions = [0.0] * len(self.counts)
        self.dispersions = [0.0] * len(self.counts)
        self.dispersion_method = None


def newCountDataSet(count_data, conditions, row_names=None, column_names=None):
    """Create a :class:`CountDataSet`, matching the DESeq constructor used here."""

    return CountDataSet(count_data, conditions, row_names, column_names)


def _mean(values):
    return sum(values) / len(values) if values else 0.0


def _median(values):
    ordered = sorted(values)
    size = len(ordered)
    if not size:
        raise ValueError("median requires at least one value")
    middle = size // 2
    if size % 2:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / 2.0


def _geometric_mean(values):
    if not values or any(value <= 0 for value in values):
        return 0.0
    return math.exp(sum(math.log(value) for value in values) / len(values))


def estimateSizeFactors(cds, method="ratio"):
    """Estimate DESeq median-ratio size factors in place and return ``cds``.

    ``method='ratio'`` implements the DESeq default.  If every feature contains
    a zero, the positive-count geometric mean is used as a sparse-count fallback
    so TE-only matrices remain analyzable without changing the command line.
    ``method='poscounts'`` requests that behavior explicitly.
    """

    if method not in ("ratio", "poscounts"):
        raise ValueError("unsupported size-factor method: %s" % method)

    use_positive = method == "poscounts"
    geometric_means = []
    for row in cds.counts:
        values = [value for value in row if value > 0] if use_positive else row
        geometric_means.append(_geometric_mean(values))

    if not any(value > 0 for value in geometric_means):
        use_positive = True
        geometric_means = [_geometric_mean([value for value in row if value > 0]) for row in cds.counts]

    factors = []
    for sample_index in range(len(cds.column_names)):
        ratios = []
        for row, geo_mean in zip(cds.counts, geometric_means):
            count = row[sample_index]
            if geo_mean > 0 and (count > 0 or not use_positive):
                ratios.append(count / geo_mean)
        if not ratios:
            if all(row[sample_index] == 0 for row in cds.counts):
                raise ValueError(
                    "cannot estimate a size factor for sample %s; "
                    "the sample library contains only zero counts"
                    % cds.column_names[sample_index]
                )
            raise ValueError("cannot estimate a size factor for sample %s" % cds.column_names[sample_index])
        factor = _median(ratios)
        if factor <= 0 or not _isfinite(factor):
            raise ValueError("estimated size factors must be positive and finite")
        factors.append(factor)

    # Size factors are identifiable only up to a common multiplier.  Centering
    # them on a geometric mean of one makes results stable across input scales.
    center = _geometric_mean(factors)
    cds.size_factors = [factor / center for factor in factors]
    _update_normalized_counts(cds)
    return cds


def _set_size_factors(cds, factors):
    if len(factors) != len(cds.column_names):
        raise ValueError("one size factor is required for each sample")
    factors = [float(value) for value in factors]
    if any(value <= 0 or not _isfinite(value) for value in factors):
        raise ValueError("size factors must be positive and finite")
    cds.size_factors = factors
    _update_normalized_counts(cds)
    return cds


def _update_normalized_counts(cds):
    cds.normalized_counts = [
        [value / factor for value, factor in zip(row, cds.size_factors)]
        for row in cds.counts
    ]
    cds.base_means = [_mean(row) for row in cds.normalized_counts]


def _condition_indices(cds):
    result = {}
    for index, condition in enumerate(cds.conditions):
        result.setdefault(condition, []).append(index)
    return result


def _raw_dispersion(row, factors, indices, residual_df=None, group_indices=None):
    if residual_df is None:
        residual_df = len(indices) - 1
    if residual_df <= 0:
        return None

    normalized = [row[index] / factors[index] for index in indices]
    mean_value = _mean(normalized)
    if mean_value <= 0:
        return 0.0

    if group_indices:
        residual_sum = 0.0
        for group in group_indices:
            group_values = [row[index] / factors[index] for index in group]
            group_mean = _mean(group_values)
            residual_sum += sum((value - group_mean) ** 2 for value in group_values)
    else:
        residual_sum = sum((value - mean_value) ** 2 for value in normalized)

    variance = residual_sum / residual_df
    poisson_variance = mean_value * _mean([1.0 / factors[index] for index in indices])
    return max((variance - poisson_variance) / max(mean_value ** 2, _EPSILON), 0.0)


def _local_dispersion_fit(means, raw_dispersions):
    valid = [
        (math.log(mean_value), math.log(dispersion))
        for mean_value, dispersion in zip(means, raw_dispersions)
        if mean_value > 0 and dispersion is not None and dispersion > 0 and _isfinite(dispersion)
    ]
    if not valid:
        fallback = 0.1
        return [fallback if mean_value > 0 else 0.0 for mean_value in means]

    valid.sort()
    window = min(len(valid), max(5, int(round(math.sqrt(len(valid))))))
    global_fit = math.exp(_median([item[1] for item in valid]))
    fitted = []
    for mean_value in means:
        if mean_value <= 0:
            fitted.append(0.0)
            continue
        target = math.log(mean_value)
        nearest = sorted(valid, key=lambda item: abs(item[0] - target))[:window]
        estimate = math.exp(_median([item[1] for item in nearest])) if nearest else global_fit
        fitted.append(max(estimate, 1e-8))
    return fitted


def _parametric_dispersion_fit(means, raw_dispersions):
    """Fit the DESeq-style ``a0 + a1 / mean`` dispersion relationship."""

    valid = [
        (1.0 / mean_value, dispersion)
        for mean_value, dispersion in zip(means, raw_dispersions)
        if mean_value > 0 and dispersion is not None and dispersion > 0 and _isfinite(dispersion)
    ]
    if len(valid) < 2:
        return _local_dispersion_fit(means, raw_dispersions)

    x_mean = _mean([item[0] for item in valid])
    y_mean = _mean([item[1] for item in valid])
    denominator = sum((item[0] - x_mean) ** 2 for item in valid)
    if denominator <= 0:
        return _local_dispersion_fit(means, raw_dispersions)
    a1 = sum((x_value - x_mean) * (y_value - y_mean) for x_value, y_value in valid) / denominator
    a0 = y_mean - a1 * x_mean
    a0 = max(a0, 1e-8)
    a1 = max(a1, 0.0)
    return [
        max(a0 + a1 / mean_value, 1e-8) if mean_value > 0 else 0.0
        for mean_value in means
    ]


def estimateDispersions(cds, method="pooled", sharingMode="maximum", fitType="parametric"):
    """Estimate feature dispersions using the modes TEtranscripts relied on.

    Supported methods are ``blind``, ``pooled``, and ``per-condition``.  The
    historical ``fit-only`` and ``maximum`` sharing modes are preserved.
    ``shrinkage`` supplies a lightweight log-scale compromise between gene-wise
    and fitted estimates for the modern DESeq2-compatible path. Both the
    parametric mean relationship and robust local smoothing are available.
    """

    if method not in ("blind", "pooled", "per-condition"):
        raise ValueError("unsupported dispersion method: %s" % method)
    if sharingMode not in ("maximum", "fit-only", "gene-est-only", "shrinkage"):
        raise ValueError("unsupported sharing mode: %s" % sharingMode)
    if fitType not in ("local", "parametric", "mean"):
        raise ValueError("unsupported dispersion fit type: %s" % fitType)

    groups = _condition_indices(cds)
    all_indices = list(range(len(cds.column_names)))
    raw = []
    for row in cds.counts:
        if method == "blind":
            estimate = _raw_dispersion(row, cds.size_factors, all_indices)
        elif method == "pooled":
            residual_df = len(all_indices) - len(groups)
            estimate = _raw_dispersion(
                row,
                cds.size_factors,
                all_indices,
                residual_df=residual_df,
                group_indices=list(groups.values()),
            )
        else:
            estimates = [
                _raw_dispersion(row, cds.size_factors, indices)
                for indices in groups.values()
                if len(indices) > 1
            ]
            estimates = [value for value in estimates if value is not None]
            estimate = max(estimates) if estimates else None
        raw.append(estimate)

    fit_input = [0.0 if value is None else value for value in raw]
    if fitType == "mean":
        positive = [value for value in fit_input if value > 0]
        common = _median(positive) if positive else 0.1
        fitted = [common if mean_value > 0 else 0.0 for mean_value in cds.base_means]
    elif fitType == "parametric":
        fitted = _parametric_dispersion_fit(cds.base_means, fit_input)
    else:
        fitted = _local_dispersion_fit(cds.base_means, fit_input)

    if sharingMode == "fit-only":
        final = fitted[:]
    elif sharingMode == "gene-est-only":
        final = [fit if value is None else max(value, 1e-8) for value, fit in zip(raw, fitted)]
    elif sharingMode == "shrinkage":
        final = [
            fit if value is None or value <= 0 else math.sqrt(value * fit)
            for value, fit in zip(raw, fitted)
        ]
    else:
        final = [fit if value is None else max(value, fit) for value, fit in zip(raw, fitted)]

    cds.raw_dispersions = raw
    cds.fitted_dispersions = fitted
    cds.dispersions = final
    cds.dispersion_method = method
    return cds


def _benjamini_hochberg(pvalues):
    adjusted = [None] * len(pvalues)
    valid = [(value, index) for index, value in enumerate(pvalues) if value is not None and _isfinite(value)]
    valid.sort(reverse=True)
    total = len(valid)
    running = 1.0
    for reverse_rank, (value, index) in enumerate(valid):
        rank = total - reverse_rank
        running = min(running, value * total / rank)
        adjusted[index] = min(max(running, 0.0), 1.0)
    return adjusted


def _log_choose(n, k):
    return math.lgamma(n + 1.0) - math.lgamma(k + 1.0) - math.lgamma(n - k + 1.0)


def _log_beta_binomial_pmf(k, n, alpha, beta):
    return (
        _log_choose(n, k)
        + math.lgamma(k + alpha)
        + math.lgamma(n - k + beta)
        - math.lgamma(n + alpha + beta)
        + math.lgamma(alpha + beta)
        - math.lgamma(alpha)
        - math.lgamma(beta)
    )


def _log_binomial_pmf(k, n, probability):
    if probability <= 0:
        return 0.0 if k == 0 else float("-inf")
    if probability >= 1:
        return 0.0 if k == n else float("-inf")
    return _log_choose(n, k) + k * math.log(probability) + (n - k) * math.log1p(-probability)


def _two_sided_count_test(observed, total, probability, dispersion):
    if total <= 0:
        return 1.0
    observed = int(round(observed))
    total = int(round(total))
    probability = min(max(probability, _EPSILON), 1.0 - _EPSILON)

    if dispersion > 1e-10:
        concentration = max(1.0 / dispersion, 1e-6)
        alpha = probability * concentration
        beta = (1.0 - probability) * concentration
        log_pmf = lambda value: _log_beta_binomial_pmf(value, total, alpha, beta)
        variance = total * alpha * beta * (alpha + beta + total) / (
            (alpha + beta) ** 2 * (alpha + beta + 1.0)
        )
    else:
        log_pmf = lambda value: _log_binomial_pmf(value, total, probability)
        variance = total * probability * (1.0 - probability)

    # Exact two-sided enumeration is inexpensive for the low/moderate feature
    # counts where discreteness matters.  The beta-binomial can remain strongly
    # skewed or heavy-tailed at large totals, so a normal approximation is not
    # reliable there; retain exact recurrence enumeration for dispersed data.
    # The continuity-corrected normal tail is limited to the ordinary binomial
    # case above the historical cutoff.
    if total <= 10000 or dispersion > 1e-10:
        observed_log_probability = log_pmf(observed)
        probability_sum = 0.0
        if dispersion > 1e-10:
            current_log_probability = _log_beta_binomial_pmf(
                0, total, alpha, beta
            )
        else:
            current_log_probability = total * math.log1p(-probability)

        # Walk the PMF using P(k + 1) / P(k).  This preserves exact
        # enumeration while avoiding thousands of repeated lgamma calls near
        # the historical 10,000-count cutoff.
        for value in range(total + 1):
            if current_log_probability <= observed_log_probability + 1e-12:
                probability_sum += math.exp(current_log_probability)
            if value == total:
                break
            if dispersion > 1e-10:
                current_log_probability += (
                    math.log(total - value)
                    - math.log(value + 1)
                    + math.log(value + alpha)
                    - math.log(total - value - 1 + beta)
                )
            else:
                current_log_probability += (
                    math.log(total - value)
                    - math.log(value + 1)
                    + math.log(probability)
                    - math.log1p(-probability)
                )
        return min(max(probability_sum, 0.0), 1.0)

    if variance <= 0:
        return 1.0
    distance = max(abs(observed - total * probability) - 0.5, 0.0)
    return min(math.erfc(distance / math.sqrt(2.0 * variance)), 1.0)


def _safe_fold_change(numerator, denominator):
    if denominator == 0:
        if numerator == 0:
            return float("nan"), float("nan")
        return float("inf"), float("inf")
    fold_change = numerator / denominator
    if fold_change == 0:
        return 0.0, float("-inf")
    return fold_change, math.log(fold_change, 2)


def nbinomTest(cds, condition_a, condition_b):
    """Test two conditions and return DESeq-compatible result dictionaries."""

    groups = _condition_indices(cds)
    if condition_a not in groups or condition_b not in groups:
        raise ValueError("both requested conditions must be present")
    if not any(cds.dispersions):
        raise ValueError("estimateDispersions must be called before nbinomTest")

    a_indices = groups[condition_a]
    b_indices = groups[condition_b]
    exposure_a = sum(cds.size_factors[index] for index in a_indices)
    exposure_b = sum(cds.size_factors[index] for index in b_indices)
    null_probability_b = exposure_b / (exposure_a + exposure_b)

    rows = []
    pvalues = []
    for name, counts, normalized, dispersion in zip(
        cds.row_names, cds.counts, cds.normalized_counts, cds.dispersions
    ):
        mean_a = _mean([normalized[index] for index in a_indices])
        mean_b = _mean([normalized[index] for index in b_indices])
        base_mean = _mean(normalized)
        fold_change, log2_fold_change = _safe_fold_change(mean_b, mean_a)
        count_a = sum(counts[index] for index in a_indices)
        count_b = sum(counts[index] for index in b_indices)
        pvalue = _two_sided_count_test(count_b, count_a + count_b, null_probability_b, dispersion)
        pvalues.append(pvalue)
        rows.append({
            "id": name,
            "baseMean": base_mean,
            "baseMeanA": mean_a,
            "baseMeanB": mean_b,
            "foldChange": fold_change,
            "log2FoldChange": log2_fold_change,
            "pval": pvalue,
            "padj": None,
        })

    for row, adjusted in zip(rows, _benjamini_hochberg(pvalues)):
        row["padj"] = adjusted
    return rows


def _wald_results(cds, control_condition, treatment_condition):
    groups = _condition_indices(cds)
    control_indices = groups[control_condition]
    treatment_indices = groups[treatment_condition]
    rows = []
    pvalues = []

    for name, normalized, dispersion in zip(cds.row_names, cds.normalized_counts, cds.dispersions):
        control_mean = _mean([normalized[index] for index in control_indices])
        treatment_mean = _mean([normalized[index] for index in treatment_indices])
        base_mean = _mean(normalized)

        # Half a normalized count provides a finite boundary estimate when one
        # group is all zero, analogous to a weak count-scale prior.
        pseudo = 0.5 / max(_mean(cds.size_factors), _EPSILON)
        log2_fold_change = math.log((treatment_mean + pseudo) / (control_mean + pseudo), 2)
        control_variance = (
            _mean([1.0 / cds.size_factors[index] for index in control_indices])
            / (max(control_mean, pseudo) * len(control_indices))
            + dispersion / len(control_indices)
        )
        treatment_variance = (
            _mean([1.0 / cds.size_factors[index] for index in treatment_indices])
            / (max(treatment_mean, pseudo) * len(treatment_indices))
            + dispersion / len(treatment_indices)
        )
        lfc_se = math.sqrt(max(control_variance + treatment_variance, _EPSILON)) / _LOG_2
        statistic = log2_fold_change / lfc_se
        pvalue = min(math.erfc(abs(statistic) / math.sqrt(2.0)), 1.0)
        pvalues.append(pvalue)
        rows.append({
            "id": name,
            "baseMean": base_mean,
            "log2FoldChange": log2_fold_change,
            "lfcSE": lfc_se,
            "stat": statistic,
            "pvalue": pvalue,
            "padj": None,
        })

    for row, adjusted in zip(rows, _benjamini_hochberg(pvalues)):
        row["padj"] = adjusted
    return rows


def _quantile_normalize(counts):
    row_count = len(counts)
    column_count = len(counts[0])
    sorted_columns = []
    for column in range(column_count):
        sorted_columns.append(sorted((counts[row][column], row) for row in range(row_count)))
    rank_means = [
        _mean([sorted_columns[column][rank][0] for column in range(column_count)])
        for rank in range(row_count)
    ]
    result = [[0.0] * column_count for _ in range(row_count)]
    for column in range(column_count):
        rank = 0
        while rank < row_count:
            stop = rank + 1
            value = sorted_columns[column][rank][0]
            while stop < row_count and sorted_columns[column][stop][0] == value:
                stop += 1
            tied_mean = _mean(rank_means[rank:stop])
            for _, original_row in sorted_columns[column][rank:stop]:
                result[original_row][column] = tied_mean
            rank = stop
    return result


def _read_count_table(filename, treatment_count, control_count, min_read):
    names = []
    counts = []
    with open(filename, "r") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError("count table is empty")
        expected = treatment_count + control_count
        if len(header) != expected + 1:
            raise ValueError("count table has %d samples; expected %d" % (len(header) - 1, expected))
        for line_number, fields in enumerate(reader, 2):
            if not fields or all(not value for value in fields):
                continue
            if len(fields) != expected + 1:
                raise ValueError("count table line %d has the wrong number of columns" % line_number)
            try:
                row = [float(value) for value in fields[1:]]
            except ValueError:
                raise ValueError("count table line %d contains a non-numeric count" % line_number)
            if max(row) > min_read:
                names.append(fields[0])
                counts.append(row)
    if not counts:
        raise ValueError("no features remain after applying --minread %s" % min_read)
    return header[1:], names, counts


def _format_value(value):
    if isinstance(value, _STRING_TYPES):
        return value
    if value is None:
        return "NA"
    if math.isnan(value):
        return "NaN"
    if value == float("inf"):
        return "Inf"
    if value == float("-inf"):
        return "-Inf"
    return "%.12g" % value


def _write_results(filename, rows, columns, row_names=True):
    with open(filename, "w") as handle:
        handle.write("\t".join(([""] if row_names else []) + columns) + "\n")
        for row in rows:
            values = [_format_value(row[column]) for column in columns]
            handle.write("\t".join(([row["id"]] if row_names else []) + values) + "\n")


def run_differential_analysis(
    count_table,
    treatment_count,
    control_count,
    project_name,
    legacy_deseq=False,
    normalization="DESeq_default",
    padj_cutoff=0.05,
    fold_change_cutoff=1.0,
    library_sizes=None,
    min_read=1,
):
    """Run the Python-native analysis and write the historical result files."""

    if treatment_count <= 0:
        raise ValueError("treatment_count must be positive")
    if control_count <= 0:
        raise ValueError("control_count must be positive")

    column_names, row_names, counts = _read_count_table(
        count_table, treatment_count, control_count, min_read
    )
    conditions = ["TGroup"] * treatment_count + ["CGroup"] * control_count
    cds = newCountDataSet(counts, conditions, row_names=row_names, column_names=column_names)

    if legacy_deseq and normalization == "TC":
        if not library_sizes or len(library_sizes) != len(column_names):
            raise ValueError("TC normalization requires one library size per sample")
        minimum = min(library_sizes)
        if minimum <= 0:
            raise ValueError("library sizes must be positive for TC normalization")
        _set_size_factors(cds, [value / minimum for value in library_sizes])
    else:
        estimateSizeFactors(cds)

    if treatment_count == 1 and control_count == 1:
        estimateDispersions(cds, method="blind", sharingMode="fit-only", fitType="local")
    elif legacy_deseq and treatment_count > 1 and control_count > 1:
        estimateDispersions(cds, method="per-condition")
    else:
        estimateDispersions(cds, method="pooled", sharingMode="shrinkage")

    analysis_filename = project_name + "_gene_TE_analysis.txt"
    significant_filename = project_name + "_sigdiff_gene_TE.txt"
    threshold = math.log(fold_change_cutoff, 2)

    if not legacy_deseq:
        rows = _wald_results(cds, "CGroup", "TGroup")
        columns = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
        significant = [
            row for row in rows
            if row["padj"] is not None
            and row["padj"] < padj_cutoff
            and abs(row["log2FoldChange"]) > threshold
        ]
        _write_results(analysis_filename, rows, columns, row_names=True)
        _write_results(significant_filename, significant, columns, row_names=True)
    else:
        rows = nbinomTest(cds, "CGroup", "TGroup")
        if normalization == "quant":
            quantile_counts = _quantile_normalize(cds.counts)
            treatment_indices = list(range(treatment_count))
            control_indices = list(range(treatment_count, treatment_count + control_count))
            for row, normalized in zip(rows, quantile_counts):
                sample1_mean = _mean([normalized[index] for index in treatment_indices])
                sample2_mean = _mean([normalized[index] for index in control_indices])
                fold_change, log2_fold_change = _safe_fold_change(sample2_mean, sample1_mean)
                row.update({
                    "sample1Mean": sample1_mean,
                    "sample2Mean": sample2_mean,
                    "FoldChange": fold_change,
                    "log2FoldChange": log2_fold_change,
                })
            columns = ["id", "sample1Mean", "sample2Mean", "FoldChange", "log2FoldChange", "pval", "padj"]
            significant = [
                row for row in rows
                if row["padj"] is not None
                and row["padj"] < padj_cutoff
                and abs(row["log2FoldChange"]) > threshold
            ]
            _write_results(analysis_filename, rows, columns, row_names=False)
            _write_results(significant_filename, significant, columns, row_names=False)
        else:
            columns = ["id", "baseMean", "baseMeanA", "baseMeanB", "foldChange", "log2FoldChange", "pval", "padj"]
            significant = [
                row for row in rows
                if row["padj"] is not None
                and row["padj"] < padj_cutoff
                and abs(row["log2FoldChange"]) > threshold
            ]
            _write_results(analysis_filename, rows, columns, row_names=False)
            _write_results(significant_filename, significant, columns, row_names=False)

    return analysis_filename, significant_filename

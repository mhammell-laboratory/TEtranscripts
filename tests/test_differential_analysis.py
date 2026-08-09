from __future__ import division

import math
import os
import tempfile
import unittest

from TEToolkit.DifferentialAnalysis import (
    _benjamini_hochberg,
    _log_beta_binomial_pmf,
    _log_binomial_pmf,
    _quantile_normalize,
    _set_size_factors,
    _two_sided_count_test,
    estimateDispersions,
    estimateSizeFactors,
    nbinomTest,
    newCountDataSet,
    run_differential_analysis,
)


class DifferentialAnalysisTests(unittest.TestCase):
    def test_count_dataset_rejects_invalid_inputs(self):
        invalid = (
            ([], [], "at least one feature"),
            ([[1, 2], [3]], ["A", "B"], "rectangular"),
            ([[1, -1]], ["A", "B"], "non-negative"),
            ([[1, float("inf")]], ["A", "B"], "finite"),
            ([[1, 2]], ["A"], "condition"),
        )
        for counts, conditions, message in invalid:
            with self.subTest(message=message):
                with self.assertRaisesRegex(ValueError, message):
                    newCountDataSet(counts, conditions)

        with self.assertRaisesRegex(ValueError, "row name"):
            newCountDataSet([[1, 2]], ["A", "B"], row_names=["one", "two"])
        with self.assertRaisesRegex(ValueError, "column name"):
            newCountDataSet([[1, 2]], ["A", "B"], column_names=["only-one"])

    def test_median_ratio_size_factors(self):
        cds = newCountDataSet(
            [[10, 20], [30, 60], [5, 10]],
            ["A", "B"],
            row_names=["one", "two", "three"],
        )

        estimateSizeFactors(cds)

        self.assertAlmostEqual(cds.size_factors[0], 2 ** -0.5)
        self.assertAlmostEqual(cds.size_factors[1], 2 ** 0.5)
        for row in cds.normalized_counts:
            self.assertAlmostEqual(row[0], row[1])

    def test_size_factors_support_sparse_te_matrices(self):
        cds = newCountDataSet([[0, 10], [20, 0]], ["A", "B"])
        estimateSizeFactors(cds)
        self.assertEqual(len(cds.size_factors), 2)
        self.assertTrue(all(value > 0 for value in cds.size_factors))

    def test_size_factors_reject_an_all_zero_matrix(self):
        cds = newCountDataSet(
            [[0, 0], [0, 0]],
            ["T", "C"],
            column_names=["zero_treatment", "zero_control"],
        )
        with self.assertRaisesRegex(ValueError, "only zero counts"):
            estimateSizeFactors(cds)

    def test_size_factors_reject_an_entirely_zero_sample_library(self):
        cds = newCountDataSet(
            [[0, 10, 12], [0, 25, 20], [0, 5, 8]],
            ["T", "T", "C"],
            column_names=["empty_sample", "treatment", "control"],
        )
        with self.assertRaisesRegex(ValueError, "empty_sample.*only zero counts"):
            estimateSizeFactors(cds)

    def test_size_factor_validation(self):
        cds = newCountDataSet([[1, 2], [3, 4]], ["A", "B"])
        with self.assertRaisesRegex(ValueError, "unsupported"):
            estimateSizeFactors(cds, method="unknown")
        with self.assertRaisesRegex(ValueError, "one size factor"):
            _set_size_factors(cds, [1])
        for factors in ([0, 1], [-1, 1], [float("nan"), 1]):
            with self.subTest(factors=factors):
                with self.assertRaisesRegex(ValueError, "positive and finite"):
                    _set_size_factors(cds, factors)

    def test_benjamini_hochberg_preserves_order_and_missing_values(self):
        self.assertEqual(
            _benjamini_hochberg([0.01, 0.02, 0.04, None]),
            [0.03, 0.03, 0.04, None],
        )

    def test_quantile_normalization_equalizes_column_distributions(self):
        normalized = _quantile_normalize([[5, 2], [1, 4], [3, 6]])
        self.assertEqual(normalized, [[5.5, 1.5], [1.5, 3.5], [3.5, 5.5]])
        self.assertEqual(
            sorted(row[0] for row in normalized),
            sorted(row[1] for row in normalized),
        )

    def test_quantile_normalization_averages_tied_ranks(self):
        counts = [[1, 1], [1, 2], [3, 3]]
        self.assertEqual(
            _quantile_normalize(counts),
            [[1.25, 1.0], [1.25, 1.5], [3.0, 3.0]],
        )
        permuted = [counts[index] for index in (2, 0, 1)]
        normalized = _quantile_normalize(permuted)
        self.assertEqual([normalized[index] for index in (1, 2, 0)],
                         [[1.25, 1.0], [1.25, 1.5], [3.0, 3.0]])

    def test_all_historical_dispersion_modes_are_supported(self):
        counts = [
            [10, 12, 9, 11],
            [100, 130, 20, 25],
            [5, 8, 40, 55],
            [50, 44, 48, 51],
        ]
        for method in ("blind", "pooled", "per-condition"):
            cds = newCountDataSet(counts, ["T", "T", "C", "C"])
            estimateSizeFactors(cds)
            estimateDispersions(
                cds,
                method=method,
                sharingMode="fit-only" if method == "blind" else "maximum",
            )
            self.assertEqual(cds.dispersion_method, method)
            self.assertTrue(all(value >= 0 for value in cds.dispersions))

    def test_all_dispersion_fit_and_sharing_modes_are_supported(self):
        counts = [[10, 12, 9, 11], [100, 130, 20, 25], [5, 8, 40, 55]]
        for fit_type in ("local", "parametric", "mean"):
            for sharing_mode in ("maximum", "fit-only", "gene-est-only", "shrinkage"):
                with self.subTest(fit_type=fit_type, sharing_mode=sharing_mode):
                    cds = newCountDataSet(counts, ["T", "T", "C", "C"])
                    estimateSizeFactors(cds)
                    estimateDispersions(cds, fitType=fit_type, sharingMode=sharing_mode)
                    self.assertEqual(len(cds.dispersions), len(counts))
                    self.assertTrue(all(value >= 0 for value in cds.dispersions))

    def test_shrinkage_dispersion_lies_between_gene_and_fitted_estimates(self):
        cds = newCountDataSet(
            [[10, 12, 9, 11], [100, 130, 20, 25], [5, 8, 40, 55]],
            ["T", "T", "C", "C"],
        )
        estimateSizeFactors(cds)
        estimateDispersions(cds, sharingMode="shrinkage")
        for raw, fitted, final in zip(
            cds.raw_dispersions, cds.fitted_dispersions, cds.dispersions
        ):
            if raw is None or raw <= 0:
                self.assertEqual(final, fitted)
            else:
                self.assertGreaterEqual(final, min(raw, fitted))
                self.assertLessEqual(final, max(raw, fitted))

    def test_dispersion_options_are_validated(self):
        cds = newCountDataSet([[1, 2], [2, 3]], ["A", "B"])
        estimateSizeFactors(cds)
        for keyword, value in (
            ("method", "unknown"),
            ("sharingMode", "unknown"),
            ("fitType", "unknown"),
        ):
            with self.subTest(keyword=keyword):
                with self.assertRaisesRegex(ValueError, "unsupported"):
                    estimateDispersions(cds, **{keyword: value})

    def test_nbinom_test_reports_effect_direction_and_adjusted_pvalues(self):
        cds = newCountDataSet(
            [[100, 120, 10, 12], [30, 28, 31, 29], [3, 4, 70, 80]],
            ["T", "T", "C", "C"],
            row_names=["up", "same", "down"],
        )
        estimateSizeFactors(cds)
        estimateDispersions(cds, method="per-condition")

        rows = nbinomTest(cds, "C", "T")

        self.assertGreater(rows[0]["log2FoldChange"], 0)
        self.assertAlmostEqual(rows[1]["log2FoldChange"], 0)
        self.assertLess(rows[2]["log2FoldChange"], 0)
        self.assertTrue(all(0 <= row["padj"] <= 1 for row in rows))

    def test_nbinom_test_requires_fitted_dispersions_and_known_conditions(self):
        cds = newCountDataSet([[10, 12, 8, 9]], ["T", "T", "C", "C"])
        estimateSizeFactors(cds)
        with self.assertRaisesRegex(ValueError, "estimateDispersions"):
            nbinomTest(cds, "C", "T")
        estimateDispersions(cds)
        with self.assertRaisesRegex(ValueError, "conditions"):
            nbinomTest(cds, "missing", "T")

    def test_count_test_covers_exact_and_large_count_paths(self):
        cases = (
            (7, 20, 0.5, 0.0),
            (7, 20, 0.5, 0.2),
            (7500, 20000, 0.5, 0.0),
            (7500, 20000, 0.5, 0.2),
            (0, 0, 0.5, 0.0),
        )
        for arguments in cases:
            with self.subTest(arguments=arguments):
                pvalue = _two_sided_count_test(*arguments)
                self.assertGreaterEqual(pvalue, 0.0)
                self.assertLessEqual(pvalue, 1.0)

    def test_optimized_exact_count_test_matches_direct_enumeration(self):
        cases = (
            (7, 20, 0.3, 0.0),
            (12, 35, 0.6, 0.15),
            (2, 75, 0.1, 0.4),
        )
        for observed, total, probability, dispersion in cases:
            with self.subTest(
                observed=observed,
                total=total,
                probability=probability,
                dispersion=dispersion,
            ):
                if dispersion > 1e-10:
                    concentration = 1.0 / dispersion
                    alpha = probability * concentration
                    beta = (1.0 - probability) * concentration
                    log_pmf = lambda value: _log_beta_binomial_pmf(
                        value, total, alpha, beta
                    )
                else:
                    log_pmf = lambda value: _log_binomial_pmf(
                        value, total, probability
                    )
                observed_log_probability = log_pmf(observed)
                expected = sum(
                    math.exp(log_pmf(value))
                    for value in range(total + 1)
                    if log_pmf(value) <= observed_log_probability + 1e-12
                )
                self.assertAlmostEqual(
                    _two_sided_count_test(
                        observed, total, probability, dispersion
                    ),
                    min(expected, 1.0),
                    places=10,
                )

    def test_dispersed_count_test_remains_exact_above_old_cutoff(self):
        observed, total, probability, dispersion = 1300, 10001, 0.1, 0.5
        concentration = 1.0 / dispersion
        alpha = probability * concentration
        beta = (1.0 - probability) * concentration
        log_pmf = lambda value: _log_beta_binomial_pmf(value, total, alpha, beta)
        observed_log_probability = log_pmf(observed)
        expected = sum(
            math.exp(log_pmf(value))
            for value in range(total + 1)
            if log_pmf(value) <= observed_log_probability + 1e-12
        )
        self.assertAlmostEqual(
            _two_sided_count_test(observed, total, probability, dispersion),
            min(expected, 1.0),
            places=10,
        )

    def test_nbinom_test_is_invariant_to_sample_order_and_condition_reversal(self):
        counts = [[100, 120, 10, 12], [30, 28, 31, 29], [3, 4, 70, 80]]

        def analyze(matrix, conditions, numerator="T", denominator="C"):
            cds = newCountDataSet(matrix, conditions, row_names=["up", "same", "down"])
            estimateSizeFactors(cds)
            estimateDispersions(cds, method="per-condition")
            return nbinomTest(cds, denominator, numerator)

        baseline = analyze(counts, ["T", "T", "C", "C"])
        order = (1, 0, 3, 2)
        permuted = analyze(
            [[row[index] for index in order] for row in counts],
            ["T", "T", "C", "C"],
        )
        reversed_rows = analyze(counts, ["T", "T", "C", "C"], "C", "T")
        for original, shuffled, reversed_row in zip(baseline, permuted, reversed_rows):
            self.assertAlmostEqual(original["pval"], shuffled["pval"], places=12)
            self.assertAlmostEqual(original["padj"], shuffled["padj"], places=12)
            self.assertAlmostEqual(original["log2FoldChange"],
                                   -reversed_row["log2FoldChange"], places=12)
            self.assertAlmostEqual(original["pval"], reversed_row["pval"], places=12)

    def test_default_analysis_preserves_deseq2_table_contract(self):
        with tempfile.TemporaryDirectory() as directory:
            count_table = os.path.join(directory, "input.cntTable")
            project = os.path.join(directory, "native")
            self._write_count_table(count_table)

            analysis, significant = run_differential_analysis(
                count_table,
                2,
                2,
                project,
                padj_cutoff=1.0,
                fold_change_cutoff=1.0,
            )

            with open(analysis) as handle:
                lines = handle.read().splitlines()
            self.assertEqual(
                lines[0],
                "\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj",
            )
            self.assertEqual(len(lines), 4)
            self.assertTrue(os.path.isfile(significant))

    def test_legacy_modes_preserve_table_columns(self):
        with tempfile.TemporaryDirectory() as directory:
            count_table = os.path.join(directory, "input.cntTable")
            self._write_count_table(count_table)

            legacy = os.path.join(directory, "legacy")
            analysis, _ = run_differential_analysis(
                count_table, 2, 2, legacy, legacy_deseq=True
            )
            with open(analysis) as handle:
                header = handle.readline().rstrip("\n")
            self.assertEqual(
                header,
                "id\tbaseMean\tbaseMeanA\tbaseMeanB\tfoldChange\tlog2FoldChange\tpval\tpadj",
            )

            quantile = os.path.join(directory, "quantile")
            analysis, _ = run_differential_analysis(
                count_table,
                2,
                2,
                quantile,
                legacy_deseq=True,
                normalization="quant",
            )
            with open(analysis) as handle:
                header = handle.readline().rstrip("\n")
            self.assertEqual(
                header,
                "id\tsample1Mean\tsample2Mean\tFoldChange\tlog2FoldChange\tpval\tpadj",
            )

    def test_single_replicate_and_unbalanced_designs_write_results(self):
        with tempfile.TemporaryDirectory() as directory:
            single_table = os.path.join(directory, "single.cntTable")
            with open(single_table, "w") as handle:
                handle.write("gene/TE\tt1.T\tc1.C\n")
                handle.write("up\t100\t10\n")
                handle.write("down\t5\t80\n")
            analysis, significant = run_differential_analysis(
                single_table,
                1,
                1,
                os.path.join(directory, "single"),
                padj_cutoff=1.01,
            )
            self.assertTrue(os.path.isfile(analysis))
            self.assertTrue(os.path.isfile(significant))

            unbalanced_table = os.path.join(directory, "unbalanced.cntTable")
            with open(unbalanced_table, "w") as handle:
                handle.write("gene/TE\tt1.T\tt2.T\tc1.C\n")
                handle.write("up\t100\t120\t10\n")
                handle.write("down\t5\t4\t80\n")
            analysis, _ = run_differential_analysis(
                unbalanced_table,
                2,
                1,
                os.path.join(directory, "unbalanced"),
            )
            with open(analysis) as handle:
                self.assertEqual(len(handle.read().splitlines()), 3)

    def test_tc_normalization_accepts_valid_library_sizes(self):
        with tempfile.TemporaryDirectory() as directory:
            count_table = os.path.join(directory, "input.cntTable")
            self._write_count_table(count_table)
            analysis, _ = run_differential_analysis(
                count_table,
                2,
                2,
                os.path.join(directory, "tc"),
                legacy_deseq=True,
                normalization="TC",
                library_sizes=[100, 200, 100, 400],
            )
            self.assertTrue(os.path.isfile(analysis))

    def test_tc_normalization_rejects_invalid_library_sizes(self):
        with tempfile.TemporaryDirectory() as directory:
            count_table = os.path.join(directory, "input.cntTable")
            self._write_count_table(count_table)
            for sizes, message in ((None, "one library size"), ([1, 2], "one library size"), ([1, 2, 0, 4], "positive")):
                with self.subTest(sizes=sizes):
                    with self.assertRaisesRegex(ValueError, message):
                        run_differential_analysis(
                            count_table,
                            2,
                            2,
                            os.path.join(directory, "tc"),
                            legacy_deseq=True,
                            normalization="TC",
                            library_sizes=sizes,
                        )

    def test_count_table_validation_reports_actionable_errors(self):
        cases = (
            ("", "empty"),
            ("gene/TE\tt1.T\nfeature\t1\n", "samples"),
            ("gene/TE\tt1.T\tc1.C\nfeature\t1\n", "columns"),
            ("gene/TE\tt1.T\tc1.C\nfeature\tnope\t2\n", "non-numeric"),
            ("gene/TE\tt1.T\tc1.C\nfeature\t1\t1\n", "no features remain"),
        )
        with tempfile.TemporaryDirectory() as directory:
            for index, (contents, message) in enumerate(cases):
                filename = os.path.join(directory, "invalid-%d.cntTable" % index)
                with open(filename, "w") as handle:
                    handle.write(contents)
                with self.subTest(message=message):
                    with self.assertRaisesRegex(ValueError, message):
                        run_differential_analysis(
                            filename,
                            1,
                            1,
                            os.path.join(directory, "invalid-%d" % index),
                        )

    def test_public_api_rejects_all_zero_tables_and_sample_libraries(self):
        with tempfile.TemporaryDirectory() as directory:
            all_zero = os.path.join(directory, "all-zero.cntTable")
            with open(all_zero, "w") as handle:
                handle.write("gene/TE\tt1.T\tc1.C\n")
                handle.write("zero_a\t0\t0\n")
                handle.write("zero_b\t0\t0\n")
            with self.assertRaisesRegex(ValueError, "no features remain"):
                run_differential_analysis(
                    all_zero, 1, 1, os.path.join(directory, "all-zero")
                )

            zero_library = os.path.join(directory, "zero-library.cntTable")
            with open(zero_library, "w") as handle:
                handle.write("gene/TE\tempty.T\tt2.T\tc1.C\n")
                handle.write("one\t0\t10\t8\n")
                handle.write("two\t0\t30\t25\n")
            with self.assertRaisesRegex(ValueError, "empty.T.*only zero counts"):
                run_differential_analysis(
                    zero_library,
                    2,
                    1,
                    os.path.join(directory, "zero-library"),
                )

    def test_public_api_rejects_zero_sized_groups(self):
        with tempfile.TemporaryDirectory() as directory:
            count_table = os.path.join(directory, "input.cntTable")
            self._write_count_table(count_table)
            with self.assertRaisesRegex(ValueError, "treatment_count must be positive"):
                run_differential_analysis(
                    count_table, 0, 2, os.path.join(directory, "no-treatment")
                )
            with self.assertRaisesRegex(ValueError, "control_count must be positive"):
                run_differential_analysis(
                    count_table, 2, 0, os.path.join(directory, "no-control")
                )

    def test_significant_output_applies_fold_change_threshold(self):
        with tempfile.TemporaryDirectory() as directory:
            count_table = os.path.join(directory, "input.cntTable")
            self._write_count_table(count_table)
            _, significant = run_differential_analysis(
                count_table,
                2,
                2,
                os.path.join(directory, "filtered"),
                padj_cutoff=1.01,
                fold_change_cutoff=1000,
            )
            with open(significant) as handle:
                self.assertEqual(len(handle.read().splitlines()), 1)

    @staticmethod
    def _write_count_table(filename):
        with open(filename, "w") as handle:
            handle.write("gene/TE\tt1.T\tt2.T\tc1.C\tc2.C\n")
            handle.write("up\t100\t120\t10\t12\n")
            handle.write("same\t30\t28\t31\t29\n")
            handle.write("down\t3\t4\t70\t80\n")


if __name__ == "__main__":
    unittest.main()

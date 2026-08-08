from __future__ import division

import os
import tempfile
import unittest

from TEToolkit.DifferentialAnalysis import (
    estimateDispersions,
    estimateSizeFactors,
    nbinomTest,
    newCountDataSet,
    run_differential_analysis,
)


class DifferentialAnalysisTests(unittest.TestCase):
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

    @staticmethod
    def _write_count_table(filename):
        with open(filename, "w") as handle:
            handle.write("gene/TE\tt1.T\tt2.T\tc1.C\tc2.C\n")
            handle.write("up\t100\t120\t10\t12\n")
            handle.write("same\t30\t28\t31\t29\n")
            handle.write("down\t3\t4\t70\t80\n")


if __name__ == "__main__":
    unittest.main()

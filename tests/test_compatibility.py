import unittest

from tests.compatibility.compare_legacy_and_native import (
    RESULT_COLUMNS,
    _render_report,
)


class CompatibilityGateTests(unittest.TestCase):
    def test_adjusted_rank_and_significant_overlap_are_acceptance_gates(self):
        legacy = {
            "one": {"log2FoldChange": 2.0, "padj": 0.001},
            "two": {"log2FoldChange": 1.5, "padj": 0.9},
        }
        native = {
            "one": {"log2FoldChange": 2.0, "padj": 0.9},
            "two": {"log2FoldChange": 1.5, "padj": 0.001},
        }
        for rows in (legacy, native):
            for row in rows.values():
                row.update(
                    baseMean=10.0,
                    lfcSE=0.1,
                    stat=2.0,
                    pvalue=row["padj"],
                )

        report, failures = _render_report(
            RESULT_COLUMNS,
            legacy,
            RESULT_COLUMNS,
            native,
            "legacy-commit",
            "R test",
            "fixture.cntTable",
            0.85,
            0.60,
        )

        self.assertIn("adjusted-p-value rank correlation is below 0.85", failures)
        self.assertIn("significant-call Jaccard overlap is below 0.60", failures)
        self.assertIn("Adjusted-p-value rank correlation | -1.000", report)
        self.assertIn("Jaccard 0.000", report)


if __name__ == "__main__":
    unittest.main()

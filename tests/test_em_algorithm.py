import contextlib
import io
import unittest

from TEToolkit.EMAlgorithm import EMestimate, computeAbundances


class FakeFeatures(object):
    def __init__(self, lengths):
        self.lengths = lengths

    def getLength(self, index):
        return self.lengths[index]


class EMAlgorithmTests(unittest.TestCase):
    def test_abundance_assignment_conserves_every_multiread(self):
        with contextlib.redirect_stderr(io.StringIO()):
            counts = computeAbundances([0.75, 0.25, 0.0], [[0, 1], [1, 2], [0, 2]])
        self.assertAlmostEqual(sum(counts), 3.0)
        self.assertGreater(counts[0], counts[1])

    def test_zero_evidence_uses_symmetric_allocation(self):
        with contextlib.redirect_stderr(io.StringIO()):
            counts = computeAbundances([0.0, 0.0], [[0, 1]])
        self.assertEqual(counts, [0.5, 0.5])

    def test_em_uses_unique_evidence_and_is_permutation_symmetric(self):
        def estimate(unique):
            with contextlib.redirect_stderr(io.StringIO()):
                return EMestimate(
                    FakeFeatures([1000, 1000]), [[0, 1]], unique, [0.5, 0.5], 8, 50
                )

        forward = estimate([10, 1])
        reversed_counts = estimate([1, 10])
        self.assertAlmostEqual(sum(forward), 1.0)
        self.assertGreater(forward[0], forward[1])
        self.assertAlmostEqual(forward[0], reversed_counts[1], places=10)
        self.assertAlmostEqual(forward[1], reversed_counts[0], places=10)

    def test_em_handles_a_stationary_symmetric_solution(self):
        with contextlib.redirect_stderr(io.StringIO()):
            counts = EMestimate(
                FakeFeatures([1000, 1000]), [[0, 1]], [1, 1], [0.5, 0.5], 8, 50
            )
        self.assertEqual(counts, [0.5, 0.5])


if __name__ == "__main__":
    unittest.main()

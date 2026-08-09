import math
import unittest

from TEToolkit.Prob import (
    binomial_cdf,
    binomial_pdf,
    normal_cdf,
    normal_cdf_inv,
    poisson_cdf,
    poisson_pdf,
)


class ProbabilityTests(unittest.TestCase):
    def test_normal_cdf_matches_erf_and_inverse_round_trips(self):
        for value in (-3.0, -1.0, 0.0, 1.0, 3.0):
            expected = 0.5 * (1.0 + math.erf(value / math.sqrt(2.0)))
            self.assertAlmostEqual(normal_cdf(value), expected, places=7)
        for probability in (0.001, 0.1, 0.5, 0.9, 0.999):
            self.assertAlmostEqual(normal_cdf(normal_cdf_inv(probability)),
                                   probability, places=7)
        self.assertEqual(normal_cdf_inv(0.0), float("-inf"))
        self.assertEqual(normal_cdf_inv(1.0), float("inf"))

    def test_poisson_pdf_and_complementary_tails_match_direct_sum(self):
        rate = 2.5
        for observed in range(7):
            expected = math.exp(-rate) * rate ** observed / math.factorial(observed)
            self.assertAlmostEqual(poisson_pdf(observed, rate), expected, places=14)
            self.assertAlmostEqual(poisson_cdf(observed, rate) +
                                   poisson_cdf(observed, rate, lower=False), 1.0, places=14)

    def test_binomial_pdf_cdf_and_degenerate_boundaries(self):
        trials, probability = 8, 0.3
        running = 0.0
        for observed in range(trials + 1):
            expected = (math.comb(trials, observed) * probability ** observed *
                        (1.0 - probability) ** (trials - observed))
            running += expected
            self.assertAlmostEqual(binomial_pdf(observed, trials, probability),
                                   expected, places=14)
            self.assertAlmostEqual(binomial_cdf(observed, trials, probability),
                                   running, places=14)
        self.assertEqual(binomial_cdf(8, 8, 1.0), 1.0)
        self.assertEqual(binomial_cdf(8, 8, 1.0, lower=False), 0.0)


if __name__ == "__main__":
    unittest.main()

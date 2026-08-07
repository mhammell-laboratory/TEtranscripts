import os
import struct
import tempfile
import unittest

from TEToolkit.Normalization import _linear_scale_factor, _write_scatter_png


class NormalizationTests(unittest.TestCase):
    def test_zero_intercept_regression_matches_former_r_helper(self):
        self.assertAlmostEqual(
            _linear_scale_factor([2, 4, 6], [1, 2, 3]),
            2.0,
        )

    def test_zero_library_has_safe_identity_factor(self):
        self.assertEqual(_linear_scale_factor([1, 2], [0, 0]), 1.0)

    def test_scatter_plot_is_a_valid_500_pixel_png(self):
        with tempfile.TemporaryDirectory() as directory:
            filename = os.path.join(directory, "plot.png")
            _write_scatter_png(filename, [1, 10, 100], [2, 20, 50])
            with open(filename, "rb") as handle:
                signature = handle.read(8)
                chunk_length = struct.unpack(">I", handle.read(4))[0]
                chunk_name = handle.read(4)
                width, height = struct.unpack(">II", handle.read(8))
            self.assertEqual(signature, b"\x89PNG\r\n\x1a\n")
            self.assertEqual(chunk_length, 13)
            self.assertEqual(chunk_name, b"IHDR")
            self.assertEqual((width, height), (500, 500))


if __name__ == "__main__":
    unittest.main()

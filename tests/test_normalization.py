import os
import struct
import tempfile
import unittest
import zlib
from unittest import mock

from TEToolkit.Normalization import (
    _linear_scale_factor,
    _write_scatter_png,
    bin_corr,
    normalize,
    seq_depth,
)


class FakeSample(object):
    def __init__(self, name, library_size=1):
        self.name = name
        self.library_size = library_size
        self.read_species = []
        self.clear_count = 0

    def get_name(self):
        return self.name

    def libsize(self):
        return self.library_size

    def read_in_bins(self, species):
        self.read_species.append(species)

    def clear_bins(self):
        self.clear_count += 1


class NormalizationTests(unittest.TestCase):
    def _validate_png(self, filename, expected_size=(500, 500)):
        with open(filename, "rb") as handle:
            content = handle.read()
        self.assertEqual(content[:8], b"\x89PNG\r\n\x1a\n")

        chunks = []
        offset = 8
        while offset < len(content):
            self.assertGreaterEqual(len(content) - offset, 12)
            length = struct.unpack(">I", content[offset:offset + 4])[0]
            chunk_type = content[offset + 4:offset + 8]
            data_start = offset + 8
            data_end = data_start + length
            crc_end = data_end + 4
            self.assertLessEqual(crc_end, len(content))
            chunk_data = content[data_start:data_end]
            expected_crc = zlib.crc32(chunk_type + chunk_data) & 0xffffffff
            actual_crc = struct.unpack(">I", content[data_end:crc_end])[0]
            self.assertEqual(actual_crc, expected_crc, chunk_type)
            chunks.append((chunk_type, chunk_data))
            offset = crc_end

        self.assertEqual(offset, len(content))
        self.assertEqual(chunks[0][0], b"IHDR")
        self.assertEqual(chunks[-1], (b"IEND", b""))
        self.assertEqual(sum(kind == b"IHDR" for kind, _ in chunks), 1)
        self.assertEqual(sum(kind == b"IEND" for kind, _ in chunks), 1)

        width, height, bit_depth, color_type, compression, filtering, interlace = (
            struct.unpack(">IIBBBBB", chunks[0][1])
        )
        self.assertEqual((width, height), expected_size)
        self.assertEqual((bit_depth, color_type), (8, 2))
        self.assertEqual((compression, filtering, interlace), (0, 0, 0))

        compressed = b"".join(data for kind, data in chunks if kind == b"IDAT")
        self.assertTrue(compressed)
        scanlines = zlib.decompress(compressed)
        stride = width * 3 + 1
        self.assertEqual(len(scanlines), stride * height)
        self.assertTrue(all(scanlines[row * stride] == 0 for row in range(height)))
        pixels = b"".join(
            scanlines[row * stride + 1:(row + 1) * stride]
            for row in range(height)
        )
        self.assertTrue(any(value != 255 for value in pixels))

    def test_zero_intercept_regression_matches_former_r_helper(self):
        self.assertAlmostEqual(
            _linear_scale_factor([2, 4, 6], [1, 2, 3]),
            2.0,
        )

    def test_zero_library_has_safe_identity_factor(self):
        self.assertEqual(_linear_scale_factor([1, 2], [0, 0]), 1.0)

    def test_sequence_depth_scales_every_sample_to_largest_library(self):
        treatment = [FakeSample("t1", 25), FakeSample("t2", 50)]
        treatment_input = [FakeSample("ti", 10)]
        control = [FakeSample("c1", 100)]
        control_input = [FakeSample("ci", 20)]

        treatment_factors, control_factors = seq_depth(
            treatment, treatment_input, control, control_input
        )

        self.assertEqual(treatment_factors, [4.0, 2.0, 10.0])
        self.assertEqual(control_factors, [1.0, 5.0])

    def test_normalize_routes_supported_methods(self):
        arguments = ([], [], [], [], "hg38", "project")
        with mock.patch("TEToolkit.Normalization.seq_depth", return_value=([], [])) as depth:
            self.assertEqual(normalize("sd", *arguments), ([], []))
            depth.assert_called_once_with([], [], [], [])
        with mock.patch("TEToolkit.Normalization.bin_corr", return_value=([1], [2])) as corr:
            self.assertEqual(normalize("bc", *arguments), ([1], [2]))
            corr.assert_called_once_with(*arguments)

    def test_bin_correlation_scales_samples_writes_plots_and_clears_bins(self):
        treatment = [FakeSample("t1"), FakeSample("t2")]
        treatment_input = [FakeSample("ti")]
        control = [FakeSample("c1")]
        control_input = [FakeSample("ci")]
        reads = [
            [1, 2, 3],
            [2, 4, 6],
            [1, 1, 1],
            [3, 6, 9],
            [0, 0, 0],
        ]

        with tempfile.TemporaryDirectory() as directory:
            project = os.path.join(directory, "native")
            with mock.patch(
                "TEToolkit.Normalization.join_bins",
                return_value=([], reads, []),
            ) as join:
                treatment_factors, control_factors = bin_corr(
                    treatment,
                    treatment_input,
                    control,
                    control_input,
                    "hg38",
                    project,
                )

            self.assertEqual(treatment_factors, [1.0, 0.5, 2.0])
            self.assertAlmostEqual(control_factors[0], 1.0 / 3.0)
            self.assertEqual(control_factors[1], 1.0)
            join.assert_called_once_with(
                treatment,
                treatment_input,
                control,
                control_input,
                ["t1", "t2", "ti", "c1", "ci"],
                0,
            )
            for sample in treatment + treatment_input + control + control_input:
                self.assertEqual(sample.read_species, ["hg38"])
                self.assertEqual(sample.clear_count, 1)
            for comparison in ("t1vst2", "t1vsti", "t1vsc1", "t1vsci"):
                filename = project + "_" + comparison + ".png"
                self.assertTrue(os.path.isfile(filename), filename)

    def test_scatter_plot_is_a_valid_500_pixel_png(self):
        with tempfile.TemporaryDirectory() as directory:
            filename = os.path.join(directory, "plot.png")
            _write_scatter_png(filename, [1, 10, 100], [2, 20, 50])
            self._validate_png(filename)

    def test_scatter_plot_handles_zero_only_data(self):
        with tempfile.TemporaryDirectory() as directory:
            filename = os.path.join(directory, "empty-plot.png")
            _write_scatter_png(filename, [0, 0], [0, 0])
            self._validate_png(filename)


if __name__ == "__main__":
    unittest.main()

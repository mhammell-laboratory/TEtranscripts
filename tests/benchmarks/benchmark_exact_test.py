#!/usr/bin/env python3
"""Measure the exact count-test path at representative enumeration totals."""

import argparse
import os
import statistics
import sys
import time

REPOSITORY_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if REPOSITORY_ROOT not in sys.path:
    sys.path.insert(0, REPOSITORY_ROOT)

from TEToolkit.DifferentialAnalysis import _two_sided_count_test


TOTALS = (100, 1000, 5000, 10000)
DISPERSIONS = (0.0, 0.1)


def _measure(total, dispersion, features, repeats):
    durations = []
    checksum = 0.0
    for repeat in range(repeats):
        started = time.perf_counter()
        for feature in range(features):
            probability = 0.25 + (feature % 11) * 0.045
            observed = int(total * (0.15 + (feature % 13) * 0.055))
            observed = min(observed, total)
            checksum += _two_sided_count_test(
                observed, total, probability, dispersion
            )
        durations.append(time.perf_counter() - started)
    return statistics.median(durations), checksum


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", type=int, default=50)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument(
        "--assert-max-seconds",
        type=float,
        help="fail if either 10,000-total workload exceeds this duration",
    )
    args = parser.parse_args()
    if args.features <= 0 or args.repeats <= 0:
        parser.error("--features and --repeats must be positive")

    print("| Total count | Dispersion | Features | Median seconds | Features/second |")
    print("| ---: | ---: | ---: | ---: | ---: |")
    cutoff_durations = []
    checksum = 0.0
    for total in TOTALS:
        for dispersion in DISPERSIONS:
            duration, result_checksum = _measure(
                total, dispersion, args.features, args.repeats
            )
            checksum += result_checksum
            rate = args.features / duration if duration else float("inf")
            print(
                "| %d | %.1f | %d | %.4f | %.1f |"
                % (total, dispersion, args.features, duration, rate)
            )
            if total == 10000:
                cutoff_durations.append(duration)

    print("checksum: %.12g" % checksum)
    if (
        args.assert_max_seconds is not None
        and max(cutoff_durations) > args.assert_max_seconds
    ):
        raise SystemExit(
            "10,000-count exact workload exceeded %.3f seconds"
            % args.assert_max_seconds
        )


if __name__ == "__main__":
    main()

# Exact-test benchmark

`benchmark_exact_test.py` exercises the exact two-sided binomial and
beta-binomial code used by the Python-native legacy DESeq path. The exact path
enumerates totals through 10,000; larger totals use the continuity-corrected
normal approximation.

The implementation advances from `P(k)` to `P(k + 1)` with the distribution's
probability recurrence. Runtime therefore scales linearly with the total count
and number of tested features, while avoiding the repeated `lgamma` calls that
previously dominated work near the cutoff. Memory use is constant.

Run the reproducible workload with:

```console
python tests/benchmarks/benchmark_exact_test.py --features 50 --repeats 3
```

CI runs a smaller smoke benchmark and allows five seconds for each 25-feature
workload at the exact 10,000-count boundary. This generous limit is a regression
guard rather than a cross-run performance comparison; record measurements only
on otherwise idle, identified hardware when comparing revisions.

# Paired core benchmark

Baseline `a152e29cd4c4657e2f104aebe1f91319c49bbef8` and candidate source hashes are recorded in
[raw results](results/core_v0.2.1.json), including all samples, validation,
settings, dependency versions and profiles. macOS arm64, CPython 3.13.3, one
worker and one BLAS/OpenMP thread. 7 warm repetitions per case.

| Case | Baseline median s | Candidate median s | Ratio | Python peak B (old/new) | Process RSS B (old/new) |
|---|---:|---:|---:|---:|---:|
| small | 0.000123458 | 4.4041e-05 | 2.80x | 9952 / 6448 | 74498048 / 74579968 |
| medium | 0.00270929 | 0.00121871 | 2.22x | 66800 / 30096 | 75890688 / 75661312 |
| large | 0.0107519 | 0.00463033 | 2.32x | 166984 / 73856 | 77922304 / 77627392 |

Profiles identified temporary derivative polynomials, products and Fraction operations. Direct monomial-pair contraction avoids those temporaries, using exact integer arithmetic internally when a coefficient is integral. Full coefficients are compared with independent SymPy derivatives. Public polynomial coefficients remain rational.

Setup/import, first call, warm calls, separate allocation and separate profiling
passes are distinguished. Validation is outside timed regions; process RSS
includes imports, input and validation. Python tracemalloc excludes untraced
native allocations. These finite cases do not prove an asymptotic or universal speed/memory gain.
Timing noise and hardware matter; small RSS changes are not interpreted as gains.
The remaining arithmetic/BLAS kernels dominate larger workloads.

```sh
git worktree add --detach ../mechanical-lie-baseline a152e29cd4c4657e2f104aebe1f91319c49bbef8
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python bench/benchmark_core.py \
  --repeats 7 --baseline-source ../mechanical-lie-baseline/src \
  --output /tmp/mechanical-lie-core.json
```

The output path must be new. Install the package's development dependencies first;
SVD also needs Boost >=1.89 and C++17. Each variant imports the specified source
in a fresh worker. Performance thresholds are not enforced on shared CI.

Earlier policy benchmarks: [original guide](README.previous.md).

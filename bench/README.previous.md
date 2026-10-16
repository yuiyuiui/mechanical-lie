# Policy comparisons

Install the package and the independent test reference (`pip install -e
".[test]"`), then run from a checkout:

```sh
python bench/benchmark_policies.py --output numerical_case/results/bench --repeats 3
python bench/check_results.py numerical_case/results/bench/results.json
```

`--output` is a required new directory (`Path`/path string), never overwritten.
`--repeats` is an integer from 1 to 20, default 3. It controls measurement
replicates, not mathematical accuracy. Three is a bounded local illustration;
use more replicates for timing analysis after checking resource needs. CI uses
one sample only to check reproducibility and exact outputs. Cases, dimensions,
fields and reconstruction targets are fixed in the script. Native/BLAS thread
caps are set to one before reference-library imports and recorded. There is one worker,
no random sampling, and no seed. SymPy 1.14.0 is the optional exact reference;
the installed runtime library still has no external dependency.

| Matched inputs | Default and explicit choices | Output validation and work |
| --- | --- | --- |
| Three coupled quartic potentials, d=1,2,3; the same nested bracket of order 5 | Standard assignment cap, exact minimum admitting cap, and cap 200000 | Complete coefficients equal direct Poisson differentiation and independent SymPy; exact assignment counts, output monomials and hashes recorded |
| Unit triangular integer matrices of sizes 6 and 12 | Standard field and explicit F_2/F_101 | Full rank agrees with an exact rational reference. These unimodular inputs have the same full-rank answer in each field; arbitrary matrices need not. Input nonzeros and field settings recorded |
| The same finite reconstruction through order 6 | Standard field and F_101; automatic marker size 4 in each | All 16 integer columns are identical, all declared checks run, and every generated certificate is separately replayed |

A separate one-column control `[1000003]` has rank zero in the default field
and rank one over the rationals and F_101. It is retained to demonstrate why
modular deficiency is not rational dependence. No automatic field replacement
or skipped validation conceals this behavior.

Each policy warms once. Timed trials rotate the policy order in the same
process on identical inputs. Instrumented call time includes parameter
validation, selection/planning and execution; independent reference construction,
imports and final answer comparisons are excluded. Common polynomial input setup
is recorded separately. The reconstruction pipeline also includes temporary-file
setup, certificate replay and cleanup. Python allocation peaks are measured with
`tracemalloc`; these exclude preexisting caches and are not whole-process RSS.
Timing is affected by instrumentation, shared caches and background load. The
small measurement is not a performance ranking or evidence of an optimal default.

The [published compact result](results/reference.json) contains the actual
environment, repeats, selected configurations/reasons, exact checks, work/storage
counts, trial times and allocation peaks. Source hashes identify the exact
implementation and benchmark script, including a pre-commit source checkout when
the observed Git commit precedes the release. This avoids treating a Git hash as
evidence for uncommitted changes. Final published source bytes must match those
hashes. Reproduction on another platform is expected to preserve exact outputs
and work counts, not timing/allocation bytes.

The caps only change whether work is admitted. All admitted settings compute the
same complete polynomial. A field change can change modular rank. No floating
tolerance, accuracy mode, probabilistic guarantee, all-order theorem, or
independent human/formal certification follows from these benchmarks.

# Fixed symbolic cases and figures

The release uses exact integers and rational coefficients. It contains no
runtime/error comparison with an ODE solver and no claimed speed advantage.

| Case | Actual inputs | Exact checks | Published material |
| --- | --- | --- | --- |
| Finite reconstruction | Original-letter orders 1–10, prime 1000003, comb marker with 6 leaves | 104 source columns; 31 matrix blocks reconstructed and replayed | `sparsity.csv`, `matrix_pattern.json`, `checks.json` |
| Full polynomial realization | Orders 1–6, dimensions 1 and 2, two fixed quartic potentials | 32 complete coefficientwise comparisons with direct recursive Poisson differentiation | Exact potentials and counts in `checks.json` |
| Tree example | d=1, V=q²/2+q⁴/4; V, {V,T}, {V,{V,T}} | Full rational polynomial outputs | `tree_examples.json`, bracket-trees figure |
| Historical dimensions | Previously computed orders 1–18, prime 1000003 | Original finite full ranks; fresh orders 1–10 checked for agreement | `graded_dimensions.csv`, attributed in `provenance.json` |

The historical scalar table is reused byte-for-byte from the
[fixed original public version](https://raw.githubusercontent.com/yuiyuiui/AIM/8dc12e3f97d7c2230981ad0d46ab00aa06188bfd/research/solutions/519-euclidean-mechanical-lie-algebra/numerical_case/summaries/rank_summary.csv).
Its SHA-256 is `341be7eda0e4640fbf5a89e73b600c8b7a03e3a835472a526b4d89781d527dfb`.
This release does **not** rerun order 18. The figure uses the actual order,
positive dimension and total dimension columns, rather than a hard-coded title
or a guessed prime. The associated manuscript supplies the mathematical
argument separately from these finite cases.

The two potentials for complete coefficient checks are q1²/2+q1⁴/4 and
q1²/2+q1*q2+q2⁴/4. These checks differ from the earlier manuscript's two
degree-five pointwise calibrations. They are new finite software checks, not
an expanded mathematical claim.

## Reproduce

From an installed checkout, run the README commands with a new output
directory. `compute_cases.py` accepts only `--output`: all scientific inputs
are fixed in the script and accurately recorded in its manifest. Edit the
separate worked example to try another potential. The Git revision is optional
metadata; computing from a source distribution without a Git executable works.

The raw output contains integer matrices, pivots and the compact summaries.
Source/output SHA-256 values and actual parameters are recorded. Full matrices
are local and ignored; only the compact files under `summaries/` are shipped.
The public record `provenance.json` contains only relative software paths,
content hashes and the explicitly public historical-source URL.

`plot_cases.py` consumes the compact results and the historical scalar table,
and saves three vector SVG/PDF pairs. Optional `--preview-dir PATH` saves PNG
previews for visual inspection. `check_outputs.py` compares the four newly
computed compact summaries byte-for-byte and checks both figure formats.
Plots display exact dimensions, derivative/tree structure and actual support
counts; the only floating-point quantities are drawing coordinates and a
rounded matrix-density label.

The demo and CI are bounded at order 10, use one process and cap native/BLAS
threads to one for plotting. Recursive expansion, bigger potentials and higher
dimensions can become expensive. Consult the method limits before extending
an experiment. Replayed finite matrices do not establish an all-order theorem.

# API and parameter contract

The supported surface is `mechanical_lie.__all__` plus the listed
`mechanical_lie.algebra` and `mechanical_lie.validation` functions below.
`_engine` and underscore helpers are implementation details. Use an installed
package without adding `src` to `PYTHONPATH`.

Coordinates are formal q1,...,qd,p1,...,pd, with zero-based Python indices.
No physical units or scaling are inferred or converted. Use consistent units
and the Euclidean kinetic energy T=sum(p_i**2)/2; general mass tensors are not
supported. Exact coefficients/coordinates/scalars are Python `int` or
`fractions.Fraction`; booleans and floats are rejected. Polynomials and
Hamiltonians are immutable and `.terms` returns a copy; zero terms are removed.
There is no conditioning tolerance, accuracy level, approximate backend,
parallel worker setting or random seed.

## Inputs and outputs

All parameters in this table are required, with no default or inferred problem
size. Dimension, degree, expression size and rational bit lengths can increase
work and memory sharply. Exact arithmetic avoids roundoff, not growth.

| Callable / parameter | Type, range and mathematical role | Return and advice |
| --- | --- | --- |
| `Polynomial(dimension, terms)` | Positive integer d; mapping exponent tuples to exact coefficients. Each tuple has exactly 2d nonnegative integer entries, q then p. Empty mapping is zero. | Sparse exact polynomial. Choose d and degree from the actual problem; no truncation. |
| `phase_space(dimension)` | Positive integer d. | Two tuples of d position/momentum polynomials; recommended coordinate constructor. |
| `Polynomial.constant(dimension, value)` | Positive integer d and exact scalar value. | Constant polynomial in that phase space. |
| `Polynomial.variable(dimension, index)` | Positive integer d; integer index in range(2*d). | Coordinate polynomial: indices 0,...,d-1 are positions, d,...,2d-1 momenta. |
| `Polynomial.kinetic(dimension)` | Positive integer d. | Exact Euclidean kinetic polynomial. |
| `Polynomial.derivative(index)` | Integer index in range(2*d). | Exact partial derivative; index defines the question. |
| `Polynomial.evaluate(point)` | Sequence of exactly 2d exact coordinates, q then p. | `Fraction`; incorrect coordinate counts fail rather than being truncated. |
| `poisson(left, right)` | Two same-dimension `Polynomial` objects. | Exact sum(left_q*right_p-left_p*right_q); useful direct reference. |
| `Hamiltonian(terms)` | Mapping mechanical-tree codes to exact coefficients. `None` denotes T; a tree node is `(momentum_count, children_tuple)` with nonnegative integer count and recursively valid children. | Canonical expression; input rootings are canonicalized and equal terms combined. `(0,())` denotes V. |
| `kinetic()`, `potential()` | No parameters. | Universal T or one V; dimension is specified only at realization. |
| `bracket(left, right)` | Two Hamiltonians. | Exact universal Poisson bracket, with {V,T}=p dot grad V. |
| `commutative_product(left, right)` | Two Hamiltonians with only zero momentum degree (zero allowed). | {right,{left,T}}, the commutative, generally nonassociative gradient product. |
| `binary_trees(leaves_count)` | Positive integer leaf count. | Tuple of canonical nonplanar binary trees; enumeration grows rapidly and is not silently capped. |
| `potential_tree(tree)` | Canonical commutative binary tree: leaf `()` or sorted pair of valid trees. | Its gradient-product Hamiltonian. |

Polynomial +,-,* accept same-dimension polynomials or exact scalars. Division
requires a nonzero exact scalar; powers require a nonnegative integer exponent.
Hamiltonian +,- require Hamiltonians; scalar multiplication/division require
exact scalars, with nonzero divisor. Ordinary Hamiltonian multiplication is
unsupported. These operands/powers/divisors define the mathematical inputs.
`.orders` and `.momentum_degrees` are sorted distinct tuples, empty for zero.

## Complete realization and its resource guard

`h.realize(potential, *, max_assignments=100000) -> Polynomial` is the input-only
default path. `h.realization_plan(potential, *, max_assignments=100000)` exposes
its effective configuration without differentiating or enumerating assignments.

| Parameter | Contract and rule | Recommendation and accuracy/cost effect |
| --- | --- | --- |
| `potential` | Required position-only Polynomial; d comes from this input. Momentum dependence fails. | Defines V and d; nonpolynomial potentials are unavailable. |
| `max_assignments` | Positive Python integer, default 100000 per nonzero tree. The automatic work estimate is d**(E+F), E thick edges and F momentum leaves. | Keep the rounded engineering safeguard for small work. A smaller cap is a conservative admission policy; raising it is an explicit resource decision. It never improves accuracy. |

All trees are checked before constructing polynomial workspaces. A count equal
to the cap is admitted; a larger one raises `ValueError`. There is no partial,
truncated or approximate polynomial return. Every admitted cap gives the same
complete exact output. The cap does not bound total assignments across trees,
polynomial term growth, coefficient bits, canonicalization, total memory or
wall time. Planning must itself read the supplied expression. It is not a
hardware capacity check. No unchecked automatic enlargement occurs.

The frozen `RealizationPlan` contains `dimension`, `tree_assignments` (in
`.terms` iteration order), effective `max_assignments`, `backend` and
`selection_reason`. Derived properties are `largest_assignments`,
`total_assignments` and `within_limit`; `.to_dict()` returns a JSON-serializable
copy. A kinetic term counts as zero contraction assignments, not zero work.
An empty expression has no assignments. Constant potentials retain conservative
tree counts; no derivative-zero shortcut is promised. An over-cap plan is
returned for inspection, whereas `realize` rejects it.

```python
from mechanical_lie import bracket, kinetic, phase_space, potential
q, p = phase_space(1)
v = q[0]**2/2 + q[0]**4/4
h = bracket(potential(), bracket(potential(), kinetic()))
answer = h.realize(v)
print(h.realization_plan(v).to_dict())
assert h.realize(v, max_assignments=1000) == answer
```

## Finite-field rank

`modular_rank(columns, prime=1000003) -> int` preserves its original return type.
`modular_rank_report(columns, prime=1000003) -> ModularRankReport` performs the
same exact calculation and exposes the effective field.

| Parameter | Type/default/role | Selection advice and effects |
| --- | --- | --- |
| `columns` | Required finite iterable of mappings from string rows to integer coefficients. Fractions, floats and booleans are rejected. Missing rows mean zero; empty columns/iterables are valid. | Materialized once. Shape, fill-in and coefficient bits affect time/storage. |
| `prime` | Prime Python integer >=2; default 1000003. Validated before consuming columns. | An explicit prime defines the field. For rational-independence checks the standard moderate prime is an engineering choice, not an optimal or universally good prime. Small primes such as 2/101 are cheap controls. Large primes increase integer arithmetic and trial-division primality-check cost. No matrix-dependent/random prime search or automatic retry. |

The frozen report contains exact `rank`, `column_count`, nonzero integer input
`row_count`, actual `prime`, `backend='sparse-modular'`, `selection_reason`, and
`.to_dict()`. `rational_independence_established` is true exactly when rank equals
the number of columns, including the empty family. False means this calculation
has not established full rational independence; it does not mean dependence.
Changing the prime can change rank because the field changed, not roundoff or
conditioning. Input support may count coefficients that vanish modulo p.
For example, column {'x':2} has rank 1 at the default prime and rank 0 at prime 2.

## Finite reconstruction and replay

`validation.reconstruction_plan(max_order, *, prime=1000003)` returns a frozen
`ReconstructionPlan` without enumeration/filesystem writes.
`validation.rebuild_basis(max_order, output, *, prime=1000003)` uses the same
policy and returns a summary with `effective_configuration`; saved
`audit_summary.json` and returned summary agree.

| Parameter or automatic output | Contract | Recommendation and effect |
| --- | --- | --- |
| `max_order` | Required Python integer 2–18; finite original-letter target. | Start at 6 or below for inspection; the fixed demo uses 10. Changing order changes the question. Upper supported order 18 is not a resource guarantee and is not run in CI. |
| `output` | Required string/path-like value resolved by `Path`, naming a new directory. | Parents are created; existing output raises `FileExistsError`. Invalid policy values fail before creation. Disk stores exact matrices/pivots/summaries. |
| `prime` | Same prime/default as modular rank. | Deficiency raises `ValidationError` because this full-rank finite check did not pass, not because dependence over the rationals was proved. |
| `marker_leaves` (output only) | Automatically max_order//2+1, using a comb. | Strictly exceeds the source side-label leaf budget in this construction. No independent marker override; this is mathematical bookkeeping, not a hardware guess. |
| `verification` (output only) | Always `all-stages` in this API. | No fast mode that silently skips required rank/reconstruction. |

The plan records max_order, actual prime, marker_leaves, `prime_reason` and `marker_reason`;
`.to_dict()` gives a copy. The requested target is never reduced to fit resources.
Failure can leave partial files; retry in a fresh directory. An elapsed time,
hash or complete finite summary is not an all-order theorem.

`validation.verify_certificate_block(data) -> int` requires a mapping with
`order, degree, prime, dimension, rank, columns, pivots, source_basis,
sha256_without_hash`. Counts must be valid/consistent integers, columns must map
string rows to integers, the modulus must be prime, and every pivot/hash must
match. Corruption or a dependent block raises `ValidationError` (a `ValueError`
subclass). Replay checks a supplied matrix, not its symbolic origin.
`validation.digest(payload) -> str` takes a JSON-serializable mapping with
string keys and hashes canonical JSON, excluding `sha256_without_hash`.
Recomputing a digest does not validate a matrix or its provenance.

## Advanced exact words and markers

`algebra.A=('A',)` and `algebra.B=()` are constants. A tree is `()` or a sorted
pair of trees. A word is a tuple of A or Y(tree) letters; `()` is the associative
unit. A word polynomial maps words to exact int/Fraction coefficients. These
representations are validated at the public boundary. Structural caches may grow
for the process lifetime; mapping caches return copies (see methods).
Every parameter below is required, with no default. Leaf/word sizes are
combinatorial counts, not accuracy or hardware settings.

| Callable / parameters | Valid values, return and advice |
| --- | --- |
| `Y(t)`, `leaves(t)` | Canonical tree t; return its letter or positive leaf count. |
| `star(a,b)` | Two canonical trees; return their sorted pair. Prefer this constructor for pair ordering. |
| `comb(n)` | Positive integer n; return n-leaf comb. More leaves increase structural work. |
| `lyndon(w)` | Valid word w; return bool in fixed lexicographic order. Empty word is false. |
| `lie_expansion(w)` | Nonempty Lyndon word w; return standard associative expansion. Otherwise `ValueError`. Long words can have many terms. |
| `commutator(a,b)` | Two exact word polynomials; return a*b-b*a. Product support may grow multiplicatively. |
| `adjoint_x(x,f)` | Canonical marker x, exact word polynomial f; return normal PBW mapping `(word,potential_factor) -> coefficient`, with `None` meaning no factor. No homogeneous Lie purity is inferred for arbitrary f. |
| `derivative_A(f,x)` | Exact word polynomial f, canonical marker x; replace each A in turn by Y(x). |
| `project_large_labels(f,M)` | Exact word polynomial f; integer threshold M>=0. Delete words with a Y label having more than M leaves. M=0 keeps A words only. This exact projection changes the input, not numerical accuracy. |
| `retract_marker(f,x)` | Exact word polynomial f, canonical marker x; replace exactly Y(x) by A and combine terms. |

A faithful higher-degree marker argument requires marker size greater than each
old Y label size. Degree-one separation requires it greater than the sum of
side-label leaf counts. These are promises of that use, not conditions inferred
by generic word operations. The finite reconstruction's automatic rule must not
be transferred to arbitrary words or treated as a purity check.

## Errors and compatibility

Type/representation failures raise `TypeError`/`ValueError` with a parameter or
structural explanation; zero division raises `ZeroDivisionError`. Unavailable
backend keywords raise `TypeError`, without fallback. There is no supported
None/auto sentinel for method settings or hidden optional arithmetic backend.
SymPy and Matplotlib are test/figure extras, not runtime dependencies.

v0.2.0 preserves valid v0.1.0 calls, integer defaults, complete exact outputs,
Polynomial/int return types for realize/modular_rank. Diagnostic APIs and the
reconstruction summary field are additive. Invalid primes fail earlier, all
contraction guards precede workspaces, and malformed advanced trees/words now
fail at the boundary. Old releases/tags are not rewritten. See the
[policy tests](../tests/test_parameter_policy.py), [matched benchmarks](../bench/README.md)
and [methods/limits](methods.md).

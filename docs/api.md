# API reference

The stable import surface for v0.1.0 is `mechanical_lie.__all__`.
All coefficient inputs are Python `int` or `fractions.Fraction`; booleans and
floats are rejected. Returned term dictionaries are copies of immutable state.

## Concrete polynomials

`Polynomial(dimension, terms)` stores a finite mapping from nonnegative integer
exponent tuples to rational coefficients. Tuples have length `2*dimension` in
the order `(q1,...,qd,p1,...,pd)`. The empty mapping is zero. `dimension` is a
positive integer. Ordinary `+`, `-`, `*`, division by a nonzero exact scalar,
and nonnegative integer powers are supported. Mixed dimensions raise `ValueError`.

- `phase_space(d) -> (positions, momenta)`: two tuples of length d.
- `Polynomial.constant(d,c)`, `.variable(d,index)`, `.kinetic(d)`.
- `.terms -> dict`: a fresh exact sparse representation.
- `.derivative(index)`: index in `range(2*d)`.
- `.evaluate(point) -> Fraction`: exactly `2*d` exact coordinates.
- `poisson(f,g) -> Polynomial`: direct polynomial Poisson differentiation.

```python
from fractions import Fraction
from mechanical_lie import phase_space, poisson
q, p = phase_space(2)
v = q[0]**2 / 2 + q[0]*q[1] + q[1]**4 / 4
print(poisson(v, p[0]))
print(v.evaluate((Fraction(1,2), 1, 0, 0)))
```

## Universal mechanical expressions

`kinetic()` and `potential()` create universal Hamiltonians for `T=|p|²/2` and
`V(q)`. `bracket(f,g)` expands the mechanical Poisson bracket with the same sign
as the polynomial API. Linear arithmetic and exact scalar scaling return new
`Hamiltonian` objects; ordinary multiplication of Hamiltonians is unsupported.

`Hamiltonian(terms)` maps canonical mechanical trees to rational coefficients.
A tree is recursively `(momentum_count, children_tuple)`; `None` denotes T.
Input tree rootings are canonicalized. The singleton `(0,())` denotes V.
Every momentum count must be a nonnegative integer.

- `.terms`: a fresh mapping; equal canonical terms are combined.
- `.orders`: sorted distinct original-letter orders of its nonzero terms.
- `.momentum_degrees`: sorted distinct momentum degrees; empty for zero.
- `.realize(v, max_assignments=100000)`: a complete exact polynomial for a
  position-only Polynomial v. Momentum-dependent potentials raise `ValueError`.
  A tree with more than the cap's index assignments raises `ValueError` before
  enumerating those assignments. The cap is per tree and does not bound total
  polynomial expansion size or runtime.

`commutative_product(f,g)` is `{g,{f,T}}` for degree-zero inputs. It realizes
`∇f·∇g`, so it is commutative but generally nonassociative.
`binary_trees(n)` enumerates canonical nonplanar binary trees with n leaves;
`potential_tree(tree)` evaluates that binary-tree product. Leaf `()` is V.
The enumeration grows rapidly; there is no silent order truncation.

## Advanced associative-word operations

`mechanical_lie.algebra` exposes the already implemented low-level integer
operations. A commutative tree is `()` or a sorted pair of such trees;
`star(a,b)` constructs a sorted pair, and `comb(n)` gives an n-leaf comb.
Letters are `A=('A',)` and `Y(t)=('Y',t)`. A word is a tuple of these letters.
An associative polynomial is a dictionary `word -> integer coefficient`.

`lyndon(word)` tests the fixed lexicographic order; `lie_expansion(word)` applies
standard Lyndon bracketing and returns its associative expansion. It requires
a nonempty Lyndon word. `commutator(f,g)` uses associative multiplication.

`adjoint_x(marker,f)` computes `[marker,f]` in normal PBW form, returning a
mapping `(word,potential_factor) -> integer coefficient`. `None` means no
potential factor, while `word=()` is the empty word. The purity conclusions
used for degree lowering apply to homogeneous Lie polynomials, not arbitrary
associative inputs. `derivative_A` replaces one A at a time by Y(marker).
`project_large_labels(f,M)` deletes words with a Y label of more than M leaves.
`retract_marker(f,marker)` replaces exactly Y(marker) by A.

For a faithful higher-degree marker test, marker size must strictly exceed
each old Y label size. For degree-one word separation it must exceed the sum
of side-label leaf counts. These are explicit mathematical preconditions.
The advanced tuple routines assume well-formed input representations; the
validated object API above is the recommended starting point.

## Finite validation

`modular_rank(columns, prime=1000003) -> int` accepts integer sparse columns
with string row labels, rejects composite/noninteger moduli and inexact matrix
coefficients. Full column rank modulo a prime implies rational independence;
deficiency at that prime alone does not prove rational dependence.

`mechanical_lie.validation.rebuild_basis(max_order, output, prime=1000003)`
requires order 2–18, a prime, and a new output directory. It checks all declared
stages and saves integer matrices/pivots. The small shipped case uses order 10.
The returned summary explicitly lists checked/skipped stages. A skipped stage
is never an overall PASS. Order 18 can require substantial memory and is an
explicit, user-initiated operation, not an installation or CI test.

`verify_certificate_block(record) -> rank` checks a block's canonical JSON
digest, integer matrix, dimensions and every modular pivot. Malformed or corrupt
records raise `ValidationError` (a `ValueError` subclass). Replay checks the
supplied matrix; rebuilding checks its relationship to the bracket construction.

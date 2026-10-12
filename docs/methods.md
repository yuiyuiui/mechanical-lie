# Methods and limitations

The fixed kinetic energy is Euclidean, `T=Σp_i²/2`. The potential is a single
position-dependent scalar function; the concrete interface uses finite real
rational-coefficient polynomials. General mass tensors, multiple independent
potentials, nonpolynomial functions and numerical trajectory integration are
outside this version's contract.

## Representations

Commutative nonassociative potential products are represented by sorted binary
trees. Free Lie words use a weighted Lyndon basis and exact associative
commutators. The lowering code normal-orders a potential through generator
words and retains a possible trailing potential factor.

A mechanical tree has a thick connected tree of potential vertices and an
integer count of attached momentum leaves at each vertex. Each thick edge
contracts a derivative index at both endpoints; a momentum leaf contracts one
derivative with a momentum component. Kinetic energy is a separate sentinel.
Canonicalization minimizes over tree rootings and counts equivalent outputs
with their full multiplicity. The Hamiltonian wrapper exposes immutable
canonical terms without leaking mutable cache results.

The polynomial realization enumerates every derivative-index assignment and
multiplies the resulting sparse rational derivative polynomials. It is
independent of direct recursive Poisson differentiation; the tests compare full
coefficient dictionaries and also use SymPy as an external exact reference.
No floating-point tolerance, numerical differentiation or SVD rank cutoff is
involved. Converting an arbitrary decimal to a rational is a caller decision:
the public API deliberately rejects float inputs.

## Growth and memory

For dimension d, E thick edges and F momentum leaves, contraction enumerates
`d**(E+F)` assignments. The default per-tree cap is 100000. Sparse polynomial
products and exact numerator/denominator bit lengths can grow even below that
cap. Binary-tree and Lyndon enumerations also grow quickly. These algorithms
are intended for small-to-moderate symbolic questions, not an unlimited-cost
service or a high-dimensional solver.

Sparse-result mapping caches are bounded at 4096 entries per function and
return fresh dictionaries. Some small structural enumeration/encoding caches
are process-lifetime caches with no size cap. Run independent large campaigns
in separate processes; do not assume the mapping cap bounds total memory.
The examples and CI use fixed modest orders and one worker. Plot scripts cap
native/BLAS thread counts before importing the plotting stack.

## Evidence boundaries

The exact finite reconstruction checks source dimensions, momentum grading,
lowering identities, actual tree columns and full modular rank. Replaying a
digest and pivot record alone checks a supplied matrix, not its symbolic
origin. Changing a prime can lower modular rank without lowering rational rank.
Partial runs name skipped stages; neither a summary nor a hash is a universal
proof certificate.

The included historical order-18 scalar table is an attributed previously
computed dataset. The new release reconstructs only the documented small cases.
Complete polynomial equalities for chosen potentials do not establish a
joint-kernel theorem, faithfulness for all dimensions, or a fixed-potential
injectivity theorem. The associated mathematical manuscript and its review
status remain separate from software installation, CI and release success.

# Changelog

## 0.2.0 — 2026-10-15

- Preserve `realize`'s 100000 assignment guard, the 1000003 modular prime,
  required finite-order targets and existing result types.
- Add immutable realization, modular-rank and reconstruction diagnostics, with
  exact work counts, automatic marker selection and effective-policy reasons.
- Check all contraction guards before allocating polynomial workspaces, reject
  invalid primes before consuming matrix iterators, and validate advanced
  tree/word inputs at the public boundary. Valid exact-rational calls remain
  supported; errors for malformed inputs are now earlier and clearer.
- Document every public parameter and supported alternative, add output-level
  policy regressions and matched bounded exact-output benchmarks, and explicitly
  disclose AI-created implementation, API, tests, benchmarks and documentation.
- Preserve v0.1.0 and the associated mathematical submission unchanged.

## 0.1.0 — 2026-10-13

- Initial standalone exact symbolic package with canonical Hamiltonian trees,
  immutable rational polynomials and complete polynomial realization.
- Explicit Poisson sign, dimensional validation, exact coefficient handling
  and contraction-size guards.
- Graded tree/Lie-word and finite modular-rank/certificate operations.
- Modifiable quartic example, independent exact checks, three computed SVG/PDF
  figures, source summaries and reproducible case scripts.
- Wheel/sdist packaging, Python 3.11–3.14 CI, BSD 3-Clause license and citation.

Finite validation and release success do not change the mathematical
manuscript's independent-review status.

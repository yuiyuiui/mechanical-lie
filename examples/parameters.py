"""Use defaults first, inspect actual work, and choose explicit guards/fields."""

from mechanical_lie import bracket, kinetic, modular_rank_report, phase_space, potential
from mechanical_lie.validation import reconstruction_plan

q, _ = phase_space(2)
v = q[0] ** 2 / 2 + q[0] * q[1] + q[1] ** 4 / 4
h = bracket(potential(), bracket(potential(), kinetic()))
answer = h.realize(v)  # No method tuning is required.
print("Realization:", h.realization_plan(v).to_dict())
if h.realize(v, max_assignments=2) != answer:
    raise RuntimeError("An admitted resource guard changed exact coefficients")
print("Exact polynomial:", answer)

columns = [{"x": 2}]
print("Default field:", modular_rank_report(columns).to_dict())
print("Explicit field F_2:", modular_rank_report(columns, prime=2).to_dict())
print("Finite target (no computation):", reconstruction_plan(6).to_dict())

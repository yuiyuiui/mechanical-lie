# Development

Use Python 3.11–3.14 in a virtual environment and install `-e ".[dev]"`.
The exact library has no third-party runtime dependency. SymPy is an independent
test reference; Matplotlib is used only for reproducible scientific figures.
The development extras are pinned in pyproject.toml.

Run both normal and optimized unittest discovery, Ruff lint/format checks,
the quartic example and `python -m build`. For an installation check, create
a second clean environment, install the built wheel with no dependencies,
change out of the source checkout, and run:

```sh
python -c 'from mechanical_lie import kinetic,potential,bracket; print(bracket(potential(),kinetic()).momentum_degrees)'
```

The expected output is `(1,)`. Install the `test` extra separately to run the
full suite against that wheel. Do not put `src` on PYTHONPATH when claiming
an installed-package test. CI follows this pattern and also tests an sdist
build and the README's worked example.

New symbolic operations should include a hand calculation, sign/degree checks
and an independent exact reference. Preserve dimension and float-input errors.
Treat finite rank as finite evidence. Performance claims need matched cases,
accuracy/output contracts and measured repetitions; the current figures make
no timing or solver-accuracy claims. Regenerate compact case summaries and
figures when their calculations change. Keep large raw output out of Git.

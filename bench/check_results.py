"""Check the compact policy benchmark's exact-output and evidence contract."""

import argparse
import json
from pathlib import Path


def validate(path: Path):
    result = json.loads(path.read_text())
    if result["schema"] != "mechanical-lie-policy-benchmark-v1":
        raise ValueError("Unknown benchmark schema")
    groups = {}
    for row in result["records"]:
        if (
            row["validation"] != "PASS"
            or len(row["elapsed_seconds"]) != result["measurement"]["repeats"]
        ):
            raise ValueError("Incomplete or failed benchmark")
        groups.setdefault((row["category"], row["case"]), []).append(row)
    if len(groups) != 6:
        raise ValueError("Expected three realization, two rank and one reconstruction case")
    for (category, _), rows in groups.items():
        if "default" not in {r["configuration"] for r in rows} or len(rows) < 2:
            raise ValueError("Missing default/explicit comparison")
        if category == "realization":
            if len({r["exact_output_sha256"] for r in rows}) != 1 or not all(
                r["effective_configuration"]["within_limit"] for r in rows
            ):
                raise ValueError("Exact output invariance failed")
        elif category == "reconstruction":
            if len({r["integer_matrix_sha256"] for r in rows}) != 1:
                raise ValueError("Integer matrices changed with field policy")
        elif any(r["observed_rank"] != r["rational_reference_rank"] for r in rows):
            raise ValueError("Unimodular rank reference failed")
    control = result["bad_prime_control"]
    if (
        control["default"]["rational_independence_established"]
        or control["default"]["rank"] != 0
        or control["explicit_prime_101"]["rank"] != 1
    ):
        raise ValueError("Invalid bad-prime interpretation")
    print("Policy benchmark: exact outputs, matched comparisons and bad-prime control PASS.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    validate(parser.parse_args().results)

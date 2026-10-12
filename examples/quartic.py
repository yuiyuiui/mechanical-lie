"""A modifiable exact one-dimensional mechanical Hamiltonian example."""

from mechanical_lie import Polynomial, bracket, kinetic, phase_space, poisson, potential


def main():
    q, p = phase_space(1)
    v = q[0] ** 2 / 2 + q[0] ** 4 / 4  # Change the rational polynomial here.
    t = Polynomial.kinetic(1)
    a, b = kinetic(), potential()
    cases = [
        ("{V,T}", bracket(b, a), poisson(v, t)),
        ("{V,{V,T}}", bracket(b, bracket(b, a)), poisson(v, poisson(v, t))),
        ("{T,{T,V}}", bracket(a, bracket(a, b)), poisson(t, poisson(t, v))),
    ]
    for name, expression, direct in cases:
        contracted = expression.realize(v)
        if contracted != direct:
            raise RuntimeError(f"Independent Poisson calculation disagrees for {name}")
        print(f"{name} = {contracted}")
        print(
            f"  tree terms={len(expression.terms)}, momentum degrees={expression.momentum_degrees}"
        )
    print("All comparisons are complete polynomial equalities, not pointwise tests.")


if __name__ == "__main__":
    main()

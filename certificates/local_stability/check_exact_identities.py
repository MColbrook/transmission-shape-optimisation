#!/usr/bin/env python3
"""Exact arithmetic checks for the analytic local-stability manuscript.

This standard-library-only companion verifies the finite rational identities
printed in the paper. It does not approximate transmission eigenvalues, and
the analytic proof does not depend on running it.
"""
from fractions import Fraction as Q
from math import comb, factorial


def bessel_even_sum(order, radius_squared, last):
    """J_order(sqrt(s))/(sqrt(s)/2)**order, truncated at k=last."""
    return sum(
        (-Q(radius_squared, 4))**k / (factorial(k) * factorial(k + order))
        for k in range(last + 1)
    )


def polynomial_add(*polys):
    result = [Q(0)] * max(map(len, polys))
    for poly in polys:
        for i, coefficient in enumerate(poly):
            result[i] += coefficient
    return tuple(result)


def scale(poly, scalar):
    return tuple(scalar * coefficient for coefficient in poly)


def check(name, actual, expected):
    assert actual == expected, f"{name}: got {actual}; expected {expected}"
    print(f"PASS  {name}: {actual}")


def main():
    check("J0(sqrt(6)) even upper sum", bessel_even_sum(0, 6, 4), Q(-23, 1024))
    check("J1(sqrt(15))/(sqrt(15)/2) even upper sum",
          bessel_even_sum(1, 15, 4), Q(-11, 16384))
    check("J0(sqrt(30)) even upper sum",
          bessel_even_sum(0, 30, 8), Q(-2898599, 205520896))
    check("Rayleigh tail subtraction", Q(1, 192) - Q(1, 225), Q(11, 14400))
    assert Q(11, 14400) < Q(1, 900)
    print("PASS  Rayleigh tail is strictly below 1/900")

    # Integral cos(theta)^(2r) dtheta / pi = 2*binom(2r,r)/4**r.
    angular = lambda r: Q(2 * comb(2*r, r), 4**r)
    check("gradient x angular ratio", angular(2), Q(3, 4))
    check("gradient y angular ratio", angular(1) - angular(2), Q(1, 4))
    check("Hessian xx angular ratio", angular(3), Q(5, 8))
    check("Hessian yy angular ratio", angular(1) - 2*angular(2) + angular(3), Q(1, 8))
    check("Hessian xx-yy angular ratio", angular(2) - angular(3), Q(1, 8))

    t, spectral_ratio = Q(4, 5), Q(15, 16)
    gradient_factor = (3*t + 1/t)/4
    laplacian_factor = (5*t*t + 1/(t*t) + 2)/8
    check("ellipse gradient factor", gradient_factor, Q(73, 80))
    check("ellipse Laplacian factor", laplacian_factor, Q(541, 640))
    # Polynomials are represented by ascending powers of N.
    resonance_trial = polynomial_add(
        scale((1, 0, 1), laplacian_factor/2),
        scale((1, 2, 1), -spectral_ratio*gradient_factor/2),
        (0, spectral_ratio**2, 0),
    )
    check("ellipse trial polynomial in N", resonance_trial,
          (Q(-13, 2560), Q(60, 2560), Q(-13, 2560)))
    check("negative polynomial shifted at N=5",
          polynomial_add(scale((25, -10, 1), -13), (350, -70, 0), (-38, 0, 0)),
          (-13, 60, -13))
    check("contrast derivative sign factor",
          polynomial_add(scale((1, 1), Q(73, 160)), (-Q(15, 16), 0)),
          (Q(-77, 160), Q(73, 160)))
    endpoint_value = polynomial_add(scale((1, 0, 1), Q(1, 2)),
                                    scale((1, 2, 1), -Q(1, 2)), (0, 1, 0))
    check("endpoint F(1)", endpoint_value, (0, 0, 0))
    endpoint_derivative = polynomial_add(scale((1, 0, 1), Q(1, 2)),
                                         scale((1, 2, 1), -Q(1, 4)))
    check("endpoint F'(1)=(N-1)^2/4", endpoint_derivative,
          (Q(1, 4), -Q(1, 2), Q(1, 4)))
    print("All finite rational identities passed. No floating-point arithmetic used.")


if __name__ == "__main__":
    main()

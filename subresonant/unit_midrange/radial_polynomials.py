"""Exact radial and Graeffe coefficients for the subresonant root count."""

from fractions import Fraction as Q
from math import comb, factorial


DEGREE = 24
# The first transform uses this radius. Its coefficients are then rescaled to 14.
SCALED_RADIUS = Q(61, 5)


def polynomial_add(left, right):
    result = [Q(0)] * max(len(left), len(right))
    for index, value in enumerate(left):
        result[index] += value
    for index, value in enumerate(right):
        result[index] += value
    return result


def polynomial_scale(values, factor):
    return [factor * value for value in values]


def polynomial_product(left, right):
    result = [Q(0)] * (len(left) + len(right) - 1)
    for i, value in enumerate(left):
        for j, other in enumerate(right):
            result[i + j] += value * other
    return result


def bernstein_coefficients(coefficients, upper):
    degree = len(coefficients) - 1
    return [
        sum(
            coefficients[j] * upper**j * Q(comb(k, j), comb(degree, j))
            for j in range(k + 1)
        )
        for k in range(degree + 1)
    ]


def canonical_radial_magnitude(index):
    return Q(
        factorial(2 * index),
        4**index * factorial(index)**4 * (index + 1),
    )


def graeffe_parameter_polynomials():
    """Express each transformed coefficient as a polynomial in v."""
    coefficient_polynomials = [
        [
            Q(1, 4**k * (k + 1) * factorial(j)**2 * factorial(k - j)**2)
            for j in range(k + 1)
        ]
        for k in range(DEGREE + 1)
    ]
    chebyshev = [[Q(2)], [Q(2), Q(1)]]
    for _ in range(2, DEGREE + 1):
        chebyshev.append(polynomial_add(
            polynomial_product([Q(2), Q(1)], chebyshev[-1]),
            polynomial_scale(chebyshev[-2], -1),
        ))
    result = []
    for index in range(DEGREE + 1):
        numerator = [Q(0)] * (2 * index + 1)
        first = max(0, 2 * index - DEGREE)
        last = min(DEGREE, 2 * index)
        for k in range(first, last + 1):
            product = polynomial_product(
                coefficient_polynomials[k],
                coefficient_polynomials[2 * index - k],
            )
            numerator = polynomial_add(
                numerator, polynomial_scale(product, (-1)**k)
            )
        assert numerator == numerator[::-1]
        transformed = [numerator[index]]
        for power in range(1, index + 1):
            transformed = polynomial_add(
                transformed,
                polynomial_scale(chebyshev[power], numerator[index + power]),
            )
        result.append(polynomial_scale(
            transformed, SCALED_RADIUS**(2 * index)
        ))
    return result

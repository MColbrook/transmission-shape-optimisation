"""Exact angular integrals and radial Taylor coefficients for the disk Hessian."""

from fractions import Fraction as Q
from math import factorial

from unit_endpoint_check import ComplexInterval as CI
from unit_midrange.radial_polynomials import DEGREE


def complex_polynomial_value(coefficients, value):
    result = CI.point()
    for coefficient in reversed(coefficients):
        result = result * value + coefficient
    return result


def radial_u_coefficients(point, derivative=0):
    """Group the degree-24 radial polynomial by powers of u=t+1/t."""
    powers = [CI.point(1)]
    for _ in range(DEGREE):
        powers.append(powers[-1] * point)
    result = [CI.point() for _ in range(DEGREE + 1)]
    for order in range(derivative, DEGREE + 1):
        for index in range(order // 2 + 1):
            power = order - 2 * index
            coefficient = (-1)**order * Q(
                factorial(order) // factorial(order - derivative),
                4**order * (order + 1) * factorial(order)
                * factorial(index)**2 * factorial(power),
            )
            result[power] = result[power] + (
                powers[order - derivative] * CI.point(coefficient)
            )
    return result


def angular_integral_coefficient(mode, order, t):
    return Q(factorial(mode)**2, 4**order * (mode + order + 1)) * sum(
        t**(2 * index - order) / (
            factorial(index) * factorial(order - index)
            * factorial(mode + index) * factorial(mode + order - index)
        )
        for index in range(order + 1)
    )


def coordinate_absolute_lower(interval):
    if interval.lower <= 0 <= interval.upper:
        return Q(0)
    return min(abs(interval.lower), abs(interval.upper))

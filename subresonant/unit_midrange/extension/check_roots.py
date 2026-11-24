"""Uniform exact pair isolation on an extended contrast interval."""

from fractions import Fraction
import json
from math import factorial
from pathlib import Path
import sys


sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from unit_endpoint_check import ComplexInterval, Interval, normalised_bessel
from unit_midrange import radial_polynomials


Q = Fraction
CI = ComplexInterval
I = Interval
DEGREE = radial_polynomials.DEGREE
RADIUS = Q(14)
T_UPPER = Q(21, 10)
V_UPPER = (T_UPPER - 1 / T_UPPER)**2
TAIL_ARGUMENT = Q(30)


def entire_tail(derivative=0):
    first = DEGREE + 1
    value = radial_polynomials.canonical_radial_magnitude(first)
    value *= TAIL_ARGUMENT**(first - derivative)
    value *= factorial(first) // factorial(first - derivative)
    ratio = TAIL_ARGUMENT * Q(
        2 * first + 1, 2 * (first + 1)**2 * (first + 2)
    ) * Q(first + 1, first + 1 - derivative)
    assert 0 < ratio < 1
    return value / (1 - ratio)


def transformed_coefficients():
    first = [
        radial_polynomials.polynomial_scale(
            polynomial, (RADIUS / radial_polynomials.SCALED_RADIUS)**(2 * order)
        )
        for order, polynomial in enumerate(radial_polynomials.graeffe_parameter_polynomials())
    ]
    second = []
    for order in range(DEGREE + 1):
        value = [Q(0)]
        lower = max(0, 2 * order - DEGREE)
        upper = min(DEGREE, 2 * order)
        for k in range(lower, upper + 1):
            value = radial_polynomials.polynomial_add(
                value,
                radial_polynomials.polynomial_scale(
                    radial_polynomials.polynomial_product(first[k], first[2 * order - k]),
                    (-1)**k,
                ),
            )
        second.append(value)
    return second


def rayleigh_sums(order, count):
    series = [
        (-1)**index * Q(
            factorial(order),
            4**index * factorial(index) * factorial(order + index),
        )
        for index in range(count + 1)
    ]
    result = [Q(0)]
    for index in range(1, count + 1):
        result.append(
            -index * series[index]
            -sum(series[j] * result[index - j] for j in range(1, index))
        )
    return result


def angular_one_endpoint_sign(lower, upper):
    x = CI(I(RADIUS / upper, RADIUS / lower), I.point(0))
    y = CI(I(RADIUS * lower, RADIUS * upper), I.point(0))
    first_x = CI(normalised_bessel(1, x, TAIL_ARGUMENT).real, I.point(0))
    second_x = CI(normalised_bessel(2, x, TAIL_ARGUMENT).real, I.point(0))
    first_y = CI(normalised_bessel(1, y, TAIL_ARGUMENT).real, I.point(0))
    second_y = CI(normalised_bessel(2, y, TAIL_ARGUMENT).real, I.point(0))
    determinant = (x * second_x * first_y - y * second_y * first_x) * CI.point(Q(1, 4))
    assert determinant.imag.lower == determinant.imag.upper == 0
    assert determinant.real.upper < 0, (
        "first-angular endpoint", str(lower), str(upper),
        float(determinant.real.upper),
    )
    return determinant.real.upper


def main():
    assert RADIUS * T_UPPER < TAIL_ARGUMENT
    transformed = transformed_coefficients()
    partial_margin = transformed[2]
    mixed_bound = Q(0)
    mixed_indices = []
    for index, polynomial in enumerate(transformed):
        if index == 2:
            continue
        coefficients = radial_polynomials.bernstein_coefficients(polynomial, V_UPPER)
        if min(coefficients) >= 0:
            sign = 1
        elif max(coefficients) <= 0:
            sign = -1
        else:
            mixed_bound += max(abs(value) for value in coefficients)
            mixed_indices.append(index)
            continue
        partial_margin = radial_polynomials.polynomial_add(
            partial_margin, radial_polynomials.polynomial_scale(polynomial, -sign)
        )
    coefficients = radial_polynomials.bernstein_coefficients(partial_margin, V_UPPER)
    margin = min(coefficients) - mixed_bound
    assert margin > Q(89, 100)
    norm_upper = sum(
        radial_polynomials.canonical_radial_magnitude(index) * TAIL_ARGUMENT**index
        for index in range(DEGREE + 1)
    )
    boundary_lower = margin / norm_upper**3
    assert entire_tail() < boundary_lower

    first_rayleigh = rayleigh_sums(1, 3)
    assert first_rayleigh[2] == Q(1, 192)
    assert first_rayleigh[3] == Q(1, 3072)
    assert 14**3 < 3072
    second_rayleigh = rayleigh_sums(2, 4)
    assert second_rayleigh[4] == Q(7, 3317760)
    assert 26**4 * second_rayleigh[4] < 1
    third_rayleigh = rayleigh_sums(3, 2)
    assert third_rayleigh[2] == Q(1, 1280)
    assert (RADIUS * T_UPPER)**2 < 1280

    # The first J_1 zero is below 15 and its second squared zero exceeds 36.
    first_even = sum(
        (-Q(15, 4))**j / (factorial(j) * factorial(j + 1))
        for j in range(5)
    )
    assert first_even == -Q(11, 16384)
    second_zero_inverse_square_upper = Q(1, 192) - Q(1, 225)
    assert second_zero_inverse_square_upper == Q(11, 14400)
    assert 36**2 * second_zero_inverse_square_upper < 1
    assert RADIUS * T_UPPER < 36

    # Identify the second J_0 zero above squared argument 30.
    radial_at_six = sum(
        (-Q(3, 2))**j / factorial(j)**2 for j in range(5)
    )
    radial_at_thirty = sum(
        (-Q(15, 2))**j / factorial(j)**2 for j in range(9)
    )
    assert radial_at_six < 0
    assert radial_at_thirty < 0
    assert 6 * T_UPPER**2 < 30

    lower = Q(13, 7)
    assert RADIUS * lower == 26
    endpoint_bounds = []
    while lower < T_UPPER:
        upper = min(T_UPPER, lower + Q(1, 1000))
        endpoint_bounds.append(angular_one_endpoint_sign(lower, upper))
        lower = upper
    assert lower == T_UPPER

    result = {
        "status": "passed",
        "arithmetic": "fractions.Fraction only",
        "contrast_interval": ["1", str(T_UPPER**2)],
        "sqrt_contrast_interval": ["1", str(T_UPPER)],
        "v_interval": ["0", str(V_UPPER)],
        "scaled_count_radius": str(RADIUS),
        "degree": DEGREE,
        "graeffe_steps": 2,
        "graeffe_margin_lower": str(margin),
        "mixed_sign_coefficients": mixed_indices,
        "polynomial_boundary_lower": str(boundary_lower),
        "entire_tail_upper": str(entire_tail()),
        "radial_count": 2,
        "first_angular_sign_intervals": len(endpoint_bounds),
        "first_angular_maximum_determinant_upper": str(max(endpoint_bounds)),
        "all_other_angular_orders": "Analytic Dirichlet and Bessel bounds",
        "radial_real_collision": "Excluded by alpha<6 and zeta>30",
        "active_pair": "Distinct, nonreal and algebraically simple",
        "shape_optimality": "Requires a separate Hessian check",
    }
    Path(__file__).with_name("root_certificate.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "status": "passed",
        "contrast_interval": result["contrast_interval"],
        "graeffe_margin": float(margin),
        "boundary_lower": float(boundary_lower),
        "entire_tail_upper": float(entire_tail()),
        "first_angular_sign_intervals": len(endpoint_bounds),
        "first_angular_determinant_upper": float(max(endpoint_bounds)),
    }, indent=2))


if __name__ == "__main__":
    main()

"""Exact all-mode modulus Hessian bounds on the extended interval."""

from fractions import Fraction
import json
from math import isqrt
from pathlib import Path
import sys


sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from unit_midrange.extension import check_roots as roots
from unit_midrange.hessian_series import (
    angular_integral_coefficient,
    complex_polynomial_value,
    coordinate_absolute_lower,
    radial_u_coefficients,
)
from unit_endpoint_check import ComplexInterval, Interval, normalised_bessel

Q = Fraction
CI = ComplexInterval
I = Interval
BESSEL_ARGUMENT = Q(25)


def norm_upper(real, imag):
    squared = real*real + imag*imag
    scale = 10**8
    integer = isqrt(squared.numerator * scale**2 // squared.denominator) + 1
    result = Q(integer, scale)
    assert squared < result**2
    return result


def angular_integral(mode, point, lower, upper):
    result = CI.point()
    for order in reversed(range(roots.DEGREE + 1)):
        first = angular_integral_coefficient(mode, order, lower)
        last = angular_integral_coefficient(mode, order, upper)
        assert first <= last
        coefficient = I(first, last) if order % 2 == 0 else I(-last, -first)
        result = result * point + CI(coefficient, I.point(0))
    return result.enlarged(roots.entire_tail())


def verify_tube(candidate):
    lower = Q(candidate["t_lower"])
    upper = Q(candidate["t_upper"])
    real = Q(candidate["mu_real"])
    imag = Q(candidate["mu_imag"])
    radius = Q(candidate["radius"])
    centre_norm = norm_upper(real, imag)
    local_radius = centre_norm + 2 * radius
    assert imag > radius
    assert local_radius < roots.RADIUS
    assert local_radius * upper < BESSEL_ARGUMENT
    point = CI.point(real, imag)
    u_lower = lower + 1 / lower
    u_upper = upper + 1 / upper
    u_midpoint = (u_lower + u_upper) / 2
    u_interval = CI(I(u_lower, u_upper), I.point(0))
    grouped = radial_u_coefficients(point)
    u_derivative = [
        grouped[index] * CI.point(index)
        for index in range(1, len(grouped))
    ]
    residual = complex_polynomial_value(grouped, CI.point(u_midpoint)).absolute_upper()
    residual += (u_upper - u_lower) * complex_polynomial_value(
        u_derivative, u_interval
    ).absolute_upper() / 2
    residual += roots.entire_tail()
    derivative = complex_polynomial_value(
        radial_u_coefficients(point, 1), u_interval
    )
    derivative_lower = max(
        coordinate_absolute_lower(derivative.real),
        coordinate_absolute_lower(derivative.imag),
    ) - roots.T_UPPER * roots.entire_tail(1)
    assert derivative_lower > 0
    second_upper = sum(
        order * (order - 1) * angular_integral_coefficient(0, order, upper)
        * local_radius**(order - 2)
        for order in range(2, roots.DEGREE + 1)
    ) + roots.T_UPPER**2 * roots.entire_tail(2)
    error = residual + radius*radius * second_upper / 2
    linear = radius * derivative_lower
    assert error < linear, (
        "root tube", candidate, float(error / linear)
    )

    rectangle = CI(
        I(real - radius, real + radius),
        I(imag - radius, imag + radius),
    )
    x = rectangle * CI(I(1 / upper, 1 / lower), I.point(0))
    y = rectangle * CI(I(lower, upper), I.point(0))
    first = {
        mode: normalised_bessel(mode, x, BESSEL_ARGUMENT)
        for mode in range(7)
    }
    second = {
        mode: normalised_bessel(mode, y, BESSEL_ARGUMENT)
        for mode in range(2, 7)
    }
    kappa = -x * first[1] / (CI.point(2) * first[0])
    assert kappa.real.upper < 3, ("kappa real", candidate)
    assert kappa.absolute_upper() < 10, ("kappa modulus", candidate)
    coefficient_lower = []
    coefficient_upper = []
    for mode in range(2, 7):
        integral = angular_integral(mode, rectangle, lower, upper)
        value = -kappa + CI.point(2) * first[mode] * second[mode] / integral
        scaled_lower = value.real.lower / (mode + 1)
        scaled_upper = value.absolute_upper() / (mode + 1)
        assert scaled_lower > Q(1, 4), (
            "Hessian lower", mode, candidate, float(scaled_lower)
        )
        assert scaled_upper < 14, (
            "Hessian upper", mode, candidate, float(scaled_upper)
        )
        coefficient_lower.append(scaled_lower)
        coefficient_upper.append(scaled_upper)
    return {
        "root_margin_ratio": (linear - error) / linear,
        "coefficient_lower": min(coefficient_lower),
        "coefficient_upper": max(coefficient_upper),
    }


def main():
    data = json.loads(Path(__file__).with_name("centres.json").read_text(encoding="utf-8"))
    candidates = data["centres"]
    boundary = Q(1)
    results = []
    for index, candidate in enumerate(candidates):
        assert Q(candidate["t_lower"]) == boundary
        assert boundary < Q(candidate["t_upper"]) <= roots.T_UPPER
        results.append(verify_tube(candidate))
        boundary = Q(candidate["t_upper"])
        if (index + 1) % 25 == 0:
            print("Verified extended tubes", index + 1, flush=True)
    assert boundary == roots.T_UPPER
    # For m>=7, A_{m,1}>4(m+1)*sqrt(m+2)>=12(m+1).
    relative_error = (1 - BESSEL_ARGUMENT / (12 * 8))**(-2) - 1
    assert relative_error == Q(4175, 5041)
    assert relative_error < 1
    assert 2 / (1 + relative_error) - Q(3, 8) > Q(1, 4)
    assert 2 / (1 - relative_error) + Q(10, 8) < 14
    result = {
        "status": "passed",
        "arithmetic": "fractions.Fraction only",
        "contrast_interval": ["1", str(roots.T_UPPER**2)],
        "parameter_tubes": len(candidates),
        "coverage": "Closed rational intervals with exact shared endpoints",
        "minimum_root_margin_ratio": str(min(item["root_margin_ratio"] for item in results)),
        "minimum_scaled_finite_coefficient_lower": str(min(item["coefficient_lower"] for item in results)),
        "maximum_scaled_finite_coefficient_upper": str(max(item["coefficient_upper"] for item in results)),
        "finite_modes": [2, 3, 4, 5, 6],
        "tail_modes": "every m>=7",
        "bessel_argument_upper": str(BESSEL_ARGUMENT),
        "kappa_bounds": ["Re(kappa)<3", "|kappa|<10"],
        "tail_relative_error": str(relative_error),
        "all_mode_bounds": ["Re(B_m)>(m+1)/4", "|B_m|<14(m+1)"],
    }
    Path(__file__).with_name("hessian_certificate.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "status": "passed",
        "contrast_interval": result["contrast_interval"],
        "parameter_tubes": len(candidates),
        "minimum_root_margin_ratio": float(min(item["root_margin_ratio"] for item in results)),
        "minimum_scaled_finite_coefficient_lower": float(min(item["coefficient_lower"] for item in results)),
        "maximum_scaled_finite_coefficient_upper": float(max(item["coefficient_upper"] for item in results)),
    }, indent=2))


if __name__ == "__main__":
    main()

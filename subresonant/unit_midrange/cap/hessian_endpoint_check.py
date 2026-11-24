"""Exact near-resonance root chart and two logarithmic-derivative cones."""

from fractions import Fraction as Q
import json
from math import factorial
from pathlib import Path
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import unit_endpoint_check as unit

I = unit.Interval
C = unit.ComplexInterval


def real(value):
    return C(value, I.point(0))


def scale(value, factor):
    return value * C.point(factor)


def poly_add(first, second, sign=1):
    length = max(len(first), len(second))
    result = [C.point() for _ in range(length)]
    for index, value in enumerate(first):
        result[index] = result[index] + value
    for index, value in enumerate(second):
        result[index] = result[index] + scale(value, sign)
    return result


def poly_multiply(first, second):
    result = [C.point() for _ in range(len(first) + len(second) - 1)]
    for index, value in enumerate(first):
        for other, coefficient in enumerate(second):
            result[index + other] = result[index + other] + value * coefficient
    return result


def poly_compose(coefficients, argument):
    result = [C.point()]
    for coefficient in reversed(coefficients):
        result = poly_multiply(result, argument)
        result[0] = result[0] + real(coefficient)
    return result


def series_coefficients(argument, degree):
    coefficients = [I.point(0), I.point(1)]
    for index in range(degree - 1):
        numerator = coefficients[index + 1] * I.point((index + 1)**2)
        numerator = numerator + coefficients[index] * I.point(Q(1, 4))
        denominator = argument * I.point((index + 2) * (index + 1))
        coefficients.append((-numerator) * denominator.reciprocal())
    return coefficients


def tail_constants(argument, coefficients, quotient_bound, eta_maximum):
    """Bound tails by eta^(M+1)*F and eta^M*Fprime."""
    degree = len(coefficients) - 1
    minimum = argument.lower
    assert degree - 1 > argument.upper / 4
    majorant = max(
        coefficients[degree - 1].absolute_upper() * minimum**(degree - 2),
        coefficients[degree].absolute_upper() * minimum**(degree - 1),
    )
    qmax = eta_maximum * quotient_bound / minimum
    assert qmax < 1
    field = majorant * minimum * (quotient_bound / minimum)**(degree + 1)
    field /= 1 - qmax
    derivative = majorant * (quotient_bound / minimum)**degree
    derivative *= Q(degree + 1, 1) / (1 - qmax) + qmax / (1 - qmax)**2
    return field, derivative


def polynomial_absolute_bound(coefficients, argument):
    return sum(
        coefficient.absolute_upper() * argument**index
        for index, coefficient in enumerate(coefficients)
    )


def real_bessel(argument, degree=40):
    total = sum(
        (-argument / 4)**index / factorial(index)**2
        for index in range(degree + 1)
    )
    first = degree + 1
    tail = (argument / 4)**first / factorial(first)**2
    tail /= 1 - argument / (4 * (first + 1)**2)
    return I(total - tail, total + tail)


def norm_squared_upper(value):
    return value.real.absolute_upper()**2 + value.imag.absolute_upper()**2


def main():
    alpha = I(Q("5.78318596294"), Q("5.78318596296"))
    zeta = I(Q("30.47126234365"), Q("30.47126234368"))
    assert real_bessel(alpha.lower).lower > 0
    assert real_bessel(alpha.upper).upper < 0
    assert real_bessel(zeta.lower).upper < 0
    assert real_bessel(zeta.upper).lower > 0
    contrast = zeta * alpha.reciprocal()
    eta_maximum = Q(1, 40)
    eta = I(Q(0), eta_maximum)
    centre = C.point(Q(-130655, 100000), Q(226302, 100000))
    radius = Q(1, 10)
    chart = centre.enlarged(radius)
    z_bound = chart.absolute_upper()
    degree = 12
    ca = series_coefficients(alpha, degree)
    cb = series_coefficients(zeta, degree)
    hx = [C.point(), chart]
    hy = [
        C.point(), real(contrast) * chart, C.point(), -real(alpha), -chart
    ]
    x = [real(alpha), chart]
    y = hy.copy()
    y[0] = real(zeta)
    fa = poly_compose(ca, hx)
    fb = poly_compose(cb, hy)
    dca = [coefficient * I.point(index) for index, coefficient in enumerate(ca) if index]
    dcb = [coefficient * I.point(index) for index, coefficient in enumerate(cb) if index]
    dfa = poly_compose(dca, hx)
    dfb = poly_compose(dcb, hy)
    determinant = poly_add(
        poly_multiply(poly_multiply(x, dfa), fb),
        poly_multiply(poly_multiply(y, dfb), fa),
        sign=-1,
    )
    # The first three coefficients vanish identically when zeta=N*alpha.
    for coefficient in determinant[:3]:
        assert coefficient.real.lower <= 0 <= coefficient.real.upper
        assert coefficient.imag.lower <= 0 <= coefficient.imag.upper
    leading = contrast * (contrast - I.point(1)) * I.point(Q(1, 12))
    leading_midpoint = (leading.lower + leading.upper) / 2
    alpha_midpoint = (alpha.lower + alpha.upper) / 2
    leading_width = max(
        abs(leading.lower - leading_midpoint),
        abs(leading.upper - leading_midpoint),
    )
    alpha_square_width = max(
        abs(alpha.square().lower - alpha_midpoint**2),
        abs(alpha.square().upper - alpha_midpoint**2),
    )
    parameter_error = leading_width * z_bound**3 + alpha_square_width
    polynomial_error = sum(
        coefficient.absolute_upper() * eta_maximum**(index - 3)
        for index, coefficient in enumerate(determinant) if index >= 4
    )

    hy_quotient_bound = contrast.upper * z_bound + alpha.upper * eta_maximum**2
    ta, dta = tail_constants(alpha, ca, z_bound, eta_maximum)
    tb, dtb = tail_constants(zeta, cb, hy_quotient_bound, eta_maximum)
    fa_quotient_bound = sum(
        ca[index].absolute_upper() * eta_maximum**(index - 1) * z_bound**index
        for index in range(1, degree + 1)
    )
    fb_quotient_bound = sum(
        cb[index].absolute_upper() * eta_maximum**(index - 1)
        * hy_quotient_bound**index for index in range(1, degree + 1)
    )
    dfa_bound = sum(
        index * ca[index].absolute_upper()
        * (eta_maximum * z_bound)**(index - 1)
        for index in range(1, degree + 1)
    )
    dfb_bound = sum(
        index * cb[index].absolute_upper()
        * (eta_maximum * hy_quotient_bound)**(index - 1)
        for index in range(1, degree + 1)
    )
    x_bound = alpha.upper + eta_maximum * z_bound
    y_bound = zeta.upper + eta_maximum * hy_quotient_bound
    assert contrast.upper * x_bound**2 < 14**2
    tail_error = (
        x_bound * (dta * (fb_quotient_bound + eta_maximum**degree * tb) + dfa_bound * tb)
        + y_bound * (dtb * (fa_quotient_bound + eta_maximum**degree * ta) + dfb_bound * ta)
    ) * eta_maximum**(degree - 2)

    centre_norm_squared = centre.real.lower**2 + centre.imag.lower**2
    centre_norm_upper = centre.absolute_upper()
    cubic_residual = (
        scale(centre * centre * centre, leading_midpoint)
        - C.point(alpha_midpoint**2)
    ).absolute_upper()
    contour_lower = (
        3 * leading_midpoint * centre_norm_squared * radius
        - 3 * leading_midpoint * centre_norm_upper * radius**2
        - leading_midpoint * radius**3 - cubic_residual
    )
    total_error = parameter_error + polynomial_error + tail_error
    assert contour_lower > total_error

    # f_alpha(eta*z)/eta and f'_alpha(eta*z), including eta=0.
    field_over_eta = C.point()
    derivative = C.point()
    eta_power = I.point(1)
    chart_power = C.point(1)
    for index in range(1, degree + 1):
        derivative_power = chart_power
        chart_power = chart_power * chart
        field_over_eta = field_over_eta + real(ca[index] * eta_power) * chart_power
        derivative = derivative + real(ca[index] * I.point(index) * eta_power) * (
            derivative_power
        )
        eta_power = eta_power * eta
    field_over_eta = field_over_eta.enlarged(eta_maximum**degree * ta)
    derivative = derivative.enlarged(eta_maximum**degree * dta)
    x_box = real(alpha) + real(eta) * chart
    kappa_normalised = scale(x_box * derivative / field_over_eta, 2)
    quotient_normalised = field_over_eta / scale(derivative, 2)
    assert kappa_normalised.real.upper < 0
    assert quotient_normalised.real.upper < 0
    first_cone_margin = 36 * kappa_normalised.real.upper**2 - norm_squared_upper(kappa_normalised)
    second_cone_margin = 16 * quotient_normalised.real.upper**2 - norm_squared_upper(quotient_normalised)
    square_real_upper = (
        quotient_normalised.real.square() - quotient_normalised.imag.square()
    ).upper
    assert first_cone_margin > 0
    assert second_cone_margin > 0
    assert square_real_upper < 0
    compact_end = Q(459083, 200000)**2
    endpoint_start_upper = contrast.upper - eta_maximum**3
    assert endpoint_start_upper < compact_end < contrast.lower
    assert contrast.lower - eta_maximum**3 > Q(441, 100)
    assert (Q(21, 10) + Q(10, 21))**2 * 104 > 26**2
    assert Q(13, 12) > Q(26, 25)**2
    result = {
        "status": "passed exact rational arithmetic",
        "alpha_interval": [str(alpha.lower), str(alpha.upper)],
        "zeta_interval": [str(zeta.lower), str(zeta.upper)],
        "eta_interval": ["0", str(eta_maximum)],
        "root_chart_centre": [str(centre.real.lower), str(centre.imag.lower)],
        "root_chart_radius": str(radius),
        "ODE_series_degree": degree,
        "normalised_determinant": "D/eta^3 = N(N-1)z^3/12-alpha^2+O(eta)",
        "rouche_relative_margin_lower": unit.decimal_interval(
            I.point((contour_lower - total_error) / contour_lower), 10
        )[0],
        "tail_error_upper": str(tail_error),
        "first_cone_squared_margin_lower": unit.decimal_interval(I.point(first_cone_margin), 10)[0],
        "second_cone_squared_margin_lower": unit.decimal_interval(I.point(second_cone_margin), 10)[0],
        "normalised_quotient_square_real_upper": unit.decimal_interval(I.point(square_real_upper), 10)[1],
        "normalised_kappa_box": {
            "real": unit.decimal_interval(kappa_normalised.real, 10),
            "imag": unit.decimal_interval(kappa_normalised.imag, 10),
        },
        "normalised_quotient_box": {
            "real": unit.decimal_interval(quotient_normalised.real, 10),
            "imag": unit.decimal_interval(quotient_normalised.imag, 10),
        },
        "compact_chart_endpoint": str(compact_end),
        "endpoint_chart_start_upper": str(endpoint_start_upper),
        "coverage": "Every 0<eta<=1/40, including the compact-chart overlap",
    }
    destination = Path(__file__).with_name("hessian_endpoint_certificate.json")
    destination.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in result.items() if key != "tail_error_upper"}, indent=2))


if __name__ == "__main__":
    main()

"""Exact checks for the unit-contrast disk roots and modulus Hessian."""

from dataclasses import dataclass
from fractions import Fraction
import json
from math import factorial
from pathlib import Path


Q = Fraction
DEGREE = 30
ROOT_RADIUS = Q(1, 10**10)
ROOT_REAL = Q(7245296677252, 10**12)
ROOT_IMAG = Q(7627409818058, 10**12)
COUNT_RADIUS = Q(11)


@dataclass(frozen=True)
class Interval:
    lower: Fraction
    upper: Fraction

    @staticmethod
    def point(value):
        value = Q(value)
        return Interval(value, value)

    def __add__(self, other):
        return Interval(self.lower + other.lower, self.upper + other.upper)

    def __neg__(self):
        return Interval(-self.upper, -self.lower)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        products = (
            self.lower * other.lower,
            self.lower * other.upper,
            self.upper * other.lower,
            self.upper * other.upper,
        )
        return Interval(min(products), max(products))

    def reciprocal(self):
        assert self.lower > 0 or self.upper < 0
        return Interval(1 / self.upper, 1 / self.lower)

    def square(self):
        values = (self.lower**2, self.upper**2)
        lower = 0 if self.lower <= 0 <= self.upper else min(values)
        return Interval(Q(lower), max(values))

    def absolute_upper(self):
        return max(abs(self.lower), abs(self.upper))

    def enlarged(self, error):
        return Interval(self.lower - error, self.upper + error)


@dataclass(frozen=True)
class ComplexInterval:
    real: Interval
    imag: Interval

    @staticmethod
    def point(real=0, imag=0):
        return ComplexInterval(Interval.point(real), Interval.point(imag))

    def __add__(self, other):
        return ComplexInterval(self.real + other.real, self.imag + other.imag)

    def __neg__(self):
        return ComplexInterval(-self.real, -self.imag)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        return ComplexInterval(
            self.real * other.real - self.imag * other.imag,
            self.real * other.imag + self.imag * other.real,
        )

    def reciprocal(self):
        denominator = (self.real.square() + self.imag.square()).reciprocal()
        return ComplexInterval(
            self.real * denominator, (-self.imag) * denominator
        )

    def __truediv__(self, other):
        return self * other.reciprocal()

    def enlarged(self, error):
        return ComplexInterval(
            self.real.enlarged(error), self.imag.enlarged(error)
        )

    def absolute_upper(self):
        return self.real.absolute_upper() + self.imag.absolute_upper()


def convolution(left, right):
    result = [Q(0)] * (len(left) + len(right) - 1)
    for i, value in enumerate(left):
        for j, other in enumerate(right):
            result[i + j] += value * other
    return result


def graeffe(coefficients):
    """Return E(z)^2-z O(z)^2 for P(z)=E(z^2)+z O(z^2)."""
    even = convolution(coefficients[::2], coefficients[::2])
    odd = convolution(coefficients[1::2], coefficients[1::2])
    result = [Q(0)] * len(coefficients)
    for j, value in enumerate(even):
        result[j] += value
    for j, value in enumerate(odd):
        result[j + 1] -= value
    return result


def radial_coefficient(j):
    return (-1)**j * Q(
        factorial(2 * j), 4**j * factorial(j)**4 * (j + 1)
    )


def radial_tail(radius, derivative=0):
    """Bound the derivative tail after degree DEGREE on a closed disk."""
    first = DEGREE + 1
    value = abs(radial_coefficient(first)) * radius**(first - derivative)
    value *= factorial(first) // factorial(first - derivative)
    ratio = radius * Q(2 * first + 1, 2 * (first + 1)**2 * (first + 2))
    ratio *= Q(first + 1, first + 1 - derivative)
    assert 0 < ratio < 1
    return value / (1 - ratio)


def polynomial_at(coefficients, value):
    result = ComplexInterval.point()
    for coefficient in reversed(coefficients):
        result = result * value + ComplexInterval.point(coefficient)
    return result


def normalised_bessel(order, value, radius):
    """Enclose m! (2/sqrt(lambda))^m J_m(sqrt(lambda))."""
    term = ComplexInterval.point(1)
    total = term
    for j in range(1, DEGREE + 1):
        term = term * (-value) / ComplexInterval.point(4 * j * (order + j))
        total = total + term
    first = DEGREE + 1
    next_size = radius**first * Q(
        factorial(order), 4**first * factorial(first) * factorial(order + first)
    )
    ratio = radius / (4 * (first + 1) * (order + first + 1))
    assert 0 < ratio < 1
    return total.enlarged(next_size / (1 - ratio))


def decimal_interval(interval, places=8):
    scale = 10**places
    lower = interval.lower.numerator * scale // interval.lower.denominator
    upper = -((-interval.upper.numerator * scale) // interval.upper.denominator)
    return [str(Q(lower, scale)), str(Q(upper, scale))]


def main():
    coefficients = [radial_coefficient(j) for j in range(DEGREE + 1)]
    scaled = [value * COUNT_RADIUS**j for j, value in enumerate(coefficients)]
    transformed = graeffe(scaled)
    pellet_margin = abs(transformed[2]) - sum(
        abs(value) for j, value in enumerate(transformed) if j != 2
    )
    assert pellet_margin > 0
    polynomial_upper = sum(abs(value) for value in scaled)
    polynomial_lower = pellet_margin / polynomial_upper
    entire_tail = radial_tail(COUNT_RADIUS)
    assert entire_tail < polynomial_lower

    centre = ComplexInterval.point(ROOT_REAL, ROOT_IMAG)
    residual = polynomial_at(coefficients, centre).absolute_upper() + entire_tail
    derivative_coefficients = [
        j * coefficients[j] for j in range(1, len(coefficients))
    ]
    derivative = polynomial_at(derivative_coefficients, centre)
    assert derivative.real.lower == derivative.real.upper
    assert derivative.imag.lower == derivative.imag.upper
    derivative_lower = max(
        abs(derivative.real.lower), abs(derivative.imag.lower)
    ) - radial_tail(COUNT_RADIUS, 1)
    assert derivative_lower > 0
    second_upper = sum(
        j * (j - 1) * abs(coefficients[j]) * COUNT_RADIUS**(j - 2)
        for j in range(2, len(coefficients))
    ) + radial_tail(COUNT_RADIUS, 2)
    local_error = residual + ROOT_RADIUS**2 * second_upper / 2
    linear_lower = ROOT_RADIUS * derivative_lower
    assert local_error < linear_lower
    assert ROOT_IMAG > ROOT_RADIUS
    assert ROOT_REAL**2 + ROOT_IMAG**2 < (Q(1053, 100) - ROOT_RADIUS)**2

    root_box = ComplexInterval(
        Interval(ROOT_REAL - ROOT_RADIUS, ROOT_REAL + ROOT_RADIUS),
        Interval(ROOT_IMAG - ROOT_RADIUS, ROOT_IMAG + ROOT_RADIUS),
    )
    bessel = {
        order: normalised_bessel(order, root_box, COUNT_RADIUS)
        for order in range(8)
    }
    kappa = -root_box * bessel[1] / (
        ComplexInterval.point(2) * bessel[0]
    )
    assert kappa.real.upper < Q(13, 10)
    finite_modes = {}
    for order in range(2, 7):
        p = ComplexInterval.point(order) - root_box * bessel[order + 1] / (
            ComplexInterval.point(2 * (order + 1)) * bessel[order]
        )
        derivative_denominator = ComplexInterval.point(order**2) - root_box - p*p
        coefficient = -kappa - ComplexInterval.point(2) * root_box / derivative_denominator
        assert coefficient.real.lower > Q(9, 10) * (order + 1)
        assert coefficient.absolute_upper() < 6 * (order + 1)
        finite_modes[str(order)] = {
            "real_enclosure": decimal_interval(coefficient.real),
            "imag_enclosure": decimal_interval(coefficient.imag),
        }

    relative_error = (1 - Q(53, 5) / 7**2)**(-2) - 1
    assert 0 < relative_error < 1
    tail_lower_slope = 2 / (1 + relative_error)
    assert (tail_lower_slope - 1) * 8 > Q(13, 10)
    tail_upper_slope = 2 / (1 - relative_error)
    assert tail_upper_slope + Q(4, 8) < 6

    result = {
        "arithmetic": "fractions.Fraction only",
        "degree": DEGREE,
        "radial_count_disk_radius": str(COUNT_RADIUS),
        "radial_count": 2,
        "graeffe_pellet_margin_lower": decimal_interval(
            Interval.point(pellet_margin)
        )[0],
        "polynomial_boundary_lower": decimal_interval(
            Interval.point(polynomial_lower), 12
        )[0],
        "entire_tail_upper": str(entire_tail),
        "upper_root_centre": [str(ROOT_REAL), str(ROOT_IMAG)],
        "root_enclosure_disk_radius": str(ROOT_RADIUS),
        "root_modulus_upper": "1053/100",
        "local_rouche_error": str(local_error),
        "local_rouche_linear_lower": str(linear_lower),
        "kappa_real_enclosure": decimal_interval(kappa.real),
        "finite_hessian_modes": finite_modes,
        "all_hessian_modes_m_ge_2": {
            "lower": "(9/10)(m+1)",
            "absolute_upper": "6(m+1)",
            "analytic_tail_from_m": 7,
        },
        "status": "Exact checks passed. Analytic implications are in unit_analysis.tex.",
    }
    destination = Path(__file__).with_name("unit_endpoint_certificate.json")
    destination.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "radial_root_count_inside_11": 2,
        "least_modulus_pair_upper": "10.53",
        "finite_hessian_modes": finite_modes,
        "all_mode_lower": "(9/10)(m+1)",
        "all_mode_upper": "6(m+1)",
        "status": "passed",
    }, indent=2))


if __name__ == "__main__":
    main()

"""Exact root and cone verification on a compact subresonant interval."""

from fractions import Fraction
import json
from math import factorial, isqrt
from pathlib import Path
import sys


sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
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
DEGREE = 24
TAIL_ARGUMENT = Q(31)
T_UPPER = Q(23, 10)
COUNT_RADIUS = Q(53, 4)
FINAL_T = Q(459083, 200000)


def norm_upper(real, imag):
    squared = real*real+imag*imag
    scale = 10**10
    integer = isqrt(squared.numerator*scale**2//squared.denominator)+1
    result = Q(integer, scale)
    assert squared < result*result
    return result


def rounded_bound(value, upper=False):
    denominator = 10**12
    integer = value.numerator*denominator//value.denominator
    if upper:
        integer += 1
    return Q(integer, denominator)


def entire_tail(derivative=0):
    first = DEGREE+1
    coefficient = Q(factorial(2*first),
                    4**first*factorial(first)**4*(first+1))
    value = coefficient*TAIL_ARGUMENT**(first-derivative)
    value *= factorial(first)//factorial(first-derivative)
    ratio = TAIL_ARGUMENT*Q(2*first+1, 2*(first+1)**2*(first+2))
    ratio *= Q(first+1, first+1-derivative)
    assert 0 < ratio < 1
    return value/(1-ratio)


def verify_tube(candidate):
    lower = Q(candidate["t_lower"])
    upper = Q(candidate["t_upper"])
    real = Q(candidate["mu_real"])
    imag = Q(candidate["mu_imag"])
    radius = Q(candidate["radius"])
    local_radius = norm_upper(real, imag)+2*radius
    assert imag > radius > 0
    assert local_radius < COUNT_RADIUS
    assert local_radius*upper < TAIL_ARGUMENT
    point = CI.point(real, imag)
    u_lower = lower+1/lower
    u_upper = upper+1/upper
    u_midpoint = (u_lower+u_upper)/2
    u_interval = CI(I(u_lower, u_upper), I.point(0))
    grouped = radial_u_coefficients(point)
    u_derivative = [grouped[j]*CI.point(j) for j in range(1, len(grouped))]
    residual = complex_polynomial_value(grouped, CI.point(u_midpoint)).absolute_upper()
    residual += (u_upper-u_lower)*complex_polynomial_value(
        u_derivative, u_interval
    ).absolute_upper()/2
    residual += entire_tail()
    derivative = complex_polynomial_value(radial_u_coefficients(point, 1), u_interval)
    derivative_lower = max(coordinate_absolute_lower(derivative.real),
                           coordinate_absolute_lower(derivative.imag))
    derivative_lower -= T_UPPER*entire_tail(1)
    assert derivative_lower > 0
    second_at_centre = complex_polynomial_value(
        radial_u_coefficients(point, 2), u_interval
    ).absolute_upper()+T_UPPER**2*entire_tail(2)
    third_bound = sum(
        k*(k-1)*(k-2)*angular_integral_coefficient(0, k, upper)
        *local_radius**(k-3) for k in range(3, DEGREE+1)
    )+T_UPPER**3*entire_tail(3)
    second_bound = second_at_centre+radius*third_bound
    error = residual+radius*radius*second_bound/2
    linear = radius*derivative_lower
    assert error < linear, ("root tube", candidate, float(error/linear))
    rectangle = CI(I(real-radius, real+radius), I(imag-radius, imag+radius))
    x = rectangle*CI(I(1/upper, 1/lower), I.point(0))
    first = normalised_bessel(0, x, TAIL_ARGUMENT)
    next_bessel = normalised_bessel(1, x, TAIL_ARGUMENT)
    kappa = -x*next_bessel/(CI.point(2)*first)
    assert kappa.real.upper < 0, ("kappa real", candidate)
    kappa_real = -kappa.real.upper
    kappa_imag = kappa.imag.absolute_upper()
    assert 63*kappa_real*kappa_real > kappa_imag*kappa_imag, (
        "kappa cone", candidate, float(kappa_imag/kappa_real)
    )
    quotient = -CI.point(2)*first/next_bessel
    assert quotient.real.upper < 0, ("quotient real", candidate)
    quotient_real = -quotient.real.upper
    quotient_imag = quotient.imag.absolute_upper()
    assert 8*quotient_real*quotient_real > quotient_imag*quotient_imag, (
        "quotient cone", candidate
    )
    assert (quotient*quotient).real.upper < 0, ("squared quotient", candidate)
    return {
        "relative_root_margin": (linear-error)/linear,
        "kappa_cone_ratio": kappa_imag*kappa_imag/(kappa_real*kappa_real),
        "quotient_cone_ratio": quotient_imag*quotient_imag/(quotient_real*quotient_real),
    }


def main():
    data = json.loads((HERE/"hessian_compact_centres.json").read_text(encoding="utf-8"))
    candidates = data["centres"]
    if "--probe" in sys.argv:
        chosen = sorted({0, 1, len(candidates)//4, len(candidates)//2,
                         3*len(candidates)//4, len(candidates)-1})
        for index in chosen:
            result = verify_tube(candidates[index])
            print("Probe", index, float(result["relative_root_margin"]), flush=True)
        print("Selected probes passed. No full certificate generated.")
        return
    boundary = Q(21, 10)
    results = []
    for index, candidate in enumerate(candidates):
        assert Q(candidate["t_lower"]) == boundary
        assert boundary < Q(candidate["t_upper"]) <= FINAL_T
        results.append(verify_tube(candidate))
        boundary = Q(candidate["t_upper"])
        if (index+1)%25 == 0:
            print("Verified compact cone tubes", index+1, flush=True)
    assert boundary == FINAL_T
    assert (Q(21, 10)+Q(10, 21))**2*104 > 26**2
    assert Q(13, 12) > Q(26, 25)**2
    assert 1280 > 35**2
    lower_factor = (1-Q(6, 35))*(1-Q(31, 35))/(1+Q(31, 35))**4
    upper_factor = (1-Q(31, 35))**(-2)
    assert 2*lower_factor/upper_factor**2 > Q(1, 400000)
    result = {
        "status": "passed",
        "arithmetic": "fractions.Fraction only",
        "sqrt_contrast_interval": ["21/10", str(FINAL_T)],
        "contrast_interval": ["441/100", str(FINAL_T**2)],
        "parameter_tubes": len(candidates),
        "minimum_relative_root_margin_lower": str(rounded_bound(min(r["relative_root_margin"] for r in results))),
        "maximum_kappa_squared_cone_ratio_upper": str(rounded_bound(max(r["kappa_cone_ratio"] for r in results), True)),
        "maximum_quotient_squared_cone_ratio_upper": str(rounded_bound(max(r["quotient_cone_ratio"] for r in results), True)),
        "proved_cones": ["-Re(kappa)>|kappa|/8", "-Re(x/kappa)>|x/kappa|/3",
                          "Re((x/kappa)^2)<0"],
    }
    (HERE/"hessian_compact_certificate.json").write_text(
        json.dumps(result, indent=2)+"\n", encoding="utf-8"
    )
    print(json.dumps({
        "status": "passed", "parameter_tubes": len(candidates),
        "contrast_upper": float(FINAL_T**2),
        "minimum_relative_root_margin": float(min(r["relative_root_margin"] for r in results)),
        "maximum_kappa_squared_cone_ratio": float(max(r["kappa_cone_ratio"] for r in results)),
        "maximum_quotient_squared_cone_ratio": float(max(r["quotient_cone_ratio"] for r in results)),
    }, indent=2))


if __name__ == "__main__":
    main()

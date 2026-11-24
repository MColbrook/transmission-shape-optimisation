"""Exact Rayleigh sums and rational margins for the disk product cutoff."""

from fractions import Fraction
import json
from math import factorial
from pathlib import Path


def rising(value, length):
    result = 1
    for index in range(length):
        result *= value + index
    return result


def rayleigh_sums(order, count):
    elementary = [
        Fraction(1, 4**degree * factorial(degree) * rising(order + 1, degree))
        for degree in range(count + 1)
    ]
    powers = [Fraction(0)]
    for degree in range(1, count + 1):
        total = sum(
            (-1)**(index - 1) * elementary[index] * powers[degree - index]
            for index in range(1, degree)
        )
        total += (-1)**(degree - 1) * degree * elementary[degree]
        powers.append(total)
    return powers


def main():
    first = rayleigh_sums(1, 3)[3]
    second = rayleigh_sums(2, 4)[4]
    assert first == Fraction(1, 3072)
    assert second == Fraction(7, 3317760)
    angular_dirichlet_margin = 1 - Fraction(29, 2)**3 * first
    angular_buckling_margin = 1 - 26**4 * second
    assert angular_dirichlet_margin > 0
    assert angular_buckling_margin > 0
    angular_plate = Fraction(29, 2) * 26
    assert angular_plate == 377
    assert 30**2 > angular_plate
    product = 104 * angular_plate
    assert product == 39208
    assert product > 14**4
    result = {
        "status": "passed exact rational arithmetic",
        "rayleigh_sum_j1_power6": str(first),
        "rayleigh_sum_j2_power8": str(second),
        "dirichlet_comparison_margin": str(angular_dirichlet_margin),
        "buckling_comparison_margin": str(angular_buckling_margin),
        "angular_plate_lower_bound": str(angular_plate),
        "disk_plate_product_lower_bound": str(product),
        "cutoff_fourth_power": 14**4,
        "product_margin": str(product - 14**4),
        "imported_analytic_inputs": {
            "first_disk_plate_lower_bound": 104,
            "second_radial_dirichlet_lower_bound": 30,
            "source": "Manuscript Section 4.4 for the first disk plate value and Section 2.1 for the second radial Dirichlet value",
            "scope": "The checker verifies the consequences of these separately justified bounds.",
        },
    }
    destination = Path(__file__).with_name("plate_root_count_certificate.json")
    destination.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

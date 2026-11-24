"""Exact coarse Bessel inequalities for full subresonant continuation."""

from fractions import Fraction
import json
from math import factorial
from pathlib import Path


Q = Fraction


def partial(order, argument, degree):
    return sum(
        (-argument / 4)**j * factorial(order)
        / (factorial(j) * factorial(order + j))
        for j in range(degree + 1)
    )


def main():
    alpha_upper = partial(0, Q(6), 4)
    zeta_lower = partial(0, Q(30), 8)
    zeta_upper = partial(0, Q(31), 9)
    first_angular_upper = partial(1, Q(15), 4)
    assert alpha_upper < 0
    assert zeta_lower < 0
    assert zeta_upper > 0
    assert first_angular_upper < 0
    # All tails alternate with decreasing terms beyond these partial sums.
    assert Q(6, 4) < 6**2
    assert Q(30, 4) < 10**2
    assert Q(31, 4) < 11**2
    assert Q(15, 4) < 6 * 7
    # The first angular squared zero is below 15. Its second exceeds 36.
    remainder_sum = Q(1, 192) - Q(1, 15**2)
    assert remainder_sum > 0
    assert remainder_sum < Q(1, 36**2)
    assert 6 * 31 < 14**2
    result = {
        "status": "passed",
        "arithmetic": "fractions.Fraction only",
        "J0_squared_argument_6_even_partial": str(alpha_upper),
        "J0_squared_argument_30_even_partial": str(zeta_lower),
        "J0_squared_argument_31_odd_partial": str(zeta_upper),
        "J1_squared_argument_15_even_partial": str(first_angular_upper),
        "second_angular_inverse_square_upper": str(remainder_sum),
        "consequences_using_Bessel_interlacing": [
            "alpha<6", "30<zeta<31", "third J0 squared zero>36",
            "alpha*zeta<186<196",
        ],
    }
    Path(__file__).with_name("coarse_constants.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

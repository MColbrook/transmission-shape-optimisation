"""Exact normalization bounds for the all-contrast Robin parametrization.

Uses the already audited directed integer interval arithmetic and Bessel
series.  No floating point value enters an assertion.
"""
from fractions import Fraction as F
from pathlib import Path
import sys
import json

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'shape_optimization'))
from calibration_interval import I, Q, j0, hj1, XL, ROOT_DEN


def certify(left, right, sign, threshold):
    pending = [(F(left), F(right), 0)]
    cells = []
    while pending:
        a, b, depth = pending.pop()
        bound = sign * (j0(I(a, b)) + hj1(I(a, b)))
        if bound.lo > threshold * Q:
            cells.append((a, b, bound.lo))
            continue
        assert depth < 24
        m = (a + b) / 2
        pending.extend([(a, m, depth + 1), (m, b, depth + 1)])
    return {
        'interval': [str(F(left)), str(F(right))],
        'sign': sign,
        'strict_lower_bound': str(threshold),
        'cells': len(cells),
        'partition': [[str(a), str(b)] for a, b, _ in sorted(cells)],
        'minimum_enclosure_lower': str(F(min(c for _, _, c in cells), Q)),
    }


def main():
    rows = [certify(0, 6, 1, F(1, 2)), certify(14, 31, -1, F(1, 4))]
    alpha = I(F(XL, ROOT_DEN), F(XL + 1, ROOT_DEN))
    phi_denominator = hj1(alpha)
    assert phi_denominator.lo > Q
    result = {
        'status': 'all_exact_assertions_passed',
        'expression': 'J0(sqrt(z)) + sqrt(z)*J1(sqrt(z))',
        'rows': rows,
        'phi_denominator_lower': str(F(phi_denominator.lo, Q)),
        'consequence': 'For Robin roots 0<=x<6,14<y<31: |eta/J0(sqrt(x))|<2, |eta/J0(sqrt(y))|<4; y-x>8 and sqrt(y)<6 imply |ubar^(j)(r)|<=6^j.',
    }
    Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': result['status'], 'cells': [r['cells'] for r in rows]}))


if __name__ == '__main__':
    main()

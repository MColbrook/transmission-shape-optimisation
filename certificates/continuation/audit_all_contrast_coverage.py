"""Check exact coverage of the radial calibration charts.

This checks the rational partitions, including the inherited regular part
of the resolved boundary chart B. It does not replace the sign verifiers
or the analytic remainder proofs. Python standard library only.
"""
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]


def rectangle_partition(rows, box):
    a0, b0, c0, d0 = map(F, box)
    rectangles = [tuple(map(F, row)) for row in rows]
    events = defaultdict(list)
    total = F(0)
    for index, (a, b, c, d) in enumerate(rectangles):
        assert a0 <= a < b <= b0 and c0 <= c < d <= d0
        events[a].append((1, index))
        events[b].append((-1, index))
        total += (b-a)*(d-c)
    xs = sorted(events)
    assert xs[0] == a0 and xs[-1] == b0
    assert total == (b0-a0)*(d0-c0)
    active = set()
    for j, x in enumerate(xs):
        for sign, index in events[x]:
            if sign == 1:
                assert index not in active
                active.add(index)
            else:
                assert index in active
                active.remove(index)
        if j == len(xs)-1:
            assert not active
            break
        cursor = c0
        for c, d in sorted((rectangles[i][2], rectangles[i][3]) for i in active):
            assert c == cursor, ('gap or overlap', x, c, cursor)
            cursor = d
        assert cursor == d0
    return {'rectangles': len(rectangles), 'area': str(total),
            'sweep_slabs': len(xs)-1, 'no_gaps_or_overlaps': True}


def main():
    specifications = [
        ('center', 'minimizer_literature/continuation/all_contrast_center_reduced.json',
         ('0', '1', '0', '9/16')),
        ('A', 'shape_optimization/continuation/all_contrast_endpoint_A_M48.json',
         ('0', '1/4', '0', '1')),
        ('E', 'shape_optimization/continuation/all_contrast_endpoint_fast_E.json',
         ('1/4', '1', '0', '1/4')),
        ('Esmall', 'shape_optimization/continuation/all_contrast_endpoint_fast_Esmall.json',
         ('1/16', '1/4', '0', '1/4')),
        ('B', 'shape_optimization/continuation/all_contrast_endpoint_B_covered.json',
         ('0', '1/4', '0', '1')),
    ]
    results = {}
    for name, relative, box in specifications:
        path = ROOT / relative
        data = json.loads(path.read_text())
        rows = data['partition']
        if name == 'center':
            assert data['status'] == 'exact_center_PSD_assertions_passed'
            assert all(F(pair[0]) > 0 for pair in data['certified_bounds'])
        if name == 'B':
            assert data['status'] == 'all_exact_resolved_B_assertions_passed'
            assert F(data['eta_cutoff']) == F(1, 16)
            for row in data['inherited_partition']:
                a, b, c, d = map(F, row)
                assert a*c >= F(1, 16) and b*d <= F(1, 4)
            rows = rows + data['inherited_partition']
            assert set(data['inherited_certificates']) == {
                'all_contrast_endpoint_fast_E.json',
                'all_contrast_endpoint_fast_Esmall.json',
            }
        results[name] = rectangle_partition(rows, box)
        results[name]['source'] = relative
    result = {
        'status': 'all_charts_cover_exactly',
        'arithmetic': 'exact rational sweep; no floating point',
        'charts': results,
        'global_cover': (
            'Center covers 0<=r<=3/4 for every eta. For s=1-r<=1/4, '
            'E and Esmall cover eta>=1/16; the original A chart covers '
            '0<=s<=eta<=1/4; B covers eta<=s via delta=eta/s. '
            'In B, eta>=1/16 is inherited from E/Esmall, and the '
            'remaining subset is sign-verified using capped Robin roots. '
            'The endpoint limits are included by the analytic chart factors.'
        ),
        'scope': 'Coverage of the calibration charts. Sign bounds are checked separately.',
    }
    output = Path(__file__).with_name('all_contrast_coverage_audit.json')
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()

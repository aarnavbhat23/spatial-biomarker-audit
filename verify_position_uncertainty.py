"""Independent exhaustive replay of saved discrete geometry certificates.

Uses scalar edge lists and rational contrasts, without importing bound/oracle code.
Checks every saved small-case graph, all admissible categorical assignments,
component count bounds and endpoint witnesses. Scaling claims are summaries only.
"""
from fractions import Fraction
from itertools import combinations, product
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent


def edges_at(positions, selected):
    xy = [positions[i][p] for i,p in enumerate(selected)]
    return tuple((i,j) for i,j in combinations(range(len(xy)),2)
                 if sum((a-b)**2 for a,b in zip(xy[i],xy[j])) <= 25)


def counts(edges, labels):
    a = sum(int(labels[i] == 1)+int(labels[j] == 1) for i,j in edges)
    b = sum(int(labels[i] == 2)+int(labels[j] == 2) for i,j in edges)
    c = sum(int(labels[i] == 1 and labels[j] == 2)+int(labels[i] == 2 and labels[j] == 1) for i,j in edges)
    m = 2*len(edges)
    return m,a,b,c


def contrast(edges, labels):
    m,a,b,c = counts(edges,labels)
    return Fraction(m*c-a*b,m*m) if m else None


def main():
    report = json.loads((ROOT/'results/position_uncertainty.json').read_text())
    assert report['protocol_commit'] == 'b96dac4'
    assert len(report['small']) == len(report['scale']) == 12
    assert {(r['n'],r['seed'],r['budget']) for r in report['small']} == set(product((5,6),(19,43,71),(0,1)))
    assert {(r['n'],r['seed'],r['budget']) for r in report['scale']} == {(n,s,b) for n in (128,512,2048) for s in (101,202) for b in (0,int(np.ceil(.005*n)))}
    worlds = 0
    for r in report['small']:
        n = r['n']; raw = r['labels']; positions = r['positions']; budget = r['budget']
        labels = [x for x in product(range(3),repeat=n)
                  if sum((a == 1) != (b == 1) for a,b in zip(x,raw)) <= budget
                  and sum((a == 2) != (b == 2) for a,b in zip(x,raw)) <= budget]
        graphs = {edges_at(positions,s) for s in product(range(3),repeat=n)}
        raw_edges = set(edges_at(positions,[0]*n)); low = high = None; bound = r['bound']
        for edges in graphs:
            deleted, added = 2*len(raw_edges-set(edges)),2*len(set(edges)-raw_edges)
            for caps in bound['edit_caps'].values():
                assert deleted <= caps[0] and added <= caps[1]
            for labeling in labels:
                m,a,b,c = counts(edges,labeling); worlds += 1
                for actual,key in ((m,'directed_edge_count_bounds'),(a,'a_degree_bounds'),(b,'b_degree_bounds'),(c,'contact_count_bounds')):
                    assert bound[key][0] <= actual <= bound[key][1],(r['n'],r['seed'],key,actual)
                value = contrast(edges,labeling)
                if value is None:
                    assert bound['may_include_empty_graph']; continue
                low = value if low is None else min(low,value)
                high = value if high is None else max(high,value)
                for key in ('outer_bounds','moment_bounds','local_edit_bounds','independent_edge_bounds'):
                    assert bound[key][0]-1e-12 <= value <= bound[key][1]+1e-12
        oracle = r['exact']
        assert oracle['unique_graphs'] == len(graphs)
        assert oracle['empty_graphs'] == int(() in graphs)
        assert oracle['admissible_labelings'] == len(labels)
        assert oracle['position_assignments'] == 3**n
        np.testing.assert_allclose(oracle['bounds'],[float(low),float(high)],atol=1e-15,rtol=0)
        for witness,expected in zip(oracle['witnesses'],(low,high)):
            assert tuple(witness['labels']) in labels
            assert all(0 <= p < 3 for p in witness['positions'])
            actual = contrast(edges_at(positions,witness['positions']),witness['labels'])
            assert actual == expected == Fraction(witness['numerator'],witness['denominator'])
    improvements = []
    for r in report['scale']:
        b = r['bound']; width = b['outer_bounds'][1]-b['outer_bounds'][0]
        local_width = b['local_edit_bounds'][1]-b['local_edit_bounds'][0]
        reduction = 1-width/local_width if local_width else 0.
        np.testing.assert_allclose(reduction,r['relative_width_reduction'],atol=1e-15,rtol=0)
        improvements.append(reduction)
    median = float(np.median(improvements))
    np.testing.assert_allclose(median,report['median_relative_width_reduction'],atol=1e-15,rtol=0)
    assert report['continuation_gate_passed'] == (median >= .1)
    output = {'verified_small_cases': 12,'distinct_graph_label_worlds_checked': worlds,
              'component_bounds_and_rational_witnesses_verified': True,
              'median_relative_width_reduction': median,'continuation_gate_passed': median >= .1,
              'limitations': 'Finite discrete contract only. Scaling summaries checked arithmetically, not exhaustively; no biological coverage calibration.'}
    (ROOT/'results/position_uncertainty_verification.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))


if __name__ == '__main__': main()

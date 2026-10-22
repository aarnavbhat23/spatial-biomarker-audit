import itertools
import unittest
import numpy as np
from categorical_identifiability import categorical_witness,repair_overlap,review_pair_bound
from contrast_count_slices import exact_numerator

STATES=np.array([[0,0],[1,0],[0,1]])


class CategoricalTests(unittest.TestCase):
    def test_one_move_against_independent_enumeration(self):
        edges=np.array([[0,1],[0,2],[1,2],[2,3],[3,1]])
        for categories in itertools.product(range(3),repeat=4):
            raw=STATES[list(categories)]
            possible=[exact_numerator(edges,raw)[0]]
            for i,state in itertools.product(range(4),STATES):
                candidate=raw.copy();candidate[i]=state
                if np.all(np.sum(candidate!=raw,axis=0)<=1):possible.append(exact_numerator(edges,candidate)[0])
            for maximize in (False,True):
                answer=categorical_witness(edges,raw,1,1,maximize,max_moves=1)
                self.assertEqual(answer['numerator'],(max if maximize else min)(possible))

    def test_feasible_multistep(self):
        rng=np.random.default_rng(29);n=20;edges=np.array([(i,(i+j)%n) for i in range(n) for j in (1,2,4)])
        raw=STATES[rng.integers(0,3,n)]
        for maximize in (False,True):
            w=categorical_witness(edges,raw,3,3,maximize);x=raw.copy();x[w['a_flips'],0]^=1;x[w['b_flips'],1]^=1
            self.assertTrue(np.all(x.sum(axis=1)<=1));self.assertTrue(np.all(np.sum(x!=raw,axis=0)<=3))
            self.assertEqual(exact_numerator(edges,x),(w['numerator'],w['denominator']))

    def test_overlap_repair(self):
        raw=STATES.copy();bad=np.ones((3,2),int);fixed=repair_overlap(raw,bad)
        self.assertTrue(np.all(fixed.sum(axis=1)<=1));self.assertTrue(np.all(np.sum(fixed!=raw,axis=0)<=np.sum(bad!=raw,axis=0)))

    def test_error_balls_exhaustively(self):
        vectors=list(itertools.product(range(3),repeat=3))
        for x,y in itertools.combinations(vectors,2):
            x=np.array(x);y=np.array(y)
            for error in (0,1):
                shared=any(np.sum(np.array(z)!=x)<=error and np.sum(np.array(z)!=y)<=error for z in vectors)
                result=review_pair_bound(STATES[x],STATES[y],np.arange(3),error)
                self.assertEqual(result['distinguishable_even_after_all_cells'],not shared)

    def test_field_cover_against_all_subsets(self):
        rng=np.random.default_rng(2119)
        for _ in range(30):
            x=STATES[rng.integers(0,3,12)];y=STATES[rng.integers(0,3,12)];fields=rng.integers(0,4,12)
            for error in (0,1,2):
                result=review_pair_bound(x,y,fields,error);different=np.any(x!=y,axis=1);costs=[]
                for subset in itertools.product((False,True),repeat=4):
                    selected=np.flatnonzero(subset)
                    if np.sum(different & np.isin(fields,selected))>=2*error+1:costs.append(len(selected))
                self.assertEqual(result['minimum_unit_cost_fields_for_pair'],min(costs) if costs else None)

    def test_invalid_identity(self):
        with self.assertRaises(ValueError):categorical_witness(np.array([[0,1]]),np.ones((2,2)),1,1,True)

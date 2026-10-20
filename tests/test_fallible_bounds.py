import itertools,unittest
import numpy as np
from fallible_bounds import linear_optimum,review_bounds,flippable


class FallibleBoundsTests(unittest.TestCase):
    def test_signed_oracle_and_coordinate_feasibility_match_exhaustive_worlds(self):
        rng=np.random.default_rng(19);worlds=np.array(list(itertools.product((0,1),repeat=5)))
        for _ in range(180):
            raw=rng.integers(0,2,5);ref=rng.integers(-1,2,5);budget=int(rng.integers(0,6));errors=int(rng.integers(0,3))
            valid=(np.sum(worlds!=raw,axis=1)<=budget)&(np.sum((worlds!=ref)&(ref>=0),axis=1)<=errors)
            weights=rng.normal(size=5)
            for maximize in (False,True):
                result=linear_optimum(weights,raw,ref,budget,errors,maximize)
                if not valid.any():self.assertIsNone(result);continue
                value,labels=result
                self.assertAlmostEqual(value,(max if maximize else min)((worlds[valid]@weights).tolist()))
                self.assertTrue(np.any(np.all(worlds[valid]==labels,axis=1)))
            if valid.any():
                anchor=np.where(ref>=0,ref,raw);possible=flippable(raw,ref,budget,errors,anchor)
                np.testing.assert_equal(possible,np.any(worlds[valid]!=anchor,axis=0))

    def test_graph_outer_bounds_and_witnesses_match_enumeration_with_fallible_reviews(self):
        rng=np.random.default_rng(29);n=4
        worlds=np.array(list(itertools.product((0,1),repeat=2*n))).reshape(-1,n,2)
        edges=np.array([(0,1),(0,2),(1,3),(2,0),(3,1),(3,2)])
        values=np.mean(worlds[:,edges[:,0],0]*worlds[:,edges[:,1],1],axis=1)
        for _ in range(120):
            raw=rng.integers(0,2,(n,2));ref=rng.integers(-1,2,(n,2));budgets=rng.integers(0,n+1,2);errors=rng.integers(0,3,2)
            valid=np.all(np.sum(worlds!=raw,axis=1)<=budgets,axis=1)&np.all(np.sum((worlds!=ref)&(ref>=0),axis=1)<=errors,axis=1)
            result=review_bounds(edges,raw,ref,budgets,errors)
            self.assertEqual(result['empty'],not valid.any())
            if not valid.any():continue
            self.assertLessEqual(result['outer_bounds'][0],values[valid].min()+1e-12)
            self.assertGreaterEqual(result['outer_bounds'][1],values[valid].max()-1e-12)
            for index,detail in enumerate(result['details']):
                labels=raw.copy();labels[detail['a_flips'],0]^=1;labels[detail['b_flips'],1]^=1
                self.assertTrue(np.any(np.all(worlds[valid]==labels,axis=(1,2))))
                self.assertAlmostEqual(np.mean(labels[edges[:,0],0]*labels[edges[:,1],1]),result['attained_bounds'][index])

    def test_systematic_reference_floor_and_empty_intersection_are_retained(self):
        raw=np.ones((3,2),int);edges=np.array([(0,1),(1,2),(2,0)])
        result=review_bounds(edges,raw,raw,(1,1),(1,1))
        self.assertGreater(result['attained_bounds'][1]-result['attained_bounds'][0],0.)
        self.assertTrue(review_bounds(edges,raw,1-raw,(0,0),(0,0))['empty'])
        exact=review_bounds(edges,raw,1-raw,(3,3),(0,0))
        self.assertEqual(exact['outer_bounds'],[0.,0.])


if __name__=='__main__':unittest.main()

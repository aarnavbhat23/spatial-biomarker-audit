import unittest
import numpy as np
from fallible_policy import fields_from_coordinates,run_policy


class FalliblePolicyTests(unittest.TestCase):
    def test_every_policy_keeps_truth_and_monotone_bounds_after_fallible_reviews(self):
        edges=np.array([(0,1),(0,2),(1,3),(2,4),(3,5),(4,0),(5,2)])
        truth=np.array([[1,0],[0,1],[1,0],[1,1],[0,1],[0,0]])
        raw=truth.copy();raw[[0,3],0]^=1;raw[[1,4],1]^=1
        reference=truth.copy();reference[2,0]^=1;reference[5,1]^=1;fields=np.arange(6)
        target=np.mean(truth[edges[:,0],0]*truth[edges[:,1],1])
        for policy in ('uniform','degree','influence','lookahead','witness'):
            queried=[]
            def review(ids):queried.extend(ids.tolist());return reference[ids]
            result=run_policy(edges,raw,(2,2),(1,1),fields,review,policy,17,target_width=0.)
            self.assertFalse(result['empty']);self.assertEqual(len(queried),len(set(queried)))
            previous=[0.,1.]
            for item in result['trace']:
                lo,hi=item['outer_bounds'];self.assertLessEqual(lo,target);self.assertGreaterEqual(hi,target)
                self.assertGreaterEqual(lo,previous[0]);self.assertLessEqual(hi,previous[1]);previous=[lo,hi]
            self.assertFalse(result['attained']) # Systematic reference errors leave a floor.

    def test_field_partition_includes_upper_bound_and_rejects_bad_observations(self):
        fields=fields_from_coordinates(np.array([[0,0],[1,1],[.3,.6]]));np.testing.assert_equal(fields,[0,15,9])
        edges=np.array([(0,1),(1,0)]);raw=np.zeros((2,2),int)
        with self.assertRaises(ValueError):
            run_policy(edges,raw,(1,1),(0,0),np.arange(2),lambda ids:np.full((len(ids),2),-1),'uniform',1,target_width=0.)


if __name__=='__main__':unittest.main()

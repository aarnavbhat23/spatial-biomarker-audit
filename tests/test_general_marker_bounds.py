import itertools,unittest
import numpy as np
from general_marker_bounds import fast_bounds,optimize_bounds


class GeneralMarkerTests(unittest.TestCase):
    def test_full_enumeration_validates_outer_bounds_witnesses_and_integer_optima(self):
        edges=np.array([(0,1),(0,2),(1,0),(1,2),(2,0)])
        worlds=np.array(list(itertools.product((0,1),repeat=6))).reshape(-1,3,2)
        targets=np.mean(worlds[:,edges[:,0],0]*worlds[:,edges[:,1],1],axis=1)
        for i,recorded in enumerate(worlds):
            for ba,bb in ((0,0),(1,1),(2,1),(2,2),(3,3)):
                valid=(np.sum(worlds[:,:,0]!=recorded[:,0],axis=1)<=ba)&(np.sum(worlds[:,:,1]!=recorded[:,1],axis=1)<=bb)
                exact=[targets[valid].min(),targets[valid].max()];bound=fast_bounds(edges,recorded,ba,bb)
                self.assertLessEqual(bound['outer_bounds'][0],exact[0]+1e-12)
                self.assertGreaterEqual(bound['outer_bounds'][1],exact[1]-1e-12)
                for j,detail in enumerate(bound['details']):
                    changed=recorded.copy();changed[detail['a_flips'],0]^=1;changed[detail['b_flips'],1]^=1
                    attained=np.mean(changed[edges[:,0],0]*changed[edges[:,1],1])
                    self.assertAlmostEqual(attained,bound['attained_bounds'][j])
                if i%8==0:
                    integer=optimize_bounds(edges,recorded,ba,bb)
                    np.testing.assert_allclose(integer['attained_bounds'],exact,atol=1e-12)
                    self.assertTrue(all(x['optimal'] for x in integer['objectives']))

    def test_fallible_reference_feasible_set_matches_enumeration(self):
        edges=np.array([(0,1),(1,2),(2,0)]);raw=np.array([[1,0],[1,1],[0,0]]);ref=np.array([[0,0],[-1,-1],[1,1]])
        worlds=np.array(list(itertools.product((0,1),repeat=6))).reshape(-1,3,2)
        valid=np.ones(len(worlds),bool)
        for k in (0,1):
            valid&=np.sum(worlds[:,:,k]!=raw[:,k],axis=1)<=2
            valid&=np.sum((worlds[:,:,k]!=ref[:,k])&(ref[:,k]>=0),axis=1)<=1
        values=np.mean(worlds[:,edges[:,0],0]*worlds[:,edges[:,1],1],axis=1)
        answer=optimize_bounds(edges,raw,2,2,reference=ref,reference_budgets=(1,1))
        np.testing.assert_allclose(answer['attained_bounds'],[values[valid].min(),values[valid].max()])

    def test_timeout_uses_dual_bound_instead_of_suboptimal_incumbent(self):
        from types import SimpleNamespace
        from unittest.mock import patch
        edges=np.array([(0,1),(1,2),(2,0)]);calls=np.zeros((3,2),int)
        timeout=SimpleNamespace(status=1,success=False,fun=0.,mip_dual_bound=-1.,x=np.zeros(9))
        with patch('general_marker_bounds.milp',return_value=timeout):
            result=optimize_bounds(edges,calls,1,1)
        self.assertEqual(result['attained_bounds'],[0.,0.])
        self.assertGreaterEqual(result['outer_bounds'][1],1/3)
        self.assertFalse(result['objectives'][1]['optimal'])

    def test_lp_is_outer_bound_and_conflicting_exact_reference_is_infeasible(self):
        edges=np.array([(0,1),(1,2),(2,0)]);calls=np.zeros((3,2),int)
        lp=optimize_bounds(edges,calls,1,1,relax=True);integer=optimize_bounds(edges,calls,1,1)
        self.assertLessEqual(lp['outer_bounds'][0],integer['attained_bounds'][0])
        self.assertGreaterEqual(lp['outer_bounds'][1],integer['attained_bounds'][1])
        self.assertIsNone(lp['attained_bounds'])
        self.assertTrue(optimize_bounds(edges,calls,0,0,reference=np.ones((3,2),int))['empty'])


if __name__=='__main__':unittest.main()

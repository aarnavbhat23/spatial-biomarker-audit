"""Joint importance sampling, full-support checks and complementarity."""
import itertools
import unittest
import numpy as np
from study_multiclaim import envelope_distribution,choose_distribution,multiclaim_trial
from study_raster_acquisition import prepare_frame_pairs
from study_cell_reuse import resolved_edges
from study_residual_acquisition import residual_wealth


class MulticlaimTests(unittest.TestCase):
    def test_envelope_minimax_and_conditional_betting_expectation(self):
        owner=np.array([0,0,1,2]);masks=np.array([[1,1,0,0],[0,1,1,1]],bool)
        totals=masks.sum(axis=1);remaining=np.stack([np.bincount(owner[m],minlength=3) for m in masks])
        p,h=envelope_distribution(remaining,totals,np.ones(2,bool))
        np.testing.assert_allclose(p,[.6,.2,.2]);self.assertAlmostEqual(h,5/3)
        for first in np.linspace(.01,.98,30):
            for second in np.linspace(.01,.99-first,20):
                candidate=np.array([first,second,1-first-second])
                self.assertGreaterEqual(np.max(remaining/totals[:,None]/candidate),h-1e-12)
        for bits in itertools.product((0,1),repeat=4):
            outcome=np.array(bits);truth=np.sum(masks*outcome,axis=1)/totals
            group_positive=np.stack([np.bincount(owner,weights=m*outcome,minlength=3) for m in masks])
            z=group_positive/(totals[:,None]*p)
            np.testing.assert_allclose(z@p,truth,atol=1e-13)
            for j in range(2):
                for positive in (False,True):
                    expectation=0.
                    for a,b in itertools.product(range(3),repeat=2):
                        w=np.exp(residual_wealth(truth[j],[z[j,a],z[j,b]],[0,0],[h,h],positive))
                        expectation+=p[a]*p[b]*w
                    self.assertAlmostEqual(expectation,1.,places=12)

    def test_focused_reuse_requires_full_support(self):
        remaining=np.array([[1,0],[1,1]]);totals=np.array([1,2]);active=np.ones(2,bool)
        _,update,_=choose_distribution(remaining,totals,active,np.array([1.,.5]),'focused_reuse')
        np.testing.assert_array_equal(update,[True,False])
        _,update,_=choose_distribution(remaining,totals,active,np.array([.5,1.]),'focused_reuse')
        np.testing.assert_array_equal(update,[True,True])

    def test_resolution_is_not_submodular(self):
        edges=np.array([[0,1]]);labels=np.ones((2,2),bool)
        empty=resolved_edges(edges,labels,np.array([0,0],bool)).sum()
        first=resolved_edges(edges,labels,np.array([1,0],bool)).sum()
        second=resolved_edges(edges,labels,np.array([0,1],bool)).sum()
        both=resolved_edges(edges,labels,np.array([1,1],bool)).sum()
        self.assertEqual(first-empty,0);self.assertEqual(both-second,1)

    def test_all_policies_reach_joint_exact_census_small_graph(self):
        coords=np.array([[120,120],[135,135],[200,140]])
        edges=np.array([[0,1],[1,0],[1,2],[2,1]])
        masks=np.array([[1,1,0,0],[1,1,1,1]],bool)
        for bits in itertools.product((False,True),repeat=6):
            data={'coords':coords,'edges':edges,'labels':np.array(bits).reshape(3,2)}
            prepared=prepare_frame_pairs(data,128)
            for policy in ('envelope','mixture','uniform','focused','focused_reuse','greedy'):
                row=multiclaim_trial(data,prepared,masks,128,9,policy,width=0.,batch=2,max_queries=64)
                self.assertTrue(row['stopped']);self.assertTrue(row['joint_covered_all_checkpoints'])
                np.testing.assert_allclose(np.array(row['intervals'])[:,0],row['truth'])
                np.testing.assert_allclose(np.array(row['intervals'])[:,1],row['truth'])


if __name__=='__main__':unittest.main()

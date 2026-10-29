import unittest
import numpy as np
from study_conserved_spatial_mass import proximity,bounds


class ConservedMolecularTests(unittest.TestCase):
    def test_rigid_motion_and_repartition_leave_score_unchanged(self):
        a=np.array([[1.,3.],[2.,7.],[8.,9.]])
        b=np.array([[3.,5.],[9.,2.],[4.,6.]])
        original=proximity(a,b,4.)
        # Regrouping molecule rows into different cell partitions changes neither
        # the original positions nor mass. Both partitions flatten to a permutation.
        regrouped=np.concatenate((a[[2]],a[[0,1]]))
        self.assertAlmostEqual(original,proximity(regrouped,b[[1,2,0]],4.),places=14)
        self.assertAlmostEqual(original,proximity(a+[17,-9],b+[17,-9],4.),places=14)
        self.assertGreater(abs(original-proximity(a+20,b,4.)),.1)

    def test_one_missing_per_channel_bound(self):
        a=np.array([[0.,0.],[1.,1.],[2.,0.]])
        b=np.array([[0.,2.],[1.,0.],[4.,0.],[3.,2.]])
        baseline=proximity(a,b,2.)
        allowed=1/(len(a)+1)+1/(len(b)+1)
        for x in (-100.,0.,1.,100.):
            for y in (-100.,0.,1.,100.):
                full=proximity(np.vstack((a,[x,0.])),np.vstack((b,[y,0.])),2.)
                self.assertLessEqual(abs(full-baseline),allowed)

    def test_zero_radius_and_extreme_admissible_movements(self):
        a=np.array([[0.,0.],[3.,2.]]); b=np.array([[5.,1.],[9.,3.]])
        zero=bounds(a,b,3.,0.)
        np.testing.assert_allclose(zero['pairwise_bounds'],[zero['score']]*2,atol=1.1e-12,rtol=0)
        interval=bounds(a,b,3.,1.)['pairwise_bounds']
        for shift in (np.array([1.,0.]),np.array([0.,-1.])):
            moved=proximity(a+shift,b-shift,3.)
            self.assertLessEqual(interval[0],moved)
            self.assertGreaterEqual(interval[1],moved)


if __name__ == '__main__': unittest.main()

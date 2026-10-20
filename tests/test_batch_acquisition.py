"""Batch martingale, duplicate-query and sharp certificate checks."""
import itertools
import unittest
import numpy as np
from study_batch_acquisition import batch_values,batch_trial
from study_raster_acquisition import prepare_frame_pairs
from study_residual_acquisition import residual_wealth


class BatchTests(unittest.TestCase):
    def test_product_expectation_for_every_small_population(self):
        # Unequal PPS groups, duplicate draws, and nontrivial known positives.
        groups=[np.array([0]),np.array([1,2])];probabilities=[1/3,2/3]
        for bits in itertools.product((0.,1.),repeat=3):
            outcomes=np.array(bits);truth=(1+sum(bits))/5
            for positive in (False,True):
                expectation=0.
                for ids in itertools.product(range(2),repeat=3):
                    values=batch_values(outcomes,[groups[i] for i in ids],1,3,5)
                    wealth=np.exp(residual_wealth(truth,values,[.2]*3,[.8]*3,positive))
                    expectation+=np.prod([probabilities[i] for i in ids])*wealth
                self.assertAlmostEqual(expectation,1.,places=12)

    def test_duplicate_draws_are_retained(self):
        values=batch_values(np.array([1,0,0]),[np.array([0]),np.array([0])],0,3,3)
        self.assertEqual(values,[1.,1.])

    def test_acquisition_and_deterministic_bounds_all_labelings(self):
        coords=np.array([[120,120],[135,135],[200,140]])
        edges=np.array([[0,1],[1,0],[1,2],[2,1]])
        for bits in itertools.product((False,True),repeat=6):
            data={'coords':coords,'edges':edges,'labels':np.array(bits).reshape(3,2)}
            prepared=prepare_frame_pairs(data,128)
            truth=float(np.mean(data['labels'][edges[:,0],0]&data['labels'][edges[:,1],1]))
            for batch,greedy_batch in ((1,1),(4,1),(None,1),(None,4)):
                row=batch_trial(data,prepared,128,7,batch,width=0.,max_draws=64,greedy_batch=greedy_batch)
                self.assertTrue(row['stopped'])
                np.testing.assert_allclose(row['interval'],[truth,truth],atol=1e-14)
                self.assertLessEqual(row['physical_fields'],2)


if __name__=='__main__':unittest.main()

"""Enumerated checks for acquisition, interval inversion and optional stopping."""
import itertools
import unittest

import numpy as np
from scipy.stats import binom

from study_region_scale import acquisition_trial,bounded_sequence,frame,logwealth


class RegionScaleTests(unittest.TestCase):
    def test_census_returns_exact_nonnegative_width(self):
        f={'counts':np.array([89613,89613]),'positive':np.array([512.,512.]),
           'rectangles':np.array([[0,0,10,10],[10,0,20,10]])}
        row=acquisition_trial(f,np.tile([0,1],4),width=.005)
        self.assertTrue(row['stopped'])
        self.assertEqual(row['interval'],[1024/179226]*2)
        self.assertFalse(row['empty_confidence_set'])

    def test_pps_expectation_and_halo_containment(self):
        coords=np.array([[127,127],[130,128],[250,130],[260,131]])
        edges=np.array([[0,1],[1,0],[1,2],[2,1],[2,3]])
        labels=np.array([[1,0],[0,1],[1,1],[0,0]],dtype=bool)
        data={'coords':coords,'edges':edges,'labels':labels}
        f=frame(data,128)
        truth=(labels[edges[:,0],0]&labels[edges[:,1],1]).mean()
        self.assertAlmostEqual(np.dot(f['counts']/f['counts'].sum(),f['positive']/f['counts']),truth)
        for edge,owner in zip(edges,f['owners']):
            x0,y0,x1,y1=f['rectangles'][owner]
            self.assertTrue(np.all((coords[edge,0]>=x0)&(coords[edge,0]<x1)&(coords[edge,1]>=y0)&(coords[edge,1]<y1)))

    def test_mixture_mean_one_and_optional_crossing(self):
        probability=.2;size=10;mean_terminal=[0.,0.];cross_probability=0.
        for path in itertools.product((0.,1.),repeat=size):
            mass=probability**sum(path)*(1-probability)**(size-sum(path))
            log_positive=logwealth(probability,path,True)
            log_negative=logwealth(probability,path,False)
            mean_terminal[0]+=mass*np.exp(log_positive)
            mean_terminal[1]+=mass*np.exp(log_negative)
            crossed=any(max(logwealth(probability,path[:t],True),logwealth(probability,path[:t],False))>=np.log(40) for t in range(1,size+1))
            cross_probability+=mass*crossed
        np.testing.assert_allclose(mean_terminal,[1.,1.],rtol=1e-12)
        self.assertLessEqual(cross_probability,.05+1e-12)

    def test_interval_inversion_exact_binomial_coverage(self):
        probability=.2;size=12;coverage=0.
        for k in range(size+1):
            low,high=bounded_sequence([1.]*k+[0.]*(size-k))
            self.assertLessEqual(low,high)
            coverage+=binom.pmf(k,size,probability)*(low<=probability<=high)
        self.assertGreaterEqual(coverage,.95-1e-12)


if __name__=='__main__':unittest.main()

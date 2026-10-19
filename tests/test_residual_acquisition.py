"""Exact small-population checks under outcome-adaptive scale selection."""
import itertools
import unittest

import numpy as np

from study_residual_acquisition import residual_interval, residual_trial, residual_wealth


class ResidualAcquisitionTests(unittest.TestCase):
    def test_outcome_adaptive_sampling_martingale_by_enumeration(self):
        counts=np.array([1,2,3,4])
        for binary in itertools.product((0,1),repeat=4):
            positive=counts*np.array(binary);truth=positive.sum()/counts.sum()
            terminal=np.zeros(2);crossing=0.
            def visit(known, mass, values, lower, upper, crossed):
                nonlocal crossing
                unknown=~known;remaining=counts[unknown].sum()
                known_positive=positive[known].sum()
                if remaining==0:
                    terminal[:]+=mass*np.exp([residual_wealth(truth,values,lower,upper,sign) for sign in (True,False)])
                    crossing+=mass*crossed
                    return
                # The next partition explicitly depends on already revealed labels.
                groups=np.array([0,0,1,1]) if known_positive>0 else np.arange(4)
                expectation=0.
                choices=[]
                for group in np.unique(groups[unknown]):
                    members=(groups==group)&unknown
                    count=counts[members].sum();probability=count/remaining
                    z=known_positive/10+remaining/10*positive[members].sum()/count
                    expectation+=probability*z
                    choices.append((members,probability,z))
                self.assertAlmostEqual(expectation,truth,places=13)
                for members,probability,z in choices:
                    new_values=values+[z];new_lower=lower+[known_positive/10]
                    new_upper=upper+[(known_positive+remaining)/10]
                    wealth=[residual_wealth(truth,new_values,new_lower,new_upper,s) for s in (True,False)]
                    visit(known|members,mass*probability,new_values,new_lower,new_upper,
                          crossed or max(wealth)>=np.log(40))
            visit(np.zeros(4,dtype=bool),1.,[],[],[],False)
            np.testing.assert_allclose(terminal,[1.,1.],atol=2e-12)
            self.assertLessEqual(crossing,.05+1e-12)

    def test_nested_measurements_finish_without_duplicates(self):
        counts=np.array([1,2,3,4]);positive=np.array([0,1,1,0])
        frames={128:{'owner':np.arange(4),'regions':4,
                     'rectangles':np.array([[0,0,10,10],[10,0,20,10],[20,0,30,10],[30,0,40,10]])},
                2048:{'owner':np.array([0,0,1,1]),'regions':2,
                      'rectangles':np.array([[0,0,20,10],[20,0,40,10]])}}
        for policy in ('128','2048','adaptive'):
            row=residual_trial(counts,positive,frames,10,policy,0.,0.)
            self.assertTrue(row['stopped']);self.assertLessEqual(row['draws'],4)
            self.assertEqual(row['interval'],[.2,.2])
            self.assertEqual(row['area_native_pixels'],400.)
            self.assertTrue(row['covered_all_checkpoints'])

    def test_interval_keeps_current_deterministic_bounds(self):
        values=[.1]*20;lower=[0.]*10+[.02]*10;upper=[1.]*10+[.4]*10
        lo,hi=residual_interval(values,lower,upper,.03,.3)
        self.assertTrue(.03<=lo<=hi<=.3)
        self.assertEqual(residual_interval(values,lower,upper,.1,.1),(.1,.1))


if __name__=='__main__':unittest.main()

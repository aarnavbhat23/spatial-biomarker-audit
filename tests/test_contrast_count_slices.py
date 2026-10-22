import itertools
import unittest
import numpy as np
from contrast_count_slices import exact_numerator,maximize_by_counts

class CountSliceTests(unittest.TestCase):
    def test_count_reference_matches_all_binary_worlds(self):
        n=4;edges=np.array([[i,(i+1)%n] for i in range(n)])
        for raw,ba,bb in [(np.array([[1,0],[1,1],[0,0],[0,1]]),1,1),(np.zeros((n,2),int),2,2)]:
            results=[]
            for bits in itertools.product([0,1],repeat=2*n):
                x=np.array(bits).reshape(n,2)
                if np.sum(x[:,0]!=raw[:,0])<=ba and np.sum(x[:,1]!=raw[:,1])<=bb:
                    num,den=exact_numerator(edges,x);results.append(num/den)
            report=maximize_by_counts(edges,raw,ba,bb,time_limit=2.)
            self.assertTrue(report['all_slices_optimal'])
            self.assertAlmostEqual(report['upper_bound'],min(.25,max(results)+1e-7),places=9)
            for row in report['slices']:
                if row['witness']:
                    w=row['witness'];self.assertLessEqual(w['numerator']/w['denominator'],report['upper_bound'])

    def test_irregular_degree_rejected(self):
        with self.assertRaises(ValueError):maximize_by_counts(np.array([[0,1]]),np.zeros((2,2),int),1,1)

if __name__=='__main__':unittest.main()

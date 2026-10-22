import itertools
import unittest
import numpy as np
from enrichment_bounds import contrast_bounds,coupled_contrast_bounds

class EnrichmentBoundsTests(unittest.TestCase):
    def test_exhaustive_small_label_worlds(self):
        rng=np.random.default_rng(8);n=4
        pairs=np.array(list(itertools.product(range(n),repeat=2)))
        worlds=np.array(list(itertools.product([0,1],repeat=2*n))).reshape(-1,n,2)
        for _ in range(80):
            edges=pairs[rng.choice(len(pairs),size=rng.integers(1,len(pairs)+1),replace=False)]
            raw=rng.integers(0,2,(n,2));ba,bb=rng.integers(0,n+1,2)
            result=contrast_bounds(edges,raw,int(ba),int(bb))
            feasible=worlds[(np.sum(worlds[:,:,0]!=raw[:,0],axis=1)<=ba)&(np.sum(worlds[:,:,1]!=raw[:,1],axis=1)<=bb)]
            u,v=edges.T;m=len(edges);d=np.bincount(u,minlength=n);e=np.bincount(v,minlength=n)
            values=np.mean(feasible[:,u,0]*feasible[:,v,1],axis=1)-(feasible[:,:,0]@d/m)*(feasible[:,:,1]@e/m)
            lo,hi=result['outer_bounds'];self.assertLessEqual(lo,values.min()+1e-12);self.assertGreaterEqual(hi,values.max()-1e-12)
            refined=coupled_contrast_bounds(edges,raw,int(ba),int(bb),base=result)
            low,high=refined['outer_bounds'];self.assertLessEqual(low,values.min()+1e-12);self.assertGreaterEqual(high,values.max()-1e-12)
            self.assertGreaterEqual(low,lo);self.assertLessEqual(high,hi)
            c=np.zeros((n,n));c[u,v]=1/m;c-=np.outer(d,e)/m**2
            np.testing.assert_allclose(result['absolute_norms'],[abs(c).sum(axis=1).max(),abs(c).sum(axis=0).max()],atol=1e-12)
            for w in result['witnesses']:
                self.assertLessEqual(len(w['a_flips']),ba);self.assertLessEqual(len(w['b_flips']),bb)
                labels=raw.copy();labels[w['a_flips'],0]^=1;labels[w['b_flips'],1]^=1
                self.assertAlmostEqual(labels[:,0]@c@labels[:,1],w['value'],places=12)

    def test_complete_graph_with_loops_has_identically_zero_contrast(self):
        edges=np.array(list(itertools.product(range(4),repeat=2)));raw=np.array([[1,0],[1,1],[0,1],[0,0]])
        r=contrast_bounds(edges,raw,4,4)
        self.assertLess(r['outer_bounds'][1]-r['outer_bounds'][0],1e-10)

    def test_zero_budgets_and_invalid_graph(self):
        e=np.array([[0,1],[1,2]]);a=np.array([[1,0],[0,1],[1,1]])
        r=contrast_bounds(e,a,0,0);self.assertLess(r['outer_bounds'][1]-r['outer_bounds'][0],1e-10)
        with self.assertRaises(ValueError):contrast_bounds(np.array([[0,1],[0,1]]),a,1,1)
        with self.assertRaises(ValueError):contrast_bounds(e,a,-1,1)

if __name__=='__main__':unittest.main()

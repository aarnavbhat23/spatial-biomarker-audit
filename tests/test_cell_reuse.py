"""Independent visibility and logical-bound checks for shared cell labels."""
import itertools
import unittest

import numpy as np

from study_cell_reuse import cell_trial, prepare_visibility, resolved_edges


class CellReuseTests(unittest.TestCase):
    def test_spatial_lookup_matches_brute_force_and_census(self):
        coords=np.array([[100,100],[150,100],[200,100],[250,100],[300,100],[350,100]])
        edges=np.array([(i,j) for i in range(6) for j in range(6) if 0<abs(coords[i,0]-coords[j,0])<=80])
        labels=np.array([[1,0],[0,1],[0,0],[1,1],[1,0],[0,1]],dtype=bool)
        data={'coords':coords,'edges':edges,'labels':labels}
        for core in (128,2048):
            prepared=prepare_visibility(data,core)
            for rectangle,ids in zip(prepared['frame']['rectangles'],prepared['visible']):
                brute=np.flatnonzero(np.all((coords>=rectangle[:2])&(coords<rectangle[2:]),axis=1))
                np.testing.assert_array_equal(np.sort(ids),brute)
            for node,ids in enumerate(prepared['incident']):
                np.testing.assert_array_equal(np.sort(ids),np.flatnonzero(np.any(edges==node,axis=1)))
            row=cell_trial(data,prepared,41,0.)
            truth=float((labels[edges[:,0],0]&labels[edges[:,1],1]).mean())
            self.assertTrue(row['stopped']);self.assertEqual(row['interval'],[truth,truth])
            self.assertEqual(row['resolved_edges'],len(edges))

    def test_partial_marker_bounds_are_sharp_by_exhaustive_completion(self):
        edges=np.array([[0,1],[0,2],[1,2],[2,3],[3,0]])
        labels=np.array([[0,1],[1,0],[1,1],[0,1]],dtype=bool)
        for visible in itertools.product((False,True),repeat=4):
            seen=np.array(visible);known=resolved_edges(edges,labels,seen)
            outcomes=labels[edges[:,0],0]&labels[edges[:,1],1]
            lower=int(outcomes[known].sum());upper=lower+int((~known).sum())
            unknown=np.argwhere(np.repeat((~seen)[:,None],2,axis=1))
            totals=[]
            for bits in itertools.product((False,True),repeat=len(unknown)):
                completed=labels.copy()
                for position,value in zip(unknown,bits):completed[tuple(position)]=value
                # Visibility-based decisions cannot use any unobserved marker bit.
                np.testing.assert_array_equal(resolved_edges(edges,completed,seen),known)
                totals.append(int((completed[edges[:,0],0]&completed[edges[:,1],1]).sum()))
            self.assertEqual(min(totals),lower);self.assertEqual(max(totals),upper)


if __name__=='__main__':unittest.main()

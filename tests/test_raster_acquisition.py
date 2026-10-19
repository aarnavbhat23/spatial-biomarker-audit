"""Check camera-field ownership, physical cost and residual expectation."""
import itertools
import unittest

import numpy as np

from study_cell_reuse import cell_trial, resolved_edges
from study_raster_acquisition import prepare_frame_pairs, prepare_raster


class RasterAcquisitionTests(unittest.TestCase):
    def test_diagonal_pair_acquires_two_fields_not_bounding_rectangle(self):
        data={'coords':np.array([[127,127],[129,129]]),'edges':np.array([[0,1],[1,0]]),
              'labels':np.ones((2,2),dtype=bool)}
        prepared=prepare_frame_pairs(data,128)
        row=cell_trial(data,prepared,91,0.)
        self.assertEqual(row['draws'],1)
        self.assertEqual(row['area_native_pixels'],2*128**2)
        self.assertEqual(row['interval'],[1.,1.])

    def test_group_expectation_under_every_visibility_state(self):
        coords=np.array([[110,110],[140,140],[210,140],[250,170]])
        edges=np.array([(i,j) for i in range(4) for j in range(4) if 0<np.linalg.norm(coords[i]-coords[j])<=80])
        labels=np.array([[1,0],[0,1],[1,1],[0,0]],dtype=bool)
        data={'coords':coords,'edges':edges,'labels':labels}
        prepared=prepare_frame_pairs(data,128)
        outcomes=labels[edges[:,0],0]&labels[edges[:,1],1];truth=outcomes.mean();total=len(edges)
        for group,owned in enumerate(prepared['owned']):
            self.assertTrue(np.isin(edges[owned],prepared['visible'][group]).all())
        for state in itertools.product((False,True),repeat=4):
            known=resolved_edges(edges,labels,np.array(state))
            unknown=~known;remaining=unknown.sum()
            if not remaining:continue
            known_positive=outcomes[known].sum();expectation=0.
            for owned in prepared['owned']:
                query=owned[unknown[owned]]
                if not len(query):continue
                value=known_positive/total+remaining/total*outcomes[query].mean()
                expectation+=len(query)/remaining*value
            self.assertAlmostEqual(expectation,truth,places=13)

    def test_raster_visibility_matches_brute_force(self):
        data={'coords':np.array([[100,100],[150,100],[210,100],[270,100]]),
              'edges':np.array([[0,1],[1,0],[1,2],[2,3]]),
              'labels':np.ones((4,2),dtype=bool)}
        prepared=prepare_raster(data,128,256)
        for rectangle,ids in zip(prepared['frame']['rectangles'],prepared['visible']):
            self.assertTrue(np.all(rectangle%256==0))
            expected=np.flatnonzero(np.all((data['coords']>=rectangle[:2])&(data['coords']<rectangle[2:]),axis=1))
            np.testing.assert_array_equal(np.sort(ids),expected)


if __name__=='__main__':unittest.main()

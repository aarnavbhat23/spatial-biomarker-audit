import itertools
import unittest
import numpy as np
from graph_error import pair_frequency_bounds, conformal_upper, graph_errors, decode_mask, match_instances


class GraphErrorTests(unittest.TestCase):
    def test_exhaustive_add_delete_binary_mark_enclosure(self):
        # Every mark pattern, deletion set, and appended mark pattern on up to 5 edges.
        for m in range(1,6):
            for marks in itertools.product([0,1],repeat=m):
                s=sum(marks)
                for d in range(m+1):
                    for a in range(3):
                        lo,hi=pair_frequency_bounds(m,s,s,d,a)
                        for removed in itertools.chain.from_iterable(itertools.combinations(range(m),k) for k in range(d+1)):
                            kept=[marks[i] for i in range(m) if i not in removed]
                            for added_count in range(a+1):
                                for added in itertools.product([0,1],repeat=added_count):
                                    target=kept+list(added)
                                    if target:self.assertTrue(lo-1e-12<=np.mean(target)<=hi+1e-12)

    def test_label_interval_composition(self):
        lo,hi=pair_frequency_bounds(10,2,6,1,2)
        for s in range(2,7):
            a,b=pair_frequency_bounds(10,s,s,1,2)
            self.assertLessEqual(lo,a);self.assertGreaterEqual(hi,b)
        self.assertEqual(pair_frequency_bounds(0,0,0,0,1),(0.,1.))
        with self.assertRaises(ValueError):pair_frequency_bounds(3,0,4,0,0)

    def test_finite_sample_quantile(self):
        self.assertEqual(conformal_upper(np.arange(50)),48.)
        self.assertEqual(conformal_upper([0.]),float('inf'))
        with self.assertRaises(ValueError):conformal_upper([-1])

    def test_nucleus_birth_and_death_affect_graph(self):
        gt=np.zeros((12,12),int);gt[1:3,1:3]=1;gt[1:3,5:7]=2;gt[8:10,8:10]=3
        pred=gt.copy();pred[pred==2]=0;pred[pred==3]=2
        r=graph_errors(pred,gt,radius=5)
        self.assertEqual(r['matched_nuclei'],2);self.assertEqual(r['added_edges'],1)
        self.assertEqual(r['deleted_edges'],0)
        self.assertEqual(graph_errors(gt,gt)['score'],0.)

    def test_half_iou_splits_are_one_to_one(self):
        gt=np.ones((2,4),int);pred=gt.copy();pred[:,2:]=2
        mapping=match_instances(pred,gt)
        self.assertEqual(len(mapping),1);self.assertEqual(len(set(mapping.values())),1)

    def test_author_decode_first_channel(self):
        mask=np.zeros((6,6,3),np.uint8);mask[1:3,1:3,0]=255;mask[4,4,0]=255
        mask[:,:,1]=255
        self.assertEqual(decode_mask(mask).max(),2)

if __name__=='__main__':unittest.main()

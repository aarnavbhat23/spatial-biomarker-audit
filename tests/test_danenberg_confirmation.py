import unittest
import numpy as np
import pandas as pd
from danenberg_confirmation import retain,choose_images,selected_labels,knn_edges


class IndependentCohortTests(unittest.TestCase):
    def test_annotation_mapping_and_exclusions_are_explicit(self):
        frame=pd.DataFrame({'ImageNumber':[1]*4,'ObjectNumber':range(4),'metabric_id':['MB-01']*4,
            'cellPhenotype':['CD8^{+} T cells','CD4^{+} T cells & APCs','Basal','CD4^{+} T cells'],
            'is_epithelial':[0,0,1,0],'is_tumour':[1]*4,'is_hotAggregate':[0,0,0,1],
            'Location_Center_X':[1,2,3,4],'Location_Center_Y':[1,2,3,4]})
        clean=retain(frame);self.assertEqual(clean['ObjectNumber'].tolist(),[0,2])
        np.testing.assert_equal(selected_labels(clean),[[0,1],[1,0]])

    def test_patient_selection_is_lexical_then_lowest_eligible_image(self):
        counts=pd.DataFrame({'metabric_id':['MB-02','MB-01','MB-01','MB-00'],
                             'ImageNumber':[1,3,2,1],'cells':[900,700,800,499]})
        selected=choose_images(counts)
        self.assertEqual([(x['metabric_id'],x['ImageNumber']) for x in selected],[('MB-01',2),('MB-02',1)])

    def test_graph_ties_are_resolved_by_cell_index(self):
        xy=np.array([[0,0],[1,0],[0,1],[-1,0],[0,-1]])
        edges=knn_edges(xy,2);np.testing.assert_equal(edges[:2],[[0,1],[0,2]])
        self.assertEqual(len(edges),len(xy)*2)
        self.assertTrue(np.all(edges[:,0]!=edges[:,1]))
        with self.assertRaises(ValueError):knn_edges(np.array([[0,0],[0,0],[1,1]]),1)


if __name__=='__main__':unittest.main()

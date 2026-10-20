import unittest
from urllib.error import HTTPError
import numpy as np
from acquisition_emulator import ObservationService,HttpClient,http_service,replay


class EmulatorTests(unittest.TestCase):
    def test_http_observation_boundary_cache_and_review_order(self):
        service=ObservationService(np.array([[1.,2.,3.,4.]]),np.array([10.]),np.array([0,0,1,1]),
                                   np.array([0,1,1,0],bool),np.array([1,1,1,0],bool),1)
        with http_service(service) as url:
            client=HttpClient(url);meta=client.metadata()
            self.assertNotIn('reference',meta);self.assertNotIn('recorded',meta)
            with self.assertRaises(HTTPError) as rejected:client.observe('review',0)
            rejected.exception.close()
            first=client.observe('field',0);self.assertEqual(first['ids'],[0,1]);self.assertEqual(first['calls'],[False,True])
            cached=client.observe('field',0);self.assertEqual(cached['simulated_increment'],0)
            review=client.observe('review',0);self.assertEqual(review['calls'],[True])
            self.assertAlmostEqual(review['simulated_end'],1.15)
            result=replay(client,'joint',width=.01)
            truth=np.array([.6])
            for checkpoint in result['checkpoints']:
                bounds=np.array(checkpoint['intervals']);self.assertTrue(np.all((bounds[:,0]<=truth)&(truth<=bounds[:,1])))


if __name__=='__main__':unittest.main()

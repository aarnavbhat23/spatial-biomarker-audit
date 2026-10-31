import unittest
import numpy as np
from study_measured_barcodes import assign, sample_vectors


class MeasuredBarcodeTests(unittest.TestCase):
    def test_scaling_identity_and_zero_rejection(self):
        codes = np.eye(4)
        identities, accepted, distances, _ = assign(np.array([[3., 0, 0, 0], [0, 5., 0, 0], [0., 0, 0, 0]]), codes)
        self.assertEqual(identities[:2].tolist(), [0, 1])
        self.assertEqual(accepted.tolist(), [True, True, False])
        np.testing.assert_allclose(distances[:2], 0)

    def test_registration_sampling_direction(self):
        y, x = np.indices((30, 30)); signal = y+2*x
        values = sample_vectors(np.array([signal], float), np.array([[10., 12.]]), [[2., -3.]], [1])
        self.assertAlmostEqual(values[0, 0], 8+2*15)

    def test_ambiguous_vector_not_accepted_as_exact_code(self):
        _, accepted, _, _ = assign(np.array([[1., 1., 0, 0]]), np.eye(4), threshold=.5)
        self.assertFalse(accepted[0])


if __name__ == '__main__':
    unittest.main()

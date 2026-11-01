import unittest
import numpy as np
from study_photon_ceiling import calibrated_information, joint_information, information_ceiling


class PhotonCeilingTests(unittest.TestCase):
    def test_shared_channel_offset_annihilates_constant_signal(self):
        ar, ac = np.eye(2), np.ones((2, 1))
        value, weight = calibrated_information(np.ones(2), ar, ac, .1, .2)
        self.assertAlmostEqual(value, 0.)
        np.testing.assert_allclose(ac.T @ weight, 0, atol=1e-12)

    def test_two_round_difference_matches_analytic_ceiling(self):
        ar, ac = np.eye(2), np.ones((2, 1))
        h = np.array([-.5, .5])
        for vg in [1., .1, 1e-8]:
            info, _ = calibrated_information(h, ar, ac, vg, .25)
            self.assertAlmostEqual(info, .5/(vg+.25))
            self.assertAlmostEqual(info, joint_information(h, ar, ac, vg, .25), places=7)
        self.assertAlmostEqual(information_ceiling(h, ar, ac, .25), 2.)

    def test_uncalibrated_new_probe_offset_can_remove_all_added_information(self):
        # New signal dimension and new unrestricted nuisance dimension coincide.
        h = np.array([0., 1.]); nuisance = np.array([[0.], [1.]])
        residual = h-nuisance @ np.linalg.lstsq(nuisance, h, rcond=None)[0]
        self.assertEqual(residual @ residual, 0.)


if __name__ == '__main__':
    unittest.main()

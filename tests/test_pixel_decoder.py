import unittest
import numpy as np
from study_photon_ceiling import cases
from study_pixel_decoder import inputs, render, source_check


class PixelDecoderTests(unittest.TestCase):
    def test_complete_pixel_alias_and_known_registration_separation(self):
        names, codes, dictionary, mapping = inputs()
        self.assertEqual(dictionary.shape, (140, 8, 2))
        for case in cases():
            first, second = [case[k].split(' + ') for k in ['first', 'second']]
            coef = case['lp_coefficients']
            shift = np.array([coef[r]+coef[8+c] for r, c in mapping])
            for distance in [1., 3.]:
                one = render(first, names, codes, mapping, distance, -distance*shift/2)
                two = render(second, names, codes, mapping, distance, distance*shift/2)
                np.testing.assert_allclose(one, two, atol=1e-10)
                registered_one = render(first, names, codes, mapping, distance)
                registered_two = render(second, names, codes, mapping, distance)
                self.assertGreater(float(abs(registered_one-registered_two).max()), 1.)

    def test_probe_contains_information_absent_from_alias(self):
        names, codes, dictionary, mapping = inputs()
        one = render(['Adamdec1', 'Adgrg6'], names, codes, mapping, 1., extra_probe=True)
        two = render(['Ccr2', 'Pik3r1'], names, codes, mapping, 1., extra_probe=True)
        self.assertEqual(one.shape, (9, 2, 25, 25))
        self.assertAlmostEqual(float((one[8]-two[8]).sum()), 200.)
        self.assertEqual(source_check()['commit'], '264c9383a2d73d5669561bff179fd4c359d2e103')

# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,

#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
"""
Helper functions only for DNS simulation calculations
"""
import unittest
from unittest import mock

from mantid.geometry import OrientedLattice
from mantid.kernel import V3D

import numpy as np

from mantidqtinterfaces.dns_simulation.simulation_helpers import (
    two_theta_to_q, two_theta_to_d, angle_to_hkl, check_inplane,
    det_rot_number_to_two_theta,
    hkl_to_omega, ki_from_wavelength, list_to_v3d, lorentz_correction,
    max_q, omega_to_sample_rot, q_to_two_theta, rotate_ub, get_two_theta_path,
    get_omega_path, return_dns_surface_shape, get_angle_q_vs_hkl, return_qxqy,
    return_qx_qx_inplane_refl, shift_channels_below_23, two_theta_to_rot_number,
    get_omega_range_single_crystal, get_two_theta_end, get_two_theta_start,
    get_1deg_angle_linspace,
    get_two_theta_range_sc, get_fwhm, get_peak_width, get_peak, get_intensity_prof,
    get_two_theta_range, get_two_theta_bins, get_unique_refl, get_inplane_refl,
    get_unique_inplane_refl)

UB = np.asarray([[0.00000000e+00, -1.74272370e-17, 4.03582196e-02],
                 [0.00000000e+00, -2.84608379e-01, 0.00000000e+00],
                 [3.28637448e-01, 1.64318724e-01, 2.47122822e-18]])
UBINV = np.linalg.inv(UB)
TTH = 88.11
WAVELENGTH = 4.2
OEMEGA = -47.05 + 8.11  # det -8.11 chan 16 # sr -47.05
Q1 = np.asarray([1, 0, 0])
Q2 = np.asarray([1, 0, 1])
ORILAT = OrientedLattice(3.5, 3.5, 24, 90, 90, 120)


class DNSReductionGUIPresenterTest(unittest.TestCase):
    # pylint: disable=protected-access, too-many-public-methods

    def test_two_theta_to_q(self):
        testv = two_theta_to_q(TTH, WAVELENGTH)
        self.assertAlmostEqual(testv, 0.21733483834214878)

    def test_two_theta_to_d(self):
        testv = two_theta_to_d(TTH, WAVELENGTH)
        self.assertAlmostEqual(testv, 28.9101616432429)

    def test_angle_to_hkl(self):
        testv = angle_to_hkl(TTH, OEMEGA, UBINV, WAVELENGTH)
        self.assertAlmostEqual(testv[0], 1.00002826)
        self.assertAlmostEqual(testv[1], 0)
        self.assertAlmostEqual(testv[2], -1.00058571)

    def test_check_inplane(self):
        hkl = np.asarray([1, 0, 1])
        testv = check_inplane(Q1, Q2, hkl)
        self.assertTrue(testv)
        hkl = np.asarray([1, 1, 1])
        testv = check_inplane(Q1, Q2, hkl)
        self.assertFalse(testv)

    def test_det_rot_number_to_two_theta(self):
        testv = det_rot_number_to_two_theta(-20, 3)
        self.assertEqual(testv, 35)

    def test_hkl_to_omega(self):
        hkl = np.asarray([1, 0, 1])
        testv = hkl_to_omega(hkl, UB, 4.2, TTH)
        self.assertAlmostEqual(testv, -38.94412102205903)

    def test_ki_from_wavelength(self):
        testv = ki_from_wavelength(4.2)
        self.assertAlmostEqual(testv, 1.4959965017094252)

    def test_list_to_v3d(self):
        testv = list_to_v3d([1, 2, 3])
        self.assertIsInstance(testv, V3D)
        self.assertEqual(testv[0], 1.0)
        self.assertEqual(testv[1], 2.0)
        self.assertEqual(testv[2], 3.0)

    def test_lorentz_correction(self):
        testv = lorentz_correction(10000, TTH)
        self.assertAlmostEqual(testv, 10005.443087213083)

    def test_max_q(self):
        testv = max_q(WAVELENGTH)
        self.assertAlmostEqual(testv, 2.769210712401315)

    def test_omega_to_sample_rot(self):
        testv = omega_to_sample_rot(100, -5)
        self.assertEqual(testv, 95)

    def test_q_to_two_theta(self):
        testv = q_to_two_theta(2, 4)
        self.assertAlmostEqual(testv, 79.08044749562039)

    def test_rotate_ub(self):
        testv = rotate_ub(36, UB)
        rot_ub = np.asarray([[-0.19316825, -0.09658412, 0.03265049],
                             [0., -0.28460838, 0.],
                             [0.26587328, 0.13293664, 0.02372197]])
        self.assertTrue((rot_ub == testv).all)

    def test_get_two_theta_path(self):
        testv = get_two_theta_path(range(0, 10), range(2, 5))
        dpath = np.asarray([0., 0., 0., 0., 1., 2., 3., 4.,
                            5., 6., 7., 8., 9., 9., 9., 9.,
                            9., 8., 7., 6., 5., 4., 3., 2.,
                            1., 0.])
        self.assertTrue((testv == dpath).all())

    def test_get_omega_path(self):
        testv = get_omega_path(range(0, 10), range(2, 5))
        dpath = np.asarray([4., 3., 2., 2., 2., 2., 2., 2., 2., 2.,
                            2., 2., 2., 2., 3., 4., 4., 4., 4., 4.,
                            4., 4., 4., 4., 4., 4.])
        self.assertTrue((testv == dpath).all())

    def test_return_dns_surface_shape(self):
        options = {'sc_det_start': -5,
                   'sc_det_end': -6,
                   'sc_sam_start': 4,
                   'sc_sam_end': 5}
        testv = return_dns_surface_shape(ORILAT, [1, 1, 0], [0, 0, 1],
                                         WAVELENGTH, options)
        self.assertEqual(testv.shape, (238, 2))
        self.assertAlmostEqual(testv[-1, 0], -0.031210269114692483)
        self.assertAlmostEqual(testv[0, 1], -0.06506825759735697)
        self.assertAlmostEqual(testv[111, 1], 7.1228993468437665)
        self.assertAlmostEqual(testv[111, 0], -0.4049968852041988)

    def test_get_angle_q_vs_hkl(self):
        testv = get_angle_q_vs_hkl([1, 1, 0], [1, 1, 3], ORILAT)
        self.assertEqual(testv, 0.21535769969773777)

    def test_return_qxqy(self):
        testv = return_qxqy(ORILAT, [1, 1, 0], [0, 0, 1], [1, 1, 2])
        self.assertAlmostEqual(testv[0], 1)
        self.assertAlmostEqual(testv[1], 2)

    def test_return_qx_qx_inplane_refl(self):
        refl1 = mock.Mock()
        refl1.hkl = [1, 1, 2]
        refl1.h, refl1.k, refl1.l = refl1.hkl
        refl1.inplane = True
        refl1.fs_lc = 10
        refl2 = mock.Mock()
        refl2.inplane = False
        refl2.hkl = [1, 0, 1]
        refl2.h, refl2.k, refl2.l = refl2.hkl
        refl2.fs_lc = 20
        testv = return_qx_qx_inplane_refl(ORILAT, [1, 1, 0], [0, 0, 1],
                                          [refl1, refl2])
        self.assertIsInstance(testv, np.ndarray)
        self.assertEqual(testv.shape, (1, 6))
        control = np.asarray([[1., 2., 10., 1., 1., 2.]])
        self.assertTrue(np.allclose(testv, control))

    def test_shift_channels_below_23(self):
        testv = shift_channels_below_23(26, -5)
        self.assertEqual(testv, [23, -20])

    def test_two_theta_to_rot_number(self):
        testv = two_theta_to_rot_number(51)
        self.assertEqual(testv, [-6, 9])

    def test_get_omega_range_single_crystal(self):
        testv = get_omega_range_single_crystal(10, 12, -5)
        self.assertEqual(len(testv), 3)
        self.assertAlmostEqual(testv[0], 15)
        self.assertAlmostEqual(testv[1], 16)
        self.assertAlmostEqual(testv[2], 17)

    def test_get_two_theta_end(self):
        testv = get_two_theta_end(-10, 5)
        self.assertEqual(testv, 130)

    def test_get_two_theta_start(self):
        testv = get_two_theta_start(-10, 5)
        self.assertEqual(testv, 15)

    def test_get_1deg_angle_linspace(self):
        testv = get_1deg_angle_linspace(10, 14)
        self.assertTrue((testv == np.asarray([10, 11, 12, 13, 14])).all())

    def test_get_two_theta_range_sc(self):
        testv = get_two_theta_range_sc(2, 3)
        self.assertEqual(len(testv), 115)
        self.assertEqual(testv[0], -2)
        self.assertEqual(testv[1], -1)
        self.assertEqual(testv[114], 112)

    def test_get_fwhm(self):
        testv = get_fwhm(100)
        self.assertAlmostEqual(testv, 0.45887253942604783)

    def test_get_peak_width(self):
        testv = get_peak_width(0.45)
        self.assertAlmostEqual(testv, 0.1910974050648043)

    def test_get_peak(self):
        testv = get_peak(0.45, np.asarray([10, 11]), 11)
        self.assertAlmostEqual(testv[0], 0.07505256)
        self.assertAlmostEqual(testv[1], 0.8865384)

    def test_get_intensity_prof(self):
        refl = mock.Mock()
        refl.tth = 100
        refl.fs_lc = 10
        testv = get_intensity_prof(refl, 0, np.asarray([100, 101]))
        self.assertAlmostEqual(testv[0], 2.04727282e+01)
        self.assertAlmostEqual(testv[1], 3.91407935e-05)

    def test_get_two_theta_range(self):
        testv = get_two_theta_range(1, 3, 2)
        self.assertEqual(len(testv), 1131)
        self.assertEqual(testv[0], 1)
        self.assertEqual(testv[1], 1.1)
        self.assertEqual(testv[-1], 114)

    def test_get_two_theta_bins(self):
        testv = get_two_theta_bins(np.asarray([0, 1]), 0.1)
        self.assertAlmostEqual(testv[0], -0.05)
        self.assertAlmostEqual(testv[1], 0.95)
        self.assertAlmostEqual(testv[2], 1.05)

    def test_get_unique_refl(self):
        refl1 = mock.Mock()
        refl1.unique = True
        refl2 = mock.Mock()
        refl2.unique = False
        testv = get_unique_refl([refl1, refl2])
        self.assertEqual(testv[0], refl1)
        self.assertEqual(len(testv), 1)

    def test_get_inplane_refl(self):
        refl1 = mock.Mock()
        refl1.inplane = True
        refl2 = mock.Mock()
        refl2.inplane = False
        testv = get_inplane_refl([refl1, refl2])
        self.assertEqual(testv[0], refl1)
        self.assertEqual(len(testv), 1)

    def test_get_unique_inplane_refl(self):
        refl1 = mock.Mock()
        refl1.inplane = True
        refl1.hkl = [1, 1, 0]
        refl1.equivalents = [[1, 0, 0]]
        refl2 = mock.Mock()
        refl2.inplane = False
        refl3 = mock.Mock()
        refl3.inplane = True
        refl3.hkl = [1, 0, 0]
        refl3.equivalents = [[1, 1, 0]]
        testv = get_unique_inplane_refl([refl1, refl2, refl3])
        self.assertEqual(len(testv), 1)
        self.assertEqual(testv[0], refl1)


if __name__ == '__main__':
    unittest.main()

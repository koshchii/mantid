# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2019 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
"""
DNS test for sc elastic scripts
"""

import unittest
from unittest.mock import patch
from unittest.mock import call

from mantidqtinterfaces.dns_single_crystal_elastic.scripts.md_single_crystal_elastic import load_all, load_binned, vanadium_correction
from mantidqtinterfaces.dns_powder_tof.helpers.helpers_for_testing import (
    get_fake_MD_workspace_unique,
    get_fake_elastic_data_dic,
    get_fake_elastic_single_crystal_options,
    get_fake_elastic_single_crystal_binning,
)
from mantidqtinterfaces.dns_powder_tof.data_structures.dns_error import DNSError
from mantid.simpleapi import mtd
from mantid.api import IMDHistoWorkspace


class MDPowderElasticTest(unittest.TestCase):
    def setUp(self):
        mtd.clear()

    @staticmethod
    @patch("mantidqtinterfaces.dns_single_crystal_elastic.scripts.md_single_crystal_elastic.load_binned")
    def test_load_all(mock_load_binned):
        data_dict = get_fake_elastic_data_dic()
        load_all(data_dict, [0, 1, 2], {"a": 1})
        calls = [
            call("knso_x_nsf", [0, 1, 2], {"a": 1}, "C:/data", range(554574, 554634, 6), False),
            call("knso_x_sf", [0, 1, 2], {"a": 1}, "C:/data", range(554573, 554633, 6), False),
            call("knso_y_nsf", [0, 1, 2], {"a": 1}, "C:/data", range(554576, 554636, 6), False),
            call("knso_y_sf", [0, 1, 2], {"a": 1}, "C:/data", range(554575, 554635, 6), False),
            call("knso_z_nsf", [0, 1, 2], {"a": 1}, "C:/data", range(554578, 554638, 6), False),
            call("knso_z_sf", [0, 1, 2], {"a": 1}, "C:/data", range(554577, 554637, 6), False),
        ]
        mock_load_binned.assert_has_calls(calls)

    @patch("mantidqtinterfaces.dns_single_crystal_elastic.scripts.md_single_crystal_elastic.mtd")
    @patch("mantidqtinterfaces.dns_single_crystal_elastic.scripts.md_single_crystal_elastic.BinMD")
    @patch("mantidqtinterfaces.dns_single_crystal_elastic.scripts.md_single_crystal_elastic.LoadDNSSCD")
    def test_load_binned(self, mock_load, mock_binmd, mock_mtd):
        params = get_fake_elastic_single_crystal_options()
        binning = get_fake_elastic_single_crystal_binning()
        testv = load_binned("knso_x_nsf", binning, params, "C:/data", range(554574, 554578, 2), False)
        mock_load.assert_called_once_with(
            FileNames="C:/data_554574.d_dat, C:/data_554576.d_dat",
            OutputWorkspace="knso_x_nsf",
            NormalizationWorkspace="knso_x_nsf_norm",
            Normalization="monitor",
            a=2,
            b=3,
            c=4,
            alpha=78,
            beta=86,
            gamma=85,
            OmegaOffset=0,
            HKL1="1,2,3",
            HKL2="2,3,4",
            LoadAs="raw",
            SaveHuberTo="huber_x_nsf",
        )

        calls = [
            call(
                InputWorkspace="knso_x_nsf",
                OutputWorkspace="knso_x_nsf",
                AxisAligned=True,
                AlignedDim0="Theta,0.0,0.5,2",
                AlignedDim1="Omega,3,4,5",
            ),
            call(
                InputWorkspace="knso_x_nsf_norm",
                OutputWorkspace="knso_x_nsf_norm",
                AxisAligned=True,
                AlignedDim0="Theta,0.0,0.5,2",
                AlignedDim1="Omega,3,4,5",
            ),
        ]
        mock_binmd.assert_has_calls(calls)
        mock_mtd.__getitem__.assert_called_once_with("knso_x_nsf")
        self.assertEqual(testv, mock_mtd.__getitem__.return_value)
        mock_load.reset_mock()
        load_binned("knso_x_nsf", binning, params, "C:/data", range(554574, 554578, 2), True)
        mock_load.assert_called_once_with(
            FileNames="C:/data_554574.d_dat, C:/data_554576.d_dat",
            OutputWorkspace="knso_x_nsf",
            NormalizationWorkspace="knso_x_nsf_norm",
            Normalization="monitor",
            a=2,
            b=3,
            c=4,
            alpha=78,
            beta=86,
            gamma=85,
            OmegaOffset=0,
            HKL1="1,2,3",
            HKL2="2,3,4",
            LoadAs="raw",
            LoadHuberFrom="huber_x_nsf",
        )

    def test_vanadium_correction_1(self):
        get_fake_MD_workspace_unique(name="sample_x_sf", factor=5.2)
        get_fake_MD_workspace_unique(name="sample_x_sf_norm", factor=0.7)
        with self.assertRaises(DNSError) as context:
            vanadium_correction("sample_x_sf", vana_set=None, ignore_vana_fields=False, sum_vana_sf_nsf=False)
        self.assertTrue("No vanadium file for" in str(context.exception))
        get_fake_MD_workspace_unique(name="vana_x_sf", factor=1.3)
        get_fake_MD_workspace_unique(name="vana_x_sf_norm", factor=1.1)
        testv = vanadium_correction("sample_x_sf", vana_set=None, ignore_vana_fields=False, sum_vana_sf_nsf=False)
        self.assertIsInstance(testv, IMDHistoWorkspace)
        testv = mtd["sample_x_sf"]
        value = testv.getSignalArray()
        error = testv.getErrorSquaredArray()
        self.assertEqual(value.shape, (5, 5, 1))
        self.assertAlmostEqual(value[0, 0, 0], 6.040404040404041)
        self.assertAlmostEqual(error[0, 0, 0], 1.7832994702647458)
        testv = mtd["sample_x_sf_norm"]
        value = testv.getSignalArray()
        error = testv.getErrorSquaredArray()
        self.assertEqual(value.shape, (5, 5, 1))
        self.assertAlmostEqual(value[0, 0, 0], 53.30769230769231)
        self.assertAlmostEqual(error[0, 0, 0], 130.06516158397812)


if __name__ == "__main__":
    unittest.main()

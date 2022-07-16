# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
import unittest
from unittest import mock
from unittest.mock import patch
import numpy as np

from mantidqtinterfaces.dns_simulation.simulation_sub_sc_model import \
    DNSSimulationSubScModel


class DNSScriptGeneratorModelTest(unittest.TestCase):
    parent = None

    @classmethod
    def setUpClass(cls):
        cls.parent = mock.Mock()
        cls.parent.update_progress = mock.Mock()
        cls.model = DNSSimulationSubScModel(parent=cls.parent)

    def test___init__(self):
        self.assertIsInstance(self.model, DNSSimulationSubScModel)

    @patch('mantidqtinterfaces.dns_simulation.simulation_sub_sc_model.'
           'sim_help')
    def test_create_dns_surface(self, mock_sim_help):
        testv = self.model.create_dns_surface(1, 2, 3, 4, 5)
        mock_sim_help.return_dns_surface_shape.assert_called_once_with(1, 2, 3,
                                                                       4, 5)
        self.assertEqual(
            testv, mock_sim_help.return_dns_surface_shape.return_value)

    def test_get_min_max_int(self):
        refls = np.asarray([[0, 1, 5], [3, 4, 2]])
        testv = self.model.get_min_max_int(refls)
        self.assertEqual(testv, [2, 5])
        refls = np.asarray([])
        testv = self.model.get_min_max_int(refls)
        self.assertEqual(testv, [0, 1])

    def test_get_hkl_on_plot(self):
        testv = self.model.get_hkl_on_plot(1, 2, [1, 2, 3], [4, 5, 6])
        self.assertTrue((testv == np.asarray([9, 12, 15])).all())
        testv = self.model.get_hkl_on_plot(None, 2, [1, 2, 3], [4, 5, 6])
        self.assertIsNone(testv)
        testv = self.model.get_hkl_on_plot(1, None, [1, 2, 3], [4, 5, 6])
        self.assertIsNone(testv)
        testv = self.model.get_hkl_on_plot(None, None, [1, 2, 3], [4, 5, 6])
        self.assertIsNone(testv)


if __name__ == '__main__':
    unittest.main()

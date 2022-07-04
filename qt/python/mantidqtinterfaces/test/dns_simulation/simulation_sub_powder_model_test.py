# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
import unittest
from unittest import mock
from unittest.mock import patch

from mantidqtinterfaces.dns_simulation.simulation_sub_powder_model import \
    DNSSimulationSubPowderModel


class DNSSimulationSubPowderModelTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = DNSSimulationSubPowderModel()

    def test__init(self):
        self.assertIsInstance(self.model, DNSSimulationSubPowderModel)

    @patch('mantidqtinterfaces.dns_simulation.simulation_sub_powder_model.'
           'sim_help')
    @patch('mantidqtinterfaces.dns_simulation.simulation_sub_powder_model.'
           'CreateWorkspace')
    def test_create_powder_profile(self, mock_create_ws, mock_sim_help):
        mock_sim_help.get_tth_range.return_value = 3
        mock_sim_help.get_tth_bins.return_value = []
        mock_sim_help.get_intensity_prof.return_value = 10
        testv = self.model.create_powder_profile([1], 0, 5, 2)
        mock_sim_help.get_tth_range.assert_called_once_with(0, 5, 2)
        mock_sim_help.get_intensity_prof.assert_called_once_with(1, 2, 3)
        mock_sim_help.get_tth_bins.assert_called_once_with(3)
        mock_create_ws.assert_called_once_with(
            OutputWorkspace='mat_simulation',
            DataX=[],
            DataY=10,
            NSpec=1,
            UnitX='Degrees')
        self.assertEqual(testv, [3, 10])

    @patch('mantidqtinterfaces.dns_simulation.simulation_sub_powder_model.'
           'sim_help')
    def test_get_annotation_list(self, mock_sim_help):
        refl = mock.Mock()
        refl.tth = 3
        refl.hkl = [1, 2, 3]
        mock_sim_help.get_tth_end.return_value = 20
        mock_sim_help.get_tth_start.return_value = 0
        mock_sim_help.get_unique_refl.return_value = [refl]

        testv = self.model.get_annotation_list([refl], 0, 20, 2, [1] * 51)
        mock_sim_help.get_tth_end.assert_called_once_with(20, 2)
        mock_sim_help.get_tth_start.assert_called_once_with(0, 2)
        mock_sim_help.get_unique_refl.assert_called_once_with([refl])
        self.assertEqual(testv, [[5], [[1, 2, 3]], [1]])


if __name__ == '__main__':
    unittest.main()

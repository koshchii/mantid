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

from mantidqtinterfaces.dns_simulation.simulation_sub_common_presenter import \
    DNSSimulationSubCommonPresenter
from mantidqtinterfaces.dns_simulation.simulation_sub_sc_presenter import \
    DNSSimulationSubScPresenter
from mantidqtinterfaces.dns_simulation.simulation_sub_sc_model import \
    DNSSimulationSubScModel
from mantidqtinterfaces.dns_simulation.simulation_sub_sc_view import \
    DNSSimulationSubScView


class DNSSimulationSubScPresenterTest(unittest.TestCase):
    # pylint: disable=protected-access

    parent = None
    view = None
    model = None

    @classmethod
    def setUpClass(cls):
        cls.parent = mock.Mock()
        cls.parent.update_progress = mock.Mock()
        cls.view = mock.create_autospec(DNSSimulationSubScView)
        cls.view.sig_scplot_clicked.connect = mock.Mock()
        cls.view.sig_mouse_pos_changed.connect = mock.Mock()
        cls.model = mock.create_autospec(DNSSimulationSubScModel)
        cls.presenter = DNSSimulationSubScPresenter(parent=cls.parent,
                                                    view=cls.view,
                                                    model=cls.model,
                                                    )

    def setUp(self):
        self.view.reset_mock()
        self.model.reset_mock()

    def test___init__(self):
        self.assertIsInstance(self.presenter, DNSSimulationSubScPresenter)
        self.assertIsInstance(self.presenter, DNSSimulationSubCommonPresenter)
        self.assertTrue(hasattr(self.presenter, '_sub_dict'))

    @patch('mantidqtinterfaces.dns_simulation.simulation_sub_sc_presenter.'
           'DNSSimulationSubScPresenter._toggle_off_set_warning')
    @patch('mantidqtinterfaces.dns_simulation.simulation_sub_sc_presenter.'
           'DNSSimulationSubScPresenter._sc_plot')
    def test_process_request(self, mock_sc_plot, mock_oo_warning):
        self.presenter.process_request({})
        self.assertEqual(self.presenter._sub_dict, {})
        mock_sc_plot.assert_called_once()
        mock_oo_warning.assert_called_once_with(False)
        self.presenter.process_request({'oof_set': True})
        mock_oo_warning.assert_called_with(True)

    @patch('mantidqtinterfaces.dns_simulation.simulation_sub_sc_presenter.'
           'DNSSimulationSubScPresenter._scatter_plot')
    def test__sc_plot(self, mock_scatter):
        self.view.get_state.return_value = {}
        self.presenter._sub_dict = None
        self.presenter._sc_plot()
        self.model.create_dns_surface.assert_not_called()
        self.presenter._sub_dict = {'hkl1_v': [1, 0, 0],
                                    'hkl2_p_v': [0, 1, 0],
                                    'inplane_refls': [1],
                                    'wavelength': 4,
                                    'oriented_lattice': 3}
        self.presenter._sc_plot()
        self.model.create_dns_surface.assert_called_once_with(
            3, [1, 0, 0], [0, 1, 0], 4, {})
        self.view.start_sc_plot.assert_called_once_with(
            self.model.create_dns_surface.return_value)
        mock_scatter.assert_called_once_with([1])
        self.view.finish_sc_plot.assert_called_once_with(
            [1, 0, 0], '[0.000, 1.000, 0.000]')

    def test__scatter_plot(self):
        self.model.get_min_max_int.return_value = [1, 2]
        refls = np.asarray([])
        self.presenter._scatter_plot(refls)
        self.view.scatter_plot.assert_not_called()
        refls = np.asarray([[1, 2, 3, 4, 5, 6]])
        self.presenter._scatter_plot(refls)
        self.model.get_min_max_int.assert_called_with(refls)
        self.view.scatter_plot.assert_called_once_with(
            x=1, y=2, intensity=3, intensity_max_min=[1, 2])

    def test_annotate_refl(self):
        refls = np.asarray([[1, 2, 3, 4, 5, 6]])
        self.presenter._annotate_refl(refls)
        self.view.annotate_refl.assert_called_once_with(
            '    3\n    [ 4.00,  5.00,  6.00]', 1, 2)

    def test__set_hkl_pos_on_plot(self):
        self.presenter._sub_dict = None
        self.presenter._set_hkl_pos_on_plot(1, 2)
        self.model.get_hkl_on_plot.assert_not_called()
        self.presenter._sub_dict = {'hkl1_v': 3,
                                    'hkl2_p_v': 4}
        self.presenter._set_hkl_pos_on_plot(1, 2)
        self.model.get_hkl_on_plot.assert_called_once_with(1, 2, 3, 4)
        self.view.set_hkl_position_on_plot.assert_called_once_with(
            self.model.get_hkl_on_plot.return_value)

    def test__toggle_off_set_warning(self):
        self.presenter._toggle_off_set_warning(False)
        self.view.set_omega_offset_warning.assert_called_once_with(
            'Warning: omega offset not set')
        self.view.reset_mock()
        self.presenter._toggle_off_set_warning(True)
        self.view.set_omega_offset_warning.assert_called_once_with('')


if __name__ == '__main__':
    unittest.main()

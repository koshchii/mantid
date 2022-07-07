# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
"""
Presenter for DNS simulation
"""
import unittest
from unittest import mock
from unittest.mock import patch
from mantidqtinterfaces.dns_powder_tof.data_structures.dns_observer import DNSObserver
from mantidqtinterfaces.dns_simulation.simulation_presenter import DNSSimulationPresenter

from mantidqtinterfaces.dns_simulation.simulation_model import DNSSimulationModel
from mantidqtinterfaces.dns_simulation.simulation_view import DNSSimulationView


class DNSSimulationPresenterTest(unittest.TestCase):
    # pylint: disable=protected-access, too-many-public-methods, too-many-arguments

    view = None
    parent = None
    model = None
    subwidget = None

    @classmethod
    def setUpClass(cls):
        cls.subwidget = mock.Mock()
        cls.parent = mock.Mock()
        cls.view = mock.create_autospec(DNSSimulationView)
        cls.model = mock.create_autospec(DNSSimulationModel)
        cls.view.sig_cif_set.connect = mock.Mock()
        cls.view.sig_unit_cell_changed.connect = mock.Mock()
        cls.view.sig_wavelength_changed.connect = mock.Mock()
        cls.view.sig_calculate_clicked.connect = mock.Mock()
        cls.view.sig_sample_rot_changed.connect = mock.Mock()
        cls.view.sig_inplane_unique_switched.connect = mock.Mock()
        cls.view.get_state.return_value = {'hkl1': [1, 2, 3],
                                           'hkl2': [4, 5, 6],
                                           'fix_omega': True,
                                           'det_rot': -5,
                                           'sample_rot': 20,
                                           'wavelength': 4,
                                           'hkl1_v': [1, 0, 0],
                                           'hkl2_v': [1, 1, 0],
                                           'hkl2_p_v': [0, 1, 0]}

        cls.presenter = DNSSimulationPresenter(view=cls.view,
                                               model=cls.model,
                                               parent=cls.parent,
                                               subwidgets=[cls.subwidget])

    def setUp(self):
        self.view.reset_mock()
        self.model.reset_mock()
        self.parent.reset_mock()

    def test__init__(self):
        self.assertIsInstance(self.presenter, DNSSimulationPresenter)
        self.assertIsInstance(self.presenter, DNSObserver)
        self.assertTrue(hasattr(self.presenter, '_refls'))
        self.assertTrue(hasattr(self.presenter, 'sub_presenters'))
        self.assertEqual(len(self.presenter.sub_presenters), 1)
        self.assertTrue(hasattr(self.presenter, '_refls'))

    def test__get_sub_dict(self):
        self.presenter._refls = None
        self.model.return_reflections_in_map.return_value = [2]
        self.model.get_oriented_lattice.return_value = 3
        self.model.filter_refls.return_value = [1]
        self.parent.presenter.own_dict = {'inplane': True, 'unique': False}
        testv = self.presenter._get_sub_dict()
        self.model.return_reflections_in_map.assert_called_once_with(
            [1, 0, 0], [0, 1, 0], None)
        self.model.get_oriented_lattice.assert_called_once()
        self.model.filter_refls.assert_called_once()
        subdict = {'filtered_refls': [1], 'tth_limit': 5, 'refls': None,
                   'hkl1_v': [1, 0, 0], 'hkl2_p_v': [0, 1, 0],
                   'inplane_refls': [2], 'wavelength': 4, 'oriented_lattice': 3,
                   'oof_set': True}
        self.assertEqual(testv, subdict)

    @patch('mantidqtinterfaces.dns_simulation.simulation_presenter.'
           'DNSSimulationPresenter._get_sub_dict')
    def test_request_to_subwidget(self, mock_subdict):
        self.presenter.request_to_subwidget()
        mock_subdict.assert_called_once()
        self.subwidget.presenter.process_request.assert_called_once_with(
            mock_subdict.return_value)

    def test_back_call_from_tableitem_clicked(self):
        self.view.get_state.return_value['fix_omega'] = True
        self.presenter.back_call_from_table_item_clicked(-6, 25)
        self.model.get_oof_from_ident.assert_not_called()
        self.view.get_state.return_value['fix_omega'] = False
        self.presenter.back_call_from_table_item_clicked(-6, 25)
        self.model.get_oof_from_ident.assert_called_once_with(
            -6, 25, 20, -5)
        self.view.set_omega_offset.assert_called_once_with(
            self.model.get_oof_from_ident.return_value)

    def test_sample_rot_changed(self):
        self.view.get_state.return_value['fix_omega'] = True
        self.presenter._sample_rot_changed()
        self.view.set_omega_offset.assert_not_called()
        self.view.get_state.return_value['fix_omega'] = False
        self.presenter._sample_rot_changed()
        self.view.set_omega_offset.assert_called_once_with(0)

    @patch('mantidqtinterfaces.dns_simulation.simulation_presenter.'
           'DNSSimulationPresenter.raise_error')
    def test_get_and_validate(self, mock_raise):
        self.model.validate_hkl.return_value = [True, '']
        testv = self.presenter._get_and_validate()
        self.model.validate_hkl.assert_called_once_with([1, 0, 0], [1, 1, 0])
        self.assertEqual(testv, True)
        self.model.validate_hkl.return_value = [False, '123']
        self.presenter._get_and_validate()
        mock_raise.assert_called_once_with('123')

    @patch('mantidqtinterfaces.dns_simulation.simulation_presenter.'
           'DNSSimulationPresenter._get_and_validate')
    @patch('mantidqtinterfaces.dns_simulation.simulation_presenter.'
           'DNSSimulationPresenter._perp_inplane')
    @patch('mantidqtinterfaces.dns_simulation.simulation_presenter.'
           'DNSSimulationPresenter._d_tooltip')
    @patch('mantidqtinterfaces.dns_simulation.simulation_presenter.'
           'DNSSimulationPresenter.request_to_subwidget')
    @patch('mantidqtinterfaces.dns_simulation.simulation_presenter.'
           'DNSSimulationPresenter.raise_error')
    def test__calculate(self, mock_raise, mock_reqsu, mock_dtool, mock_inplane,
                        mock_val):
        mock_val.return_value = False
        self.presenter._calculate()
        self.model.get_refls_and_set_orientation.assert_not_called()
        mock_val.return_value = True
        self.model.get_refls_and_set_orientation.side_effect = ValueError()
        self.presenter.own_dict = {}
        self.presenter._calculate()
        self.model.get_refls_and_set_orientation.assert_called_once_with({})
        mock_raise.assert_called_once()
        mock_inplane.assert_not_called()
        self.model.get_refls_and_set_orientation.side_effect = None
        self.presenter._calculate()
        mock_inplane.assert_called_once()
        mock_reqsu.assert_called_once()
        mock_dtool.assert_called_once()

    def test_cif_set(self):
        self.model.load_cif.return_value = {}
        self.presenter._cif_set('123.cif')
        self.model.load_cif.assert_called_once_with('123.cif')
        self.view.set_state.assert_called_once_with({})
        self.assertEqual(self.presenter.own_dict['cif_set'], True)

    def test__d_tooltip(self):
        self.model.get_ds.return_value = 1, 2, 3
        self.presenter._d_tooltip()
        self.model.get_ds.assert_called_once_with([1, 0, 0], [1, 1, 0],
                                                  [0, 1, 0])
        self.view.set_d_tooltip.assert_called_once_with(1, 2, 3)

    def test__perp_inplane(self):
        self.model.get_hkl2_p.return_value = [0, 1, 0]
        self.presenter._perp_inplane()
        self.model.get_hkl2_p.assert_called_once()
        self.view.set_hkl2_p.assert_called_once_with([0, 1, 0])

    def test_set_ki(self):
        self.presenter._set_ki()
        self.model.get_ki.assert_called_once_with(4)
        self.view.set_ki.assert_called_once_with(
            self.model.get_ki.return_value)

    def test_set_spacegroup(self):
        self.presenter._set_spacegroup(114)
        self.view.set_spacegroup.assert_called_once_with(114)

    def test_unitcell_changed(self):
        self.presenter.own_dict['cif_set'] = True
        self.presenter._unit_cell_changed()
        self.assertFalse(self.presenter.own_dict['cif_set'])

    def test_get_option_dict(self):
        testv = self.presenter.get_option_dict()
        self.view.get_state.assert_called_once()
        self.model.get_hkl_vector_dict.assert_called_once_with([1, 2, 3],
                                                               [4, 5, 6])
        self.assertIsInstance(testv, dict)


if __name__ == '__main__':
    unittest.main()

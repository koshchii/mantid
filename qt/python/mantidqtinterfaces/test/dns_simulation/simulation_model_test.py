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

import mantidqtinterfaces.dns_simulation.simulation_helpers as sim_help
from mantidqtinterfaces.dns_powder_tof.data_structures.dns_obs_model import DNSObsModel
from mantidqtinterfaces.dns_simulation.simulation_model import DNSSimulationModel


class DNSSimulationModelTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        parent = mock.Mock()
        cls.model = DNSSimulationModel(parent)

    def test___init__(self):
        self.assertIsInstance(self.model, DNSSimulationModel)
        self.assertIsInstance(self.model, DNSObsModel)
        self.assertTrue(hasattr(self.model, '_sim_ws'))
        self.assertTrue(hasattr(self.model, '_orilat'))
        self.assertTrue(hasattr(self.model, '_refls'))
        self.assertTrue(hasattr(self.model, '_generator'))
        self.assertTrue(hasattr(self.model, '_cryst'))
        self.assertTrue(hasattr(self.model, '_non_rot_lat'))

    @patch('mantidqtinterfaces.dns_simulation.simulation_model.'
           'sim_help.get_unique_inplane_refl')
    @patch('mantidqtinterfaces.dns_simulation.simulation_model.'
           'sim_help.get_unique_refl')
    @patch('mantidqtinterfaces.dns_simulation.simulation_model.'
           'sim_help.get_inplane_refl')
    def test_filter_refls(self, mock_inplane, mock_unique, mock_both):
        mock_inplane.return_value = 1
        mock_unique.return_value = 2
        mock_both.return_value = 3
        testv = self.model.filter_refls([], True, False)
        self.assertEqual(testv, 1)
        testv = self.model.filter_refls([], False, True)
        self.assertEqual(testv, 2)
        testv = self.model.filter_refls([], True, True)
        self.assertEqual(testv, 3)
        testv = self.model.filter_refls([], False, False)
        self.assertEqual(testv, [])

    @patch('mantidqtinterfaces.dns_simulation.simulation_model.'
           'ReflectionGenerator')
    @patch('mantidqtinterfaces.dns_simulation.simulation_model.'
           'sim_help.list_to_v3d')
    def test_get_ds(self, mock_list_to_v3d, mock_gen):
        self.model._generator = mock_gen
        mock_list_to_v3d.return_value = []
        mock_gen.getDValues.return_value = 1, 2, 3
        testv = self.model.get_ds([1, 2, 3], [3, 5, 6], [7, 8, 9])
        self.assertEqual(mock_list_to_v3d.call_count, 3)
        mock_list_to_v3d.assert_any_call([1, 2, 3])
        mock_list_to_v3d.assert_any_call([3, 5, 6])
        mock_list_to_v3d.called_with([7, 8, 9])
        self.assertEqual(testv, (1, 2, 3))
        mock_gen.getDValues.assert_called_once_with([[], [], []])

    @patch('mantidqtinterfaces.dns_simulation.simulation_model.'
           'OrientedLattice')
    def test_get_hkl2_p(self, mock_ori):
        self.model._non_rot_lat = mock_ori
        mock_ori.getvVector.return_value = np.asarray([1, 2, 3])
        testv = self.model.get_hkl2_p()
        self.assertTrue(
            np.allclose(testv, [0.26726124, 0.53452248, 0.80178373]))

    @patch('mantidqtinterfaces.dns_simulation.simulation_model.'
           'sim_help.list_to_v3d')
    @patch('mantidqtinterfaces.dns_simulation.simulation_model.'
           'OrientedLattice')
    def test__set_lattice(self, mock_ori, mock_v3d):
        self.model._orilat = mock_ori
        self.model._non_rot_lat = mock_ori
        mock_v3d.return_value = 1
        self.model._set_lattice([1, 2, 3], [4, 5, 6])
        self.assertEqual(mock_ori.setUFromVectors.call_count, 2)
        self.assertEqual(mock_v3d.call_count, 4)
        mock_ori.setUFromVectors.assert_called_with(1, 1)
        mock_v3d.assert_called_with([4, 5, 6])

    @patch('mantidqtinterfaces.dns_simulation.simulation_model.'
           'ReflectionConditionFilter')
    def test_get_refl_filter(self, mock_filter):
        mock_filter.SpaceGroup = 1
        mock_filter.StructureFactor = 2
        testv = self.model._get_refl_filter(False)
        self.assertEqual(testv, 1)
        testv = self.model._get_refl_filter(True)
        self.assertEqual(testv, 2)

    @patch('mantidqtinterfaces.dns_simulation.simulation_model.'
           'DNSSimulationModel._get_refl_filter')
    @patch('mantidqtinterfaces.dns_simulation.simulation_model.'
           'ReflectionGenerator')
    @patch('mantidqtinterfaces.dns_simulation.simulation_model.'
           'sim_help.max_q')
    def test__get_filtered_hkls(self, mock_max_q, mock_gen, mock_filter):
        self.model._generator = mock_gen
        mock_filter.return_value = 2
        mock_max_q.return_value = 3
        mock_gen.getUniqueHKLsUsingFilter.return_value = [1]
        mock_gen.getHKLsUsingFilter.return_value = [2]
        testv = self.model._get_filtered_hkls(4.7, '123.cif')
        mock_filter.assert_called_once_with('123.cif')
        mock_max_q.assert_called_once_with(4.7)
        mock_gen.getUniqueHKLsUsingFilter.assert_called_once_with(1 / 3, 100,
                                                                  2)
        mock_gen.getHKLsUsingFilter.assert_called_once_with(1 / 3, 100, 2)
        self.assertEqual(testv, [[2], [1]])

    @patch('mantidqtinterfaces.dns_simulation.simulation_model.'
           'ReflectionGenerator')
    def test__get_q_d_val(self, mock_gen):
        self.model._generator = mock_gen
        mock_gen.getDValues.return_value = [1]
        testv = self.model._get_q_d_val([1])
        mock_gen.getDValues.assert_called_once_with([1])
        self.assertIsInstance(testv, list)
        self.assertIsInstance(testv[1], list)
        self.assertEqual(testv[0], [1])
        self.assertAlmostEqual(testv[1][0], 6.283185307179586)

    @patch('mantidqtinterfaces.dns_simulation.simulation_model.'
           'OrientedLattice')
    @patch('mantidqtinterfaces.dns_simulation.simulation_model.'
           'sim_help.rotate_ub')
    def test_get_ub(self, mock_rotate, mock_ori):
        self.model._orilat = mock_ori
        testv = self.model._get_ub(False, 100)
        self.assertEqual(testv, mock_ori.getUB.return_value)
        mock_ori.getUB.assert_called_once()
        mock_rotate.assert_not_called()
        testv = self.model._get_ub(True, 100)
        mock_rotate.assert_called_once_with(100, mock_ori.getUB.return_value)
        mock_ori.setUB.assert_called_once_with(mock_rotate.return_value)
        self.assertEqual(testv, mock_rotate.return_value)

    def test_get_refls_and_set_orientation(self):
        options = {'hkl1_v': [1, 0, 0],
                   'hkl2_v': [0, 1, 0],
                   'wavelength': 4,
                   'cif_filename': '',
                   'det_rot': -10,
                   'det_number': 3,
                   'fix_omega': True,
                   'omega_offset': 3,
                   'spacegroup': '148',
                   'a': 3,
                   'b': 3,
                   'c': 20,
                   'alpha': 90,
                   'beta': 90,
                   'gamma': 120,
                   'cifset': False,
                   }
        testv = self.model.get_refls_and_set_orientation(options)
        self.assertEqual(len(testv), 29)
        first_dict = {'hkl': sim_help.list_to_v3d([1, 0, 4]), 'unique': True,
                      'q': 2.7253974323168175, 'fs': 1, 'd': 2.30541983810352,
                      'tth': 120.34364613905412, 'fs_lc': 1, 'h': 1.0,
                      'k': 0.0, 'l': 4.0,
                      'equivalents': [sim_help.list_to_v3d([1, 0, 4]),
                                      sim_help.list_to_v3d([1, -1, -4]),
                                      sim_help.list_to_v3d([0, 1, -4]),
                                      sim_help.list_to_v3d([0, -1, 4]),
                                      sim_help.list_to_v3d([-1, 1, 4]),
                                      sim_help.list_to_v3d([-1, 0, -4])],
                      'mult': 6, 'diff': 95.34364613905412,
                      'det_rot': -5.3436461390541155, 'channel': 23.0,
                      'inplane': False, 'omega': -32.48998764184763,
                      'sample_rot': -37.83363378090174}
        self.assertEqual(testv[0], first_dict)

    @patch('mantidqtinterfaces.dns_simulation.simulation_model.'
           'CreateWorkspace')
    def test__create_sim_ws(self, mock_cws):
        testv = self.model._create_sim_ws()
        mock_cws.assert_called_once_with(OutputWorkspace='__sim_ws',
                                         DataX=[0],
                                         DataY=[0],
                                         NSpec=1,
                                         UnitX='Degrees')
        self.assertEqual(testv, mock_cws.return_value)

    def test_get_orilat(self):
        self.model._orilat = 1
        testv = self.model.get_orilat()
        self.assertEqual(testv, 1)

    @patch('mantidqtinterfaces.dns_simulation.simulation_model.'
           'OrientedLattice')
    @patch('mantidqtinterfaces.dns_simulation.simulation_model.'
           'LoadCIF')
    @patch('mantidqtinterfaces.dns_simulation.simulation_model.'
           'DNSSimulationModel._create_sim_ws')
    def test_load_cif(self, mock_cws, mock_loadcif, mock_orilat):
        testv = self.model.load_cif('123.cif')
        mock_cws.assert_called_once()
        mock_loadcif.assert_called_once_with(mock_cws.return_value, '123.cif')
        mock_cs = mock_cws.return_value.sample.return_value.getCrystalStructure
        mock_sg = mock_cs.return_value.getSpaceGroup
        mock_cs.assert_called_once()
        mock_uc = mock_cs.return_value.getUnitCell
        mock_uc.assert_called_once()
        self.assertEqual(mock_orilat.call_count, 2)
        mock_orilat.assert_called_with(mock_uc.return_value)
        self.assertIsNotNone(self.model._orilat)
        self.assertIsNotNone(self.model._non_rot_lat)
        self.assertIsInstance(testv, dict)
        self.assertEqual(
            testv['a'], mock_uc.return_value.a.return_value)
        self.assertEqual(
            testv['b'], mock_uc.return_value.b.return_value)
        self.assertEqual(
            testv['c'], mock_uc.return_value.c.return_value)
        self.assertEqual(
            testv['alpha'], mock_uc.return_value.alpha.return_value)
        self.assertEqual(
            testv['beta'], mock_uc.return_value.beta.return_value)
        self.assertEqual(
            testv['gamma'], mock_uc.return_value.gamma.return_value)
        self.assertEqual(
            testv['spacegroup'], mock_sg.return_value.getHMSymbol.return_value)

    @patch('mantidqtinterfaces.dns_simulation.simulation_model.'
           'sim_help.return_qx_qx_inplane_refl')
    def test_return_reflections_in_map(self, mock_return_q):
        self.model._orilat = 1
        testv = self.model.return_reflections_in_map([1, 2, 3], [4, 5, 6], [])
        mock_return_q.assert_called_once_with(1, [1, 2, 3], [4, 5, 6], [])
        self.assertEqual(testv, mock_return_q.return_value)

    @patch('mantidqtinterfaces.dns_simulation.simulation_model.'
           'SpaceGroupFactory')
    def test__get_hm_spacegroup(self, mock_sgf):
        mock_sgf.subscribedSpaceGroupSymbols.return_value = ['a']
        testv = self.model._get_hm_spacegroup('abc')
        self.assertEqual(testv, 'abc')
        mock_sgf.subscribedSpaceGroupSymbols.assert_not_called()
        testv = self.model._get_hm_spacegroup('1')
        mock_sgf.subscribedSpaceGroupSymbols.assert_called_once_with(1)
        self.assertEqual(testv, 'a')

    def test__get_atom_str(self):
        self.model._cryst = mock.Mock()
        self.model._cryst.getScatterers.return_value = ["a", "b"]
        testv = self.model._get_atom_str('')
        self.assertEqual(testv, "Si 0 0 1 1.0 0.01")
        self.model._cryst.getScatterers.assert_not_called()
        testv = self.model._get_atom_str('123')
        self.model._cryst.getScatterers.assert_called_once()
        self.assertEqual(testv, "a;b")

    def test__get_cell_str(self):
        mock_unitcell = mock.Mock()
        mock_unitcell.a.return_value = 'a'
        mock_unitcell.b.return_value = 'b'
        mock_unitcell.c.return_value = 'c'
        mock_unitcell.alpha.return_value = 'alpha'
        mock_unitcell.beta.return_value = 'beta'
        mock_unitcell.gamma.return_value = 'gamma'
        testv = self.model._get_cell_str(mock_unitcell)
        self.assertEqual(testv, 'a b c alpha beta gamma')

    @patch('mantidqtinterfaces.dns_simulation.simulation_model.'
           'DNSSimulationModel._get_cell_str')
    @patch('mantidqtinterfaces.dns_simulation.simulation_model.'
           'DNSSimulationModel._get_hm_spacegroup')
    @patch('mantidqtinterfaces.dns_simulation.simulation_model.'
           'UnitCell')
    @patch('mantidqtinterfaces.dns_simulation.simulation_model.'
           'CrystalStructure')
    @patch('mantidqtinterfaces.dns_simulation.simulation_model.'
           'OrientedLattice')
    def test__setcellfromparameters(self, mock_orilat, mock_cs, mock_uc,
                                    mock_hm, mock_cellstr):
        options = {'spacegroup': 148,
                   'a': 1,
                   'b': 2,
                   'c': 3,
                   'alpha': 90,
                   'beta': 90,
                   'gamma': 120,
                   'cifset': True,
                   'cif_filename': ''}
        self.model._setcellfromparameters(options)
        mock_orilat.assert_not_called()
        options['cifset'] = False
        self.model._setcellfromparameters(options)
        mock_uc.assert_called_once_with(1, 2, 3, 90, 90, 120, Unit=0)
        self.assertEqual(mock_orilat.call_count, 2)
        mock_orilat.assert_called_with(mock_uc.return_value)
        mock_cs.assert_called_once_with(mock_cellstr.return_value,
                                        mock_hm.return_value,
                                        "Si 0 0 1 1.0 0.01")
        self.assertEqual(self.model._orilat, mock_orilat.return_value)
        self.assertEqual(self.model._non_rot_lat, mock_orilat.return_value)
        self.assertEqual(self.model._cryst, mock_cs.return_value)

    @patch('mantidqtinterfaces.dns_simulation.simulation_model.'
           'sim_help.ki_from_wavelength')
    def test_get_ki(self, mock_wl):
        testv = self.model.get_ki(4.74)
        mock_wl.assert_called_once_with(4.74)
        self.assertEqual(testv, mock_wl.return_value)

    def test_validate_hkl(self):
        testv = self.model.validate_hkl([1, 2, 3], [4, 5, 6])
        self.assertEqual(testv, [True, ''])
        testv = self.model.validate_hkl([None, 2, 3], [4, 5, 6])
        self.assertEqual(testv, [False, 'Could not parse hkl, enter 3 '
                                        'comma seperated numbers.'])
        testv = self.model.validate_hkl([1, 2, 3], [None, 5, 6])
        self.assertEqual(testv, [False, 'Could not parse hkl, enter 3 '
                                        'comma seperated numbers.'])
        testv = self.model.validate_hkl([0, 0, 1], [0, 0, 2])
        self.assertEqual(testv, [False, 'hkl1 cannot be paralell hkl2'])

    @patch('mantidqtinterfaces.dns_simulation.simulation_model.'
           'convert_hkl_string_to_float')
    def test_get_hkl_vector_dict(self, mock_hkl):
        testv = self.model.get_hkl_vector_dict([1, 2, 3], [4, 5, 6])
        self.assertIsInstance(testv, dict)
        self.assertEqual(len(testv), 2)
        self.assertEqual(mock_hkl.call_count, 2)
        mock_hkl.assert_called_with([4, 5, 6])
        self.assertEqual(testv['hkl1_v'], mock_hkl.return_value)
        self.assertEqual(testv['hkl2_v'], mock_hkl.return_value)

    def test_get_oof_from_ident(self):
        testv = self.model.get_oof_from_ident(1, 2, 3, 4)
        self.assertEqual(testv, 2)


if __name__ == '__main__':
    unittest.main()

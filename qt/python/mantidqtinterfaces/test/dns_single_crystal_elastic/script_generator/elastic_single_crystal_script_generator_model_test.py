# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
"""
Common Presenter for DNS Script generators
"""
import unittest
from unittest import mock
from unittest.mock import patch

from mantidqtinterfaces.dns_powder_elastic.data_structures.dns_dataset \
    import DNSDataset
from mantidqtinterfaces.dns_powder_tof.script_generator. \
    common_script_generator_model import DNSScriptGeneratorModel
from mantidqtinterfaces.dns_powder_tof.data_structures.dns_obs_model \
    import DNSObsModel
from mantidqtinterfaces.dns_single_crystal_elastic.script_generator. \
    elastic_single_crystal_script_generator_model import DNSElasticSCScriptGeneratorModel
from mantidqtinterfaces.dns_powder_tof.helpers.helpers_for_testing import \
    get_fake_elastic_datadic, get_elastic_standard_datadic, \
    get_fake_elastic_single_crystal_options


class DNSElasticSCScriptGeneratorModelTest(unittest.TestCase):
    # pylint: disable=protected-access, too-many-public-methods
    parent = None
    model = None
    sample_data = None
    standard_data = None
    fake_workspace = None

    @classmethod
    def setUpClass(cls):
        cls.parent = mock.Mock()
        cls.model = DNSElasticSCScriptGeneratorModel(parent=cls.parent)
        cls.sample_data = mock.create_autospec(DNSDataset)
        cls.model._sample_data = cls.sample_data
        cls.sample_data.datadic = get_fake_elastic_datadic()
        cls.sample_data.create_subtract.return_value = ['knso_x_sf']
        cls.sample_data.format_dataset.return_value = '12345'
        cls.sample_data.fields = []
        cls.sample_data.ttheta = mock.Mock()
        cls.sample_data.ttheta.bin_edge_min = 1
        cls.sample_data.ttheta.bin_edge_max = 6
        cls.sample_data.ttheta.nbins = 5
        cls.sample_data.ttheta.range = [1, 2, 3]
        cls.sample_data.omega = mock.Mock()
        cls.sample_data.omega.bin_edge_min = 3
        cls.sample_data.omega.bin_edge_max = 5
        cls.sample_data.omega.nbins = 3
        cls.sample_data.omega.range = [4, 5, 6]
        cls.sample_data.script_name = 'script.py'
        cls.sample_data.banks = [1, 2, 3]
        cls.standard_data = mock.create_autospec(DNSDataset)
        cls.model._standard_data = cls.standard_data
        cls.standard_data.datadic = get_elastic_standard_datadic()
        cls.standard_data.format_dataset.return_value = '123456'
        cls.standard_data.fields = []
        cls.standard_data.ttheta = mock.Mock()
        cls.standard_data.ttheta.bin_edge_min = 1
        cls.standard_data.ttheta.bin_edge_max = 6
        cls.standard_data.ttheta.nbins = 5
        cls.standard_data.banks = [1, 2, 3]
        cls.standard_data.omega = mock.Mock()
        cls.standard_data.omega.bin_edge_min = 3
        cls.standard_data.omega.bin_edge_max = 5
        cls.standard_data.omega.nbins = 3
        cls.standard_data.omega.range = [4, 5, 6]
        cls.standard_data.script_name = '123.txt'

        cls.fake_workspace = mock.Mock()
        cls.fake_workspace.getErrorSquaredArray.return_value = 1
        cls.fake_workspace.getSignalArray.return_value = 4

    def test___init__(self):
        self.assertIsInstance(self.model, DNSElasticSCScriptGeneratorModel)
        self.assertIsInstance(self.model, DNSScriptGeneratorModel)
        self.assertIsInstance(self.model, DNSObsModel)
        self.assertTrue(hasattr(self.model, '_script'))
        self.assertTrue(hasattr(self.model, '_data_arrays'))
        self.assertTrue(hasattr(self.model, '_plot_list'))
        self.assertTrue(hasattr(self.model, '_sample_data'))
        self.assertTrue(hasattr(self.model, '_standard_data'))
        self.assertTrue(hasattr(self.model, '_loop'))
        self.assertTrue(hasattr(self.model, '_spacing'))
        self.assertTrue(hasattr(self.model, '_vana_correction'))
        self.assertTrue(hasattr(self.model, '_nicr_correction'))
        self.assertTrue(hasattr(self.model, '_sample_background_correction'))
        self.assertTrue(hasattr(self.model, '_background_factor'))
        self.assertTrue(hasattr(self.model, '_ignore_vana'))
        self.assertTrue(hasattr(self.model, '_sum_sf_nsf'))
        self.assertTrue(hasattr(self.model, '_non_magnetic'))
        self.assertTrue(hasattr(self.model, '_corrections'))
        self.assertTrue(hasattr(self.model, '_export_path'))
        self.assertTrue(hasattr(self.model, '_ascii'))
        self.assertTrue(hasattr(self.model, '_nexus'))
        self.assertTrue(hasattr(self.model, '_norm'))

    @patch(
        'mantidqtinterfaces.dns_single_crystal_elastic.'
        'script_generator.elastic_single_crystal_script_generator_mod'
        'el.DNSDataset')
    def test_script_maker(self, mock_dns_dataset):
        mock_dns_dataset.return_value = self.standard_data
        options = {
            'corrections': 1,
            'det_efficency': 1,
            'flipping_ratio': 1,
            'substract_background_from_sample': 1,
            'background_factor': 1,
            'ignore_vana_fields': 1,
            'sum_vana_sf_nsf': 1,
            'separation': 1,
            'separation_xyz': 1,
            'separation_coh_inc': 1,
            'norm_monitor': 1
        }
        options.update(get_fake_elastic_single_crystal_options())
        fselector = {'full_data': [], 'standard_data': []}
        paths = {
            'data_dir': '12',
            'standards_dir': '13',
            'export_dir': '14',
            'ascii': True,
            'nexus': True,
            'export': True
        }
        testv = self.model.script_maker(options, paths, fselector)
        self.assertTrue(self.model._vana_correction)
        self.assertTrue(self.model._nicr_correction)
        self.assertTrue(self.model._sample_background_correction)
        self.assertTrue(self.model._background_factor)
        self.assertTrue(self.model._non_magnetic)
        self.assertTrue(self.model._xyz)
        self.assertTrue(self.model._corrections)
        self.assertTrue(self.model._ascii)
        self.assertTrue(self.model._nexus)
        self.assertEqual(self.model._ignore_vana, '1')
        self.assertEqual(self.model._sum_sf_nsf, '1')
        self.assertEqual(self.model._export_path, '14')
        self.assertEqual(self.model._norm, 'monitor')
        self.assertIsInstance(testv, tuple)
        self.assertIsInstance(testv[0], list)
        lines = [76, 113, 91, 65, 49, 0, 20, 22, 0, 307, 0, 109, 51, 71, 0, 45,
                 130, 0,
                 21, 250, 129]
        for i, tv in enumerate(testv[0]):
            self.assertEqual(len(tv), lines[i])
        options = {
            'corrections': 0,
            'det_efficency': 0,
            'flipping_ratio': 0,
            'substract_background_from_sample': 0,
            'background_factor': 0,
            'ignore_vana_fields': 0,
            'sum_vana_sf_nsf': 0,
            'separation': 0,
            'separation_xyz': 0,
            'separation_coh_inc': 0,
            'norm_monitor': 0
        }
        options.update(get_fake_elastic_single_crystal_options())
        fselector = {'full_data': [], 'standard_data': []}
        paths = {
            'data_dir': '',
            'standards_dir': '',
            'export_dir': '',
            'ascii': False,
            'nexus': False,
            'export': False
        }
        testv = self.model.script_maker(options, paths, fselector)
        lines = [76, 113, 91, 65, 49, 0, 20, 0, 0, 304, 0, 109, 51, 0]
        for i, tv in enumerate(testv[0]):
            self.assertEqual(len(tv), lines[i])
        options = {
            'corrections': 0,
            'det_efficency': 0,
            'flipping_ratio': 0,
            'substract_background_from_sample': 0,
            'background_factor': 0,
            'ignore_vana_fields': 0,
            'sum_vana_sf_nsf': 0,
            'separation': 0,
            'separation_xyz': 0,
            'separation_coh_inc': 0,
            'norm_monitor': 0
        }
        options.update(get_fake_elastic_single_crystal_options())
        testv = self.model.script_maker(options, paths, fselector)
        self.assertFalse(self.model._vana_correction)
        self.assertFalse(self.model._nicr_correction)
        self.assertFalse(self.model._sample_background_correction)
        self.assertFalse(self.model._background_factor)
        self.assertFalse(self.model._non_magnetic)
        self.assertFalse(self.model._xyz)
        self.assertFalse(self.model._corrections)
        self.assertFalse(self.model._ascii)
        self.assertFalse(self.model._nexus)
        self.assertEqual(self.model._ignore_vana, '0')
        self.assertEqual(self.model._sum_sf_nsf, '0')
        self.assertEqual(self.model._export_path, '')
        self.assertEqual(self.model._norm, 'time')
        lines = [76, 113, 91, 65, 49, 0, 20, 0, 0, 304, 0, 109, 51, 0]
        for i, tv in enumerate(testv[0]):
            self.assertEqual(len(tv), lines[i])
        options = {
            'corrections': 1,
            'det_efficency': 0,
            'flipping_ratio': 0,
            'substract_background_from_sample': 0,
            'background_factor': 1,
            'ignore_vana_fields': 1,
            'sum_vana_sf_nsf': 1,
            'separation': 1,
            'separation_xyz': 0,
            'separation_coh_inc': 0,
            'norm_monitor': 1
        }
        options.update(get_fake_elastic_single_crystal_options())
        fselector = {'full_data': [], 'standard_data': []}
        paths = {
            'data_dir': '12',
            'standards_dir': '13',
            'export_dir': '14',
            'ascii': False,
            'nexus': False,
            'export': True
        }
        testv = self.model.script_maker(options, paths, fselector)
        self.assertEqual(len(testv[0][0]), 76)
        self.assertFalse(self.model._vana_correction)
        self.assertFalse(self.model._nicr_correction)
        self.assertFalse(self.model._sample_background_correction)
        self.assertFalse(self.model._non_magnetic)
        self.assertFalse(self.model._xyz)
        self.assertFalse(self.model._corrections)
        self.assertFalse(self.model._ascii)
        self.assertFalse(self.model._nexus)
        lines = [76, 113, 91, 65, 49, 0, 20, 0, 0, 307, 0, 109, 51, 0]
        for i, tv in enumerate(testv[0]):
            self.assertEqual(len(tv), lines[i])
        options = {
            'corrections': 0,
            'det_efficency': 1,
            'flipping_ratio': 1,
            'substract_background_from_sample': 1,
            'background_factor': 1,
            'ignore_vana_fields': 1,
            'sum_vana_sf_nsf': 1,
            'separation': 0,
            'separation_xyz': 1,
            'separation_coh_inc': 1,
            'norm_monitor': 1
        }
        options.update(get_fake_elastic_single_crystal_options())
        fselector = {'full_data': [], 'standard_data': []}
        paths = {
            'data_dir': '12',
            'standards_dir': '13',
            'export_dir': '14',
            'ascii': True,
            'nexus': True,
            'export': False
        }
        testv = self.model.script_maker(options, paths, fselector)
        self.assertFalse(self.model._vana_correction)
        self.assertFalse(self.model._nicr_correction)
        self.assertFalse(self.model._sample_background_correction)
        self.assertFalse(self.model._non_magnetic)
        self.assertFalse(self.model._xyz)
        self.assertFalse(self.model._corrections)
        self.assertFalse(self.model._ascii)
        self.assertFalse(self.model._nexus)
        lines = [76, 113, 91, 65, 49, 0, 20, 0, 0, 307, 0, 109, 51, 0]
        for i, tv in enumerate(testv[0]):
            self.assertEqual(len(tv), lines[i])
        paths = {
            'data_dir': '12',
            'standards_dir': '13',
            'export_dir': '',
            'ascii': True,
            'nexus': True,
            'export': True
        }
        testv = self.model.script_maker(options, paths, fselector)
        self.assertFalse(self.model._ascii)
        self.assertFalse(self.model._nexus)
        lines = [76, 113, 91, 65, 49, 0, 20, 0, 0, 307, 0, 109, 51, 0]
        for i, tv in enumerate(testv[0]):
            self.assertEqual(len(tv), lines[i])

    @patch(
        'mantidqtinterfaces.dns_single_crystal_elastic.'
        'script_generator.elastic_single_crystal_script_generator_model'
        '.DNSDataset')
    def test_setup_sample_data(self, mock_dns_dataset):
        self.model._sample_data = None
        mock_dns_dataset.return_value = self.sample_data
        self.model._setup_sample_data({'data_dir': '123'}, {'full_data': []})
        self.assertEqual(self.model._sample_data, self.sample_data)
        self.assertEqual(self.model._plot_list, ['knso_x_sf'])

    @patch(
        'mantidqtinterfaces.dns_single_crystal_elastic.'
        'script_generator.elastic_single_crystal_script_generator_model'
        '.DNSDataset')
    def test_setup_standard_data(self, mock_dns_dataset):
        self.model._standard_data = None
        self.model._corrections = True

        mock_dns_dataset.return_value = self.standard_data
        self.model._setup_standard_data({'standards_dir': '123'},
                                        {'standard_data': []})
        self.assertEqual(self.model._standard_data, self.standard_data)

    def test_set_loop(self):
        self.model._loop = None
        self.model._spacing = None
        self.model._set_loop()
        self.assertEqual(
            self.model._loop,
            "for sample, workspacelist in wss_sample.items(): \n    for work"
            "space in workspacelist:")
        self.assertEqual(self.model._spacing, "\n" + " " * 8)
        self.model._sample_data = {'123': 1}
        self.model._set_loop()
        self.assertEqual(self.model._loop,
                         "for workspace in wss_sample['123']:")
        self.assertEqual(self.model._spacing, "\n" + " " * 4)
        self.model._sample_data = self.sample_data

    def test_get_header_lines(self):
        testv = self.model._get_header_lines()
        self.assertIsInstance(testv, list)
        self.assertEqual(len(testv), 6)
        self.assertEqual(
            testv[0],
            'from mantidqtinterfaces.dns_single_crystal_elastic.'
            'scripts.md_single_crystal_elastic import load_all'
        )
        self.assertEqual(
            testv[1],
            'from mantidqtinterfaces.dns_single_crystal_elastic.'
            'scripts.md_single_crystal_elastic import '
            'vanadium_correction, fliping_ratio_correction')
        self.assertEqual(
            testv[2],
            'from mantidqtinterfaces.dns_single_crystal_elastic.'
            'scripts.md_single_crystal_elastic import'
            ' background_subtraction', )
        self.assertEqual(
            testv[3],
            "from mantid.simpleapi import ConvertMDHistoToMatrixWorkspace, mtd"
        )
        self.assertEqual(testv[4],
                         "from mantid.simpleapi import SaveAscii, SaveNexus")
        self.assertEqual(testv[5], "")

    def test_get_sample_data_lines(self):
        testv = self.model._get_sample_data_lines()
        self.assertEqual(testv, ['sample_data = 12345'])

    def test_get_standard_data_lines(self):
        self.model._corrections = True
        testv = self.model._get_standard_data_lines()
        self.assertEqual(testv, ['standard_data = 123456'])
        self.model._corrections = False
        testv = self.model._get_standard_data_lines()
        self.assertEqual(testv, [''])

    def test__get_param_lines(self):
        self.model._norm = 'time'
        testv = self.model._get_param_lines(get_fake_elastic_single_crystal_options())
        self.assertIsInstance(testv, list)
        self.assertEqual(len(testv), 3)
        comparestring = ("params = {'a' : 2, \n          'b' : 3,\n          "
                         "'c' : 4,\n          'alpha' : 78,\n          'beta'"
                         "  : 86,\n          'gamma' : 85,\n          'hkl1' "
                         " : '1,2,3',\n          'hkl2'  : '2,3,4',\n        "
                         "  'omega_offset' : 0,\n          'norm_to' : "
                         "'time',\n          'dx' : ' 1.0000',\n        "
                         "  'dy' : ' 2.0000',}")
        self.assertEqual(testv[1], comparestring)

    def test__get_binning_lines(self):
        testv = self.model._get_binning_lines()
        testl = ["binning = {'twoTheta' : [1.000, 6.000, 5],\n       "
                 "    'Omega':  [3.000, 5.000, 3]} # min, max, number_of_bins"]
        self.assertEqual(testv, testl)

    def test_get_load_data_lines(self):
        self.model._corrections = True
        testv = self.model._get_load_data_lines()
        self.assertEqual(testv, ['wss_sample = load_all(sample_data, binning,'
                                 ' params)', 'wss_standard = load_all(standar'
                                             'd_data, binning, params, standar'
                                             'd=True,)', ''])
        self.model._corrections = False
        testv = self.model._get_load_data_lines()
        self.assertEqual(testv, [
            "wss_sample = load_all(sample_data, binning, params)", ''
        ])

    def test__get_bg_corr_lines(self):
        self.model._vana_correction = True
        testv = self.model._get_bg_corr_lines()
        self.assertEqual(testv, [
            '# substract background from vanadium and nicr',
            'for sample, workspacelist in wss_standard.items(): \n    '
            'for workspace in workspacelist:\n        background_substr'
            'action(workspace)', ''
        ])
        self.model._vana_correction = False
        self.model._nicr_correction = False
        testv = self.model._get_bg_corr_lines()
        self.assertEqual(testv, [])

    def test__return_sample_bg_string(self):
        self.model._spacing = '  '
        self.model._background_factor = '123'
        testv = self.model._return_sample_bg_string()
        self.assertEqual(testv, '  background_subtraction(workspace,'
                                ' factor=123)')

    def test__return_sample_vanac_strinf(self):
        self.model._spacing = '  '
        self.model._sum_sf_nsf = 1
        self.model._ignore_vana = 2
        testv = self.model._return_sample_vanac_strinf()
        self.assertEqual(testv, "  vanadium_correction(workspace,  vana_set=s"
                                "tandard_data['vana'], ignore_vana_fields=2,"
                                " sum_vana_sf_nsf=1)")

    def test__get_vanac_lines(self):
        self.model._background_factor = '123'
        self.model._loop = 'fo:'
        self.model._spacing = '  '
        self.model._sum_sf_nsf = 1
        self.model._ignore_vana = 2
        self.model._sample_background_correction = False
        self.model._vana_correction = False
        compare = ['# correct sample data',
                   "fo:  background_subtraction(workspace, factor=123)  vanad"
                   "ium_correction(workspace,  vana_set=standard_data['vana'],"
                   " ignore_vana_fields=2, sum_vana_sf_nsf=1)"]

        testv = self.model._get_vanac_lines()
        self.assertEqual(testv, [])
        self.model._sample_background_correction = True
        testv = self.model._get_vanac_lines()
        self.assertEqual(testv, compare)
        self.model._sample_background_correction = False
        self.model._vana_correction = True
        testv = self.model._get_vanac_lines()
        self.assertEqual(testv, compare)

    def test__get_nicrc_lines(self):
        self.model._nicr_correction = False
        self.model._loop = 'fo:'
        self.model._spacing = '  '
        testv = self.model._get_nicrc_lines()
        self.assertEqual(testv, [])
        self.model._nicr_correction = True
        testv = self.model._get_nicrc_lines()
        self.assertEqual(testv, ['fo:  fliping_ratio_correction(workspace)'])

    @patch(
        'mantidqtinterfaces.dns_single_crystal_elastic.script_generator.'
        'elastic_single_crystal_script_generator_model.mtd')
    def test_get_plotlist(self, mtd):
        mtd.__getitem__.return_value = self.fake_workspace
        self.model._plot_list = ['4p1K_map']
        testv = self.model.get_plotlist()
        self.assertEqual(
            testv, (['4p1K_map'],
                    {'4p1K_map': {'ttheta': [1, 2, 3],
                                  'omega': [4, 5, 6],
                                  'intensity': 4,
                                  'error': 1.0}}))


if __name__ == '__main__':
    unittest.main()

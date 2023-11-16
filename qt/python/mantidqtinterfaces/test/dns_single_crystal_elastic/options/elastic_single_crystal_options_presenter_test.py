# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
import unittest
from unittest import mock
from unittest.mock import patch

from mantidqtinterfaces.dns_powder_tof.data_structures.dns_observer \
    import DNSObserver
from mantidqtinterfaces.dns_powder_tof.options.common_options_presenter \
    import DNSCommonOptionsPresenter
from mantidqtinterfaces.dns_single_crystal_elastic.options.elastic_single_crystal_options_presenter \
    import DNSElasticSCOptionsPresenter


class DNSElasticSCOptionsPresenterTest(unittest.TestCase):
    # pylint: disable=protected-access, too-many-public-methods
    model = None
    view = None
    parent = None
    presenter = None

    @classmethod
    def setUpClass(cls):
        cls.parent = mock.Mock()
        cls.view = mock.Mock()
        cls.model = mock.Mock()
        cls.model.get_dx_dy.return_value = [1, 2]
        cls.view.sig_get_wavelength = mock.Mock()
        cls.view.sig_get_wavelength.connect = mock.Mock()

        cls.view.get_state = mock.Mock(
            return_value={
                'wavelength': 4.74,
                'get_wavelength': False,
                'hkl1': '1,2,3',
                'hkl2': '2,3,4',
                'use_dx_dy': True
            })

        cls.presenter = DNSElasticSCOptionsPresenter(
            view=cls.view,
            model=cls.model,
            parent=cls.parent,
            name='elastic_single_crystal_options')

    def test___init__(self):
        self.assertIsInstance(self.presenter, DNSElasticSCOptionsPresenter)
        self.assertIsInstance(self.presenter, DNSCommonOptionsPresenter)
        self.assertIsInstance(self.presenter, DNSObserver)
        self.view.sig_get_wavelength.connect.assert_called_once_with(
            self.presenter._determine_wavelength)

    @patch('mantidqtinterfaces.dns_single_crystal_elastic.'
           'options.elastic_single_crystal_options_presenter.'
           'DNSElasticSCOptionsPresenter._determine_wavelength')
    def test_process_request(self, mock_wavel):
        self.view.get_state.return_value = {
            'wavelength': 4.74,
            'get_wavelength': False,
            'hkl1': '1,2,3',
            'hkl2': '2,3,4',
            'use_dx_dy': True
        }
        self.presenter.process_request()
        mock_wavel.assert_not_called()
        self.view.get_state.return_value = {
            'wavelength': 4.74,
            'get_wavelength': True,
            'hkl1': '1,2,3',
            'hkl2': '2,3,4',
            'use_dx_dy': True
        }
        self.presenter.process_request()
        mock_wavel.assert_called_once()

    def test_process_commandline_request(self):
        self.view.reset_mock()
        self.presenter.process_commandline_request({'omega_offset': 1})
        self.view.set_single_state_by_name.assert_called_with(
            'omega_offset', 1)

    def test_get_option_dict(self):
        self.view.get_state.return_value = {
            'wavelength': 4.74,
            'get_wavelength': True,
            'hkl1': '[1 2 3]',
            'hkl2': '2,3,4',
            'use_dx_dy': False
        }
        testv = self.presenter.get_option_dict()
        self.assertTrue('dx' in testv)
        self.assertTrue('dy' in testv)
        self.assertTrue(testv['hkl1'], '1,2,3')
        self.assertEqual(len(testv), 7)


if __name__ == '__main__':
    unittest.main()

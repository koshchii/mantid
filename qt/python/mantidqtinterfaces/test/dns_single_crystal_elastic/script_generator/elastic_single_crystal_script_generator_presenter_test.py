# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
"""
DNS script generator for elastic powder data
"""
import unittest
from unittest import mock
from collections import OrderedDict
from mantidqtinterfaces.dns_powder_tof.script_generator.\
    common_script_generator_presenter import DNSScriptGeneratorPresenter
from mantidqtinterfaces.dns_powder_tof.data_structures.dns_observer import \
    DNSObserver
from mantidqtinterfaces.dns_single_crystal_elastic.script_generator.\
    elastic_single_crystal_script_generator_presenter \
    import DNSElasticSCScriptGeneratorPresenter
from mantidqtinterfaces.dns_powder_tof.script_generator.\
    common_script_generator_view import DNSScriptGeneratorView
from mantidqtinterfaces.dns_single_crystal_elastic.script_generator.\
    elastic_single_crystal_script_generator_model import DNSElasticSCScriptGeneratorModel


class DNSElasticSCScriptGeneratorPresenterTest(unittest.TestCase):
    # pylint: disable=protected-access, too-many-public-methods
    view = None
    model = None
    parent = None

    @classmethod
    def setUpClass(cls):
        cls.parent = mock.Mock()
        cls.view = mock.create_autospec(DNSScriptGeneratorView)
        cls.view.sig_progress_canceled.connect = mock.Mock()
        cls.view.sig_generate_script.connect = mock.Mock()
        cls.view.get_state.return_value = {}
        cls.model = mock.create_autospec(DNSElasticSCScriptGeneratorModel)
        cls.model.get_plotlist.return_value = [[1], {2: 3}]
        cls.presenter = DNSElasticSCScriptGeneratorPresenter(
            view=cls.view,
            model=cls.model,
            name='elastic_single_crystal_script_generator',
            parent=cls.parent)

    def setUp(self):
        self.view.reset_mock()
        self.model.reset_mock()

    def test___init__(self):
        self.assertIsInstance(self.presenter,
                              DNSElasticSCScriptGeneratorPresenter)
        self.assertIsInstance(self.presenter, DNSScriptGeneratorPresenter)
        self.assertIsInstance(self.presenter, DNSObserver)
        self.assertTrue(hasattr(self.presenter, '_plotlist'))
        self.assertTrue(hasattr(self.presenter, '_data_arrays'))

    def test__finish_script_run(self):
        self.presenter._finish_script_run()
        self.model.get_plotlist.assert_called_once()
        self.assertEqual(self.presenter._plotlist, [1])
        self.assertEqual(self.presenter._data_arrays, {2: 3})

    def test_get_option_dict(self):
        self.presenter._plotlist = []
        self.presenter._data_arrays = {}
        testv = self.presenter.get_option_dict()
        self.view.get_state.assert_called_once()
        self.assertEqual(testv, OrderedDict([('script_path', ''),
                                             ('script_number', 0),
                                             ('script_text', ''),
                                             ('plotlist', []),
                                             ('data_arrays', {})]))


if __name__ == '__main__':
    unittest.main()

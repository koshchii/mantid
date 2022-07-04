# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
import unittest
from unittest import mock
from mantidqtinterfaces.dns_simulation.simulation_sub_common_presenter import \
    DNSSimulationSubCommonPresenter


class DNSSimulationSubCommonPresenterTest(unittest.TestCase):
    model = None
    view = None

    @classmethod
    def setUpClass(cls):
        cls.view = mock.Mock()
        cls.model = mock.Mock()
        cls.presenter = DNSSimulationSubCommonPresenter(view=cls.view,
                                                        model=cls.model)

    def test___init__(self):
        self.assertIsInstance(self.presenter, DNSSimulationSubCommonPresenter)
        self.assertEqual(self.presenter.view, self.view)
        self.assertEqual(self.presenter.model, self.model)

    # def test_process_request(self):
    #     pass

    def test_get_option_dict(self):
        testv = self.presenter.get_option_dict()
        self.view.get_state.assert_called_once()
        self.assertEqual(testv, self.view.get_state.return_value)
        self.presenter.view = None
        testv = self.presenter.get_option_dict()
        self.presenter.view = self.view
        self.assertEqual(testv, {})


if __name__ == '__main__':
    unittest.main()

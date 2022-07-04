# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
import unittest
from unittest import mock
from unittest.mock import patch

from mantidqtinterfaces.dns_simulation.simulation_sub_powder_presenter import \
    DNSSimulationSubPowderPresenter

from mantidqtinterfaces.dns_simulation.simulation_sub_powder_model import \
    DNSSimulationSubPowderModel

from mantidqtinterfaces.dns_simulation.simulation_sub_powder_view import \
    DNSSimulationSubPowderView

from mantidqtinterfaces.dns_simulation.simulation_sub_common_presenter import \
    DNSSimulationSubCommonPresenter


class DNSSimulationSubPowderPresenterTest(unittest.TestCase):
    # pylint: disable=protected-access

    view = None
    model = None
    parent = None

    @classmethod
    def setUpClass(cls):
        cls.parent = mock.Mock()
        cls.view = mock.create_autospec(DNSSimulationSubPowderView)
        cls.view.sig_powderplot_clicked.connect = mock.Mock()

        cls.model = mock.create_autospec(DNSSimulationSubPowderModel)
        cls.presenter = DNSSimulationSubPowderPresenter(parent=cls.parent,
                                                        view=cls.view,
                                                        model=cls.model,
                                                        )

    def setUp(self):
        self.view.reset_mock()
        self.model.reset_mock()

    def test__init(self):
        self.assertIsInstance(self.presenter, DNSSimulationSubPowderPresenter)
        self.assertIsInstance(self.presenter, DNSSimulationSubCommonPresenter)
        self.assertTrue(hasattr(self.presenter, '_sub_dict'))

    @patch('mantidqtinterfaces.dns_simulation.simulation_sub_powder_presenter.'
           'DNSSimulationSubPowderPresenter._powderplot')
    def test_process_request(self, mock_powderplot):
        self.presenter.process_request({})
        self.assertEqual(self.presenter._sub_dict, {})
        mock_powderplot.assert_called_once()

    @patch('mantidqtinterfaces.dns_simulation.simulation_sub_powder_presenter.'
           'DNSSimulationSubPowderPresenter._annotate_reflections')
    def test_powderplot(self, mock_anno):
        self.view.get_state.return_value = {'powder_start': 0,
                                            'powder_end': 3,
                                            'shift': 1}
        self.model.create_powder_profile.return_value = 4, 5
        self.presenter._sub_dict = None
        self.presenter._powderplot()
        self.model.create_powder_profile.assert_not_called()
        self.presenter._sub_dict = {'refls': []}
        self.presenter._powderplot()
        self.model.create_powder_profile.assert_called_once_with([], 0, 3, 1)
        self.model.get_annotation_list.assert_called_once_with([], 0, 3, 1, 5)
        self.view.start_powderplot.assert_called_once_with(4, 5)
        mock_anno.assert_called_once_with(
            self.model.get_annotation_list.return_value)
        self.view.finish_powderplot.assert_called_once()

    def test__annotate_reflections(self):
        self.view.get_state.return_value = {'labels': False}
        refl_to_annotate = [[5], [[1, 2, 3]], [6]]
        self.presenter._annotate_reflections(refl_to_annotate)
        self.view.annotate_reflection.assert_not_called()
        self.view.get_state.return_value = {'labels': True}
        self.presenter._annotate_reflections(refl_to_annotate)
        self.view.annotate_reflection.assert_called_once_with(
            '  [1.00, 2.00, 3.00]', 5, 6)


if __name__ == '__main__':
    unittest.main()

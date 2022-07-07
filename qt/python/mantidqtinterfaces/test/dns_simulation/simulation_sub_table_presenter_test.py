# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +

import unittest
from unittest import mock
from unittest.mock import call
from unittest.mock import patch
from mantidqtinterfaces.dns_simulation.simulation_sub_common_presenter import \
    DNSSimulationSubCommonPresenter
from mantidqtinterfaces.dns_simulation.simulation_sub_table_model import \
    DNSSimulationSubTableModel
from mantidqtinterfaces.dns_simulation.simulation_sub_table_presenter import \
    DNSSimulationSubTablePresenter
from mantidqtinterfaces.dns_simulation.simulation_sub_table_view import \
    DNSSimulationSubTableView


class DNSSimulationSubTablePresenterTest(unittest.TestCase):
    # pylint: disable=protected-access

    parent = None
    view = None
    model = None

    @classmethod
    def setUpClass(cls):
        cls.parent = mock.Mock()
        cls.parent.update_progress = mock.Mock()
        cls.view = mock.create_autospec(DNSSimulationSubTableView)
        cls.view.sig_table_item_clicked.connect = mock.Mock()
        cls.model = mock.create_autospec(DNSSimulationSubTableModel)
        cls.presenter = DNSSimulationSubTablePresenter(parent=cls.parent,
                                                       view=cls.view,
                                                       model=cls.model,
                                                       )

    def test___init__(self):
        self.assertIsInstance(self.presenter, DNSSimulationSubTablePresenter)
        self.assertIsInstance(self.presenter, DNSSimulationSubCommonPresenter)
        self.assertTrue(hasattr(self.presenter, '_sub_dict'))

    @patch('mantidqtinterfaces.dns_simulation.simulation_sub_table_presenter.'
           'DNSSimulationSubTablePresenter._write_table')
    def test_process_request(self, mock_wt):
        sub_dict = {'filtered_refls': [],
                    'tth_limit': [1, 2]}
        self.presenter.process_request(sub_dict)
        mock_wt.assert_called_once_with([], [1, 2])

    def test__tableitemdclicked(self):
        self.presenter._table_item_double_clicked(-5, 10)
        parent_presenter = self.parent.parent.presenter
        backcall = parent_presenter.back_call_from_table_item_clicked
        backcall.assert_called_once_with(-5, 10)

    @patch('mantidqtinterfaces.dns_simulation.simulation_sub_table_presenter.'
           'DNSSimulationSubTablePresenter._add_mult_tooltip')
    @patch('mantidqtinterfaces.dns_simulation.simulation_sub_table_presenter.'
           'DNSSimulationSubTablePresenter._color_identified')
    def test__writetable(self, mock_color, mock_tool):
        refl = mock.Mock()
        refl.h = 1
        refl.k = 2
        refl.l = 3
        refl.q = 4
        refl.d = 5
        refl.tth = 6
        refl.fs = 8
        refl.mult = 9
        refl.diff = 10
        refl.det_rot = 11
        refl.channel = 12
        refl.sample_rot = 13
        self.presenter._write_table([refl], 5)
        self.view.start_table.assert_called_once_with(1, 12)
        testcalls = [call(' 1 '), call(' 2 '), call(' 3 '), call(' 4.00 '),
                     call(' 5.00 '), call(' 6.00 '), call(' 8 '), call(' 9 '),
                     call(' 10.00 '), call(' 11.00 '), call(' 12 '),
                     call(' 13.00 ')]
        self.view.create_tableitem.assert_has_calls(testcalls)
        for i in range(12):
            mock_tool.assert_any_call(refl, i)
            self.view.set_tableitem.assert_any_call(0, i)
        self.assertEqual(mock_color.call_count, 12)
        mock_color.assert_called_with(refl, 5)
        self.view.finish_table.assert_called_once()

    def test__color_identified(self):
        refl = mock.Mock()
        refl.diff = 1
        self.presenter._color_identified(refl, 5)
        self.view.set_bg_color.assert_called_once_with(True)
        refl.diff = 10
        self.presenter._color_identified(refl, 5)
        self.view.set_bg_color.assert_called_with(False)

    def test_add_mult_tooltip(self):
        refl = mock.Mock()
        refl.equivalents = [1, 2]
        self.presenter._add_mult_tooltip(refl, 0)
        self.view.set_mult_tooltip.assert_not_called()
        self.presenter._add_mult_tooltip(refl, 7)
        self.view.set_mult_tooltip.assert_called_once_with('[1, 2]')


if __name__ == '__main__':
    unittest.main()

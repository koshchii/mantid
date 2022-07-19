# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
import unittest
from unittest import mock

from mantidqt.gui_helper import get_qapplication

from mantidqtinterfaces.dns_simulation.simulation_sub_table_model import \
    DNSSimulationSubTableModel
from mantidqtinterfaces.dns_simulation.simulation_sub_table_presenter import \
    DNSSimulationSubTablePresenter
from mantidqtinterfaces.dns_simulation.simulation_sub_table_view import \
    DNSSimulationSubTableView
from mantidqtinterfaces.dns_simulation.simulation_sub_table_widget import \
    DNSSimulationSubTableWidget

app, within_mantid = get_qapplication()


class DNSSimulationSubPowderWidgetTest(unittest.TestCase):
    parent = None

    @classmethod
    def setUpClass(cls):
        cls.parent = mock.Mock()
        cls.parent.view = None
        cls.widget = DNSSimulationSubTableWidget(cls.parent)

    def test___init__(self):
        self.assertIsInstance(self.widget, DNSSimulationSubTableWidget)
        self.assertIsInstance(self.widget.view, DNSSimulationSubTableView)
        self.assertIsInstance(self.widget.model, DNSSimulationSubTableModel)
        self.assertIsInstance(self.widget.presenter,
                              DNSSimulationSubTablePresenter)
        self.assertEqual(self.widget.name, "Reflection Table")


if __name__ == '__main__':
    unittest.main()

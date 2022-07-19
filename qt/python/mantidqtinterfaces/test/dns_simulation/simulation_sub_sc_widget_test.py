# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
import unittest
from unittest import mock

from mantidqt.gui_helper import get_qapplication

from mantidqtinterfaces.dns_simulation.simulation_sub_sc_model import \
    DNSSimulationSubScModel
from mantidqtinterfaces.dns_simulation.simulation_sub_sc_presenter import \
    DNSSimulationSubScPresenter
from mantidqtinterfaces.dns_simulation.simulation_sub_sc_view import \
    DNSSimulationSubScView
from mantidqtinterfaces.dns_simulation.simulation_sub_sc_widget import \
    DNSSimulationSubScWidget

app, within_mantid = get_qapplication()


class DNSSimulationSubPowderWidgetTest(unittest.TestCase):
    parent = None

    @classmethod
    def setUpClass(cls):
        cls.parent = mock.Mock()
        cls.parent.view = None
        cls.widget = DNSSimulationSubScWidget(cls.parent)

    def test___init__(self):
        self.assertIsInstance(self.widget, DNSSimulationSubScWidget)
        self.assertIsInstance(self.widget.view, DNSSimulationSubScView)
        self.assertIsInstance(self.widget.model, DNSSimulationSubScModel)
        self.assertIsInstance(self.widget.presenter,
                              DNSSimulationSubScPresenter)
        self.assertEqual(self.widget.name, "Single Crystal Map")


if __name__ == '__main__':
    unittest.main()

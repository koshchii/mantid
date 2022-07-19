# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
import unittest
from unittest import mock

from mantidqt.gui_helper import get_qapplication

from mantidqtinterfaces.dns_powder_tof.data_structures.dns_widget import \
    DNSWidget
from mantidqtinterfaces.dns_simulation.simulation_model import \
    DNSSimulationModel
from mantidqtinterfaces.dns_simulation.simulation_presenter import \
    DNSSimulationPresenter
from mantidqtinterfaces.dns_simulation.simulation_view import DNSSimulationView
from mantidqtinterfaces.dns_simulation.simulation_widget import \
    DNSSimulationWidget

app, within_mantid = get_qapplication()


class DNSSimulationWidgetTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        parent = mock.Mock()
        parent.view = None
        cls.widget = DNSSimulationWidget('simulation', parent)

    def test___init__(self):
        self.assertIsInstance(self.widget, DNSSimulationWidget)
        self.assertIsInstance(self.widget, DNSWidget)
        self.assertIsInstance(self.widget.view, DNSSimulationView)
        self.assertIsInstance(self.widget.model, DNSSimulationModel)
        self.assertIsInstance(self.widget.presenter, DNSSimulationPresenter)


if __name__ == '__main__':
    unittest.main()

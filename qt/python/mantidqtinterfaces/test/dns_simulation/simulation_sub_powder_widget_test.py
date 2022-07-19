# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
import unittest
from unittest import mock

from mantidqt.gui_helper import get_qapplication

from mantidqtinterfaces.dns_simulation.simulation_sub_powder_presenter import \
    DNSSimulationSubPowderPresenter
from mantidqtinterfaces.dns_simulation.simulation_sub_powder_model import \
    DNSSimulationSubPowderModel
from mantidqtinterfaces.dns_simulation.simulation_sub_powder_view import \
    DNSSimulationSubPowderView
from mantidqtinterfaces.dns_simulation.simulation_sub_powder_widget import \
    DNSSimulationSubPowderWidget

app, within_mantid = get_qapplication()


class DNSSimulationSubPowderWidgetTest(unittest.TestCase):
    parent = None

    @classmethod
    def setUpClass(cls):
        cls.parent = mock.Mock()
        cls.parent.view = None
        cls.widget = DNSSimulationSubPowderWidget(cls.parent)

    def test___init__(self):
        self.assertIsInstance(self.widget, DNSSimulationSubPowderWidget)
        self.assertIsInstance(self.widget.view, DNSSimulationSubPowderView)
        self.assertIsInstance(self.widget.model, DNSSimulationSubPowderModel)
        self.assertIsInstance(self.widget.presenter,
                              DNSSimulationSubPowderPresenter)
        self.assertEqual(self.widget.name, "Powder diffraktorgram")


if __name__ == '__main__':
    unittest.main()

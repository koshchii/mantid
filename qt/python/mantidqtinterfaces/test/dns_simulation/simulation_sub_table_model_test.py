# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
import unittest
from unittest import mock

from mantidqtinterfaces.dns_simulation.simulation_sub_table_model import \
    DNSSimulationSubTableModel


class DNSSimulationSubTableModelTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        parent = mock.Mock()
        cls.model = DNSSimulationSubTableModel(parent)

    def test__init(self):
        self.assertIsInstance(self.model, DNSSimulationSubTableModel)


if __name__ == '__main__':
    unittest.main()

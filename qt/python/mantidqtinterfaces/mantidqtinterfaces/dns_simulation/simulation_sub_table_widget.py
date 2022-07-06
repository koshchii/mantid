# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
"""
DNS Simulation widget.
"""

from mantidqtinterfaces.dns_simulation.simulation_sub_table_model import \
    DNSSimulationSubTableModel
from mantidqtinterfaces.dns_simulation.simulation_sub_table_presenter import \
    DNSSimulationSubTablePresenter
from mantidqtinterfaces.dns_simulation.simulation_sub_table_view import \
    DNSSimulationSubTableView


class DNSSimulationSubTableWidget:
    def __init__(self, parent):
        super().__init__()
        self.name = "Reflection Table"
        self.view = DNSSimulationSubTableView(parent=parent.view)
        self.model = DNSSimulationSubTableModel(parent=self)
        self.parent = parent
        self.presenter = DNSSimulationSubTablePresenter(parent=self,
                                                        view=self.view,
                                                        model=self.model,
                                                        )

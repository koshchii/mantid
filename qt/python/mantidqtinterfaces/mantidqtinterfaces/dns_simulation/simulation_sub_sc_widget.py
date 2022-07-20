# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +

"""
DNS Simulation widget.
"""

from mantidqtinterfaces.dns_simulation.simulation_sub_sc_model import \
    DNSSimulationSubScModel
from mantidqtinterfaces.dns_simulation.simulation_sub_sc_presenter import \
    DNSSimulationSubScPresenter
from mantidqtinterfaces.dns_simulation.simulation_sub_sc_view import \
    DNSSimulationSubScView


class DNSSimulationSubScWidget:
    def __init__(self, parent):
        super().__init__()
        self.name = "Single Crystal Map"
        self.view = DNSSimulationSubScView(parent=parent.view)
        self.model = DNSSimulationSubScModel(parent=self)
        self.presenter = DNSSimulationSubScPresenter(parent=self,
                                                     view=self.view,
                                                     model=self.model)

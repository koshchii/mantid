# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
"""
DNS Simulation widget
"""
from mantidqtinterfaces.dns_simulation.simulation_sub_powder_presenter import \
    DNSSimulationSubPowderPresenter

from mantidqtinterfaces.dns_simulation.simulation_sub_powder_model import \
    DNSSimulationSubPowderModel

from mantidqtinterfaces.dns_simulation.simulation_sub_powder_view import \
    DNSSimulationSubPowderView


class DNSSimulationSubPowderWidget:
    def __init__(self, parent):
        super().__init__()
        self.name = "Powder diffraktorgram"
        self.view = DNSSimulationSubPowderView(parent=parent.view)
        self.model = DNSSimulationSubPowderModel(parent=self)
        self.presenter = DNSSimulationSubPowderPresenter(parent=self,
                                                         view=self.view,
                                                         model=self.model,
                                                         )

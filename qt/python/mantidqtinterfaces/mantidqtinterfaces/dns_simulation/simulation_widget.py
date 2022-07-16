# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
"""
DNS Simulation widget
"""
from mantidqtinterfaces.dns_powder_tof.data_structures.dns_widget import \
    DNSWidget
from mantidqtinterfaces.dns_simulation.simulation_model import \
    DNSSimulationModel
from mantidqtinterfaces.dns_simulation.simulation_presenter import \
    DNSSimulationPresenter
from mantidqtinterfaces.dns_simulation.simulation_view import DNSSimulationView
from mantidqtinterfaces.dns_simulation.simulation_sub_table_widget import \
    DNSSimulationSubTableWidget
from mantidqtinterfaces.dns_simulation.simulation_sub_powder_widget import \
    DNSSimulationSubPowderWidget
from mantidqtinterfaces.dns_simulation.simulation_sub_sc_widget import \
    DNSSimulationSubScWidget


class DNSSimulationWidget(DNSWidget):
    def __init__(self, name, parent):
        super().__init__(name, parent)
        self.sub_widgets = [DNSSimulationSubTableWidget(parent=self),
                            DNSSimulationSubPowderWidget(parent=self),
                            DNSSimulationSubScWidget(parent=self)
                            ]
        self.view = DNSSimulationView(parent=parent.view)
        self.model = DNSSimulationModel(parent=self)
        self.presenter = DNSSimulationPresenter(parent=self,
                                                view=self.view,
                                                model=self.model,
                                                name=name,
                                                subwidgets=self.sub_widgets)

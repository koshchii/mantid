# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +

"""
Definition of DNS measurement modi and functions for switching.
"""

from collections import OrderedDict

from mantidqtinterfaces.dns_powder_tof.file_selector.file_selector_widget \
    import DNSFileSelectorWidget
from mantidqtinterfaces.dns_powder_tof.options.tof_powder_options_widget import \
    DNSTofPowderOptionsWidget
from mantidqtinterfaces.dns_powder_tof.paths.path_widget import DNSPathWidget
from mantidqtinterfaces.dns_powder_tof.plot.tof_powder_plot_widget import \
    DNSTofPowderPlotWidget
from mantidqtinterfaces.dns_powder_tof.script_generator. \
    tof_powder_script_generator_widget import \
    DNSTofPowderScriptGeneratorWidget
from mantidqtinterfaces.dns_powder_tof.xml_dump.xml_dump_widget import \
    DNSXMLDumpWidget

# powder elastic
from mantidqtinterfaces.dns_powder_elastic.options.elastic_powder_options_widget import \
    DNSElasticPowderOptionsWidget
from mantidqtinterfaces.dns_powder_elastic.plot.elastic_powder_plot_widget import \
    DNSElasticPowderPlotWidget
from mantidqtinterfaces.dns_powder_elastic.script_generator.elastic_powder_script_generator_widget import \
    DNSElasticPowderScriptGeneratorWidget

# single crystal elastic
from mantidqtinterfaces.dns_sc_elastic.options.elastic_sc_options_widget import \
    DNSElasticSCOptionsWidget
from mantidqtinterfaces.dns_sc_elastic.plot.elastic_single_crystal_plot_widget import \
    DNSElasticSCPlotWidget
from mantidqtinterfaces.dns_sc_elastic.script_generator.elastic_sc_script_generator_widget import \
    DNSElasticSCScriptGeneratorWidget

# simulation
from mantidqtinterfaces.dns_simulation.simulation_widget import \
    DNSSimulationWidget


class DNSModus:
    """
    Defines the different reduction modes and which widgets are used in
    each mode.
    """

    def __init__(self, name, parent):
        super().__init__()
        self.parent = parent  # should be main widget
        self.name = name  # only names in the mapping below are allowed
        self.widgets = OrderedDict()

        self._mode_map = {
            'powder_tof': [
                'paths', 'file_selector', 'tof_powder_options',
                'tof_powder_script_generator', 'xml_dump', 'plot_tof_powder'
            ],
            'powder_elastic': [
                'paths', 'file_selector', 'elastic_powder_options',
                'elastic_powder_script_generator', 'xml_dump',
                'plot_elastic_powder'
            ],
            'sc_elastic': [
                'paths', 'file_selector', 'elastic_sc_options',
                'elastic_sc_script_generator', 'xml_dump',
                'plot_elastic_sc'
            ],
            'simulation': ['simulation'],
        }
        # Yapf: disable
        self._widget_map = {
            'paths': DNSPathWidget,
            'file_selector': DNSFileSelectorWidget,
            'tof_powder_options': DNSTofPowderOptionsWidget,
            'tof_powder_script_generator': DNSTofPowderScriptGeneratorWidget,
            'xml_dump': DNSXMLDumpWidget,
            'plot_tof_powder': DNSTofPowderPlotWidget,

            # powder elastic
            'elastic_powder_options': DNSElasticPowderOptionsWidget,
            'elastic_powder_script_generator':
                DNSElasticPowderScriptGeneratorWidget,
            'plot_elastic_powder': DNSElasticPowderPlotWidget,

            # powder single crystal
            'elastic_sc_options': DNSElasticSCOptionsWidget,
            'elastic_sc_script_generator': DNSElasticSCScriptGeneratorWidget,
            'plot_elastic_sc': DNSElasticSCPlotWidget,

            # simulation
            'simulation': DNSSimulationWidget,
        }
        # Yapf: enable

    def change(self, name):
        """
        Changes the mode of DNS GUI only names from _mode_map are
        allowed.
        """
        self.name = name
        self.widgets.clear()
        for key in self._mode_map.get(self.name, []):
            self.widgets[key] = self._widget_map[key](key, self.parent)
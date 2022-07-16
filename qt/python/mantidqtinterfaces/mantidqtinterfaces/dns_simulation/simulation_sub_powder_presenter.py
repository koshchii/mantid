# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
"""
DNS View for simulation elastic DNS data
"""
from mantidqtinterfaces.dns_simulation.simulation_sub_common_presenter import \
    DNSSimulationSubCommonPresenter


class DNSSimulationSubPowderPresenter(DNSSimulationSubCommonPresenter):
    """
        Sub Widget to show Table of Relfections for Simulation
    """

    def __init__(self, parent=None, view=None, model=None):
        super().__init__(parent, view, model)
        self.view.sig_powderplot_clicked.connect(self._powderplot)
        self._sub_dict = None

    def process_request(self, sub_dict):
        self._sub_dict = sub_dict
        self._powderplot()

    def _powderplot(self):
        if self._sub_dict is None:
            return
        refls = self._sub_dict['refls']
        own_dict = self.get_option_dict()
        start = own_dict['powder_start']
        end = own_dict['powder_end']
        shift = own_dict['shift']
        x, y = self.model.create_powder_profile(refls, start, end, shift)
        annotate_list = self.model.get_annotation_list(refls, start,
                                                       end, shift, y)
        self.view.start_powderplot(x, y)
        self._annotate_reflections(annotate_list)
        self.view.finish_powderplot()

    def _annotate_reflections(self, refl_to_annotate):
        own_dict = self.get_option_dict()
        if own_dict['labels']:
            for i, hkl in enumerate(refl_to_annotate[1]):
                label = f'  [{hkl[0]:4.2f}, {hkl[1]:4.2f}, {hkl[2]:4.2f}]'
                x = refl_to_annotate[0][i]
                y = refl_to_annotate[2][i]
                self.view.annotate_reflection(label, x, y)

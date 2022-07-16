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


class DNSSimulationSubScPresenter(DNSSimulationSubCommonPresenter):
    """
        Sub Widget to show Table of Relfections for Simulation
    """

    def __init__(self, parent=None, view=None, model=None):
        super().__init__(parent, view, model)
        self.view.sig_scplot_clicked.connect(self._sc_plot)
        self.view.sig_mouse_pos_changed.connect(self._set_hkl_pos_on_plot)

        self._sub_dict = None

    def process_request(self, sub_dict):
        self._sub_dict = sub_dict
        self._sc_plot()
        self._toggle_off_set_warning(sub_dict.get('oof_set', False))

    def _sc_plot(self):
        if self._sub_dict is None:
            return
        own_dict = self.get_option_dict()
        q1 = self._sub_dict['hkl1_v']
        q2 = self._sub_dict['hkl2_p_v']
        refls = self._sub_dict['inplane_refls']
        wavelength = self._sub_dict['wavelength']
        orilat = self._sub_dict['orilat']

        line = self.model.create_dns_surface(orilat, q1, q2, wavelength,
                                             own_dict)
        self.view.start_sc_plot(line)
        self._scatter_plot(refls)
        self.view.finish_sc_plot(
            q1, f"[{q2[0]:5.3f}, {q2[1]:5.3f}, {q2[2]:5.3f}]")

    def _scatter_plot(self, refls):
        inten_max_min = self.model.get_min_max_int(refls)
        if refls.any():
            self.view.scatter_plot(
                x=refls[:, 0],
                y=refls[:, 1],
                intensity=refls[:, 2],
                inten_max_min=inten_max_min)
            self._annotate_refl(refls)

    def _annotate_refl(self, refls):
        for i in range(len(refls[:, 0])):
            label = f'    {refls[i, 2]:.0f}\n    [{refls[i, 3]:5.2f},' \
                    f' {refls[i, 4]:5.2f}, {refls[i, 5]:5.2f}]'
            x = refls[i, 0]
            y = refls[i, 1]
            self.view.annotate_refl(label, x, y)

    def _set_hkl_pos_on_plot(self, x, y):
        if self._sub_dict is None:
            return
        hkl = self.model.get_hkl_on_plot(x, y, self._sub_dict['hkl1_v'],
                                         self._sub_dict['hkl2_p_v'])
        if hkl is not None:
            self.view.set_hkl_position_on_plot(hkl)

    def _toggle_off_set_warning(self, checked):
        if checked:
            self.view.set_off_warning('')
        else:
            self.view.set_off_warning('Warning: omega offset not set')

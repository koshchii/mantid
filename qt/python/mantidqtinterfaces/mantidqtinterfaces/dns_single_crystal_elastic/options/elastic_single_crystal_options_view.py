# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +

"""
DNS single crystal elastic options tab view of DNS reduction GUI.
"""

from qtpy.QtCore import Signal

from mantidqt.utils.qt import load_ui

from mantidqtinterfaces.dns_powder_tof.data_structures.dns_view import DNSView


class DNSElasticSCOptionsView(DNSView):
    """
    Widget that lets user select reduction options.
    """

    NAME = "Options"

    def __init__(self, parent):
        super().__init__(parent)
        self._content = load_ui(__file__, "elastic_single_crystal_options.ui", baseinstance=self)
        self._map = {
            "wavelength": self._content.dSB_wavelength,
            "get_wavelength": self._content.cB_get_wavelength,
            "norm_time": self._content.rB_norm_time,
            "norm_monitor": self._content.rB_norm_monitor,
            "corrections": self._content.gB_corrections,
            "det_efficiency": self._content.cB_det_efficiency,
            "sum_vana_sf_nsf": self._content.cB_sum_vana_sf_nsf,
            "ignore_vana_fields": self._content.cB_ignore_vana_fields,
            "flipping_ratio": self._content.cB_flipping_ratio,
            "subtract_background_from_sample": self._content.cB_subtract_background_from_sample,
            "background_factor": self._content.dSB_background_factor,
            "automatic_binning": self._content.cB_automatic_binning,
            "two_theta_min": self._content.dSB_two_theta_min,
            "two_theta_max": self._content.dSB_two_theta_max,
            "two_theta_bin_size": self._content.dSB_two_theta_bin_size,
            "omega_min": self._content.dSB_omega_min,
            "omega_max": self._content.dSB_omega_max,
            "omega_bin_size": self._content.dSB_omega_bin_size,
            "a": self._content.dSB_a,
            "b": self._content.dSB_b,
            "c": self._content.dSB_c,
            "alpha": self._content.dSB_alpha,
            "beta": self._content.dSB_beta,
            "gamma": self._content.dSB_gamma,
            "hkl1": self._content.lE_hkl1,
            "hkl2": self._content.lE_hkl2,
            "omega_offset": self._content.dSB_omega_offset,
            "use_dx_dy": self._content.cB_use_dx_dy,
            "dx": self._content.dSB_dx,
            "dy": self._content.dSB_dy,
        }

        # connect signals
        self._attach_signal_slots()

    # signals
    sig_get_wavelength = Signal()
    sig_two_theta_max_changed = Signal()
    sig_two_theta_min_changed = Signal()
    sig_omega_max_changed = Signal()
    sig_omega_min_changed = Signal()
    sig_auto_binning_clicked = Signal(int)

    def deactivate_get_wavelength(self):
        self._map["get_wavelength"].setCheckState(0)

    def _disable_det_efficiency(self, state):
        if state == 0:
            self._map["ignore_vana_fields"].setChecked(False)
            self._map["ignore_vana_fields"].setEnabled(False)
            self._map["sum_vana_sf_nsf"].setChecked(False)
            self._map["sum_vana_sf_nsf"].setEnabled(False)
        else:
            self._map["ignore_vana_fields"].setEnabled(True)
            self._map["sum_vana_sf_nsf"].setEnabled(True)

    def _disable_sum_vanadium(self, state):
        if state == 0:
            self._map["ignore_vana_fields"].setEnabled(True)
        else:
            self._map["ignore_vana_fields"].setEnabled(False)
            self._map["ignore_vana_fields"].setChecked(False)

    def _disable_ignore_vana(self, state):
        if state == 0:
            self._map["sum_vana_sf_nsf"].setEnabled(True)
        else:
            self._map["sum_vana_sf_nsf"].setEnabled(False)
            self._map["sum_vana_sf_nsf"].setChecked(False)

    def _disable_lattice(self, state):
        self._map["a"].setEnabled(not state)
        self._map["b"].setEnabled(not state)
        self._map["c"].setEnabled(not state)
        self._map["alpha"].setEnabled(not state)
        self._map["beta"].setEnabled(not state)
        self._map["gamma"].setEnabled(not state)
        self._map["dx"].setEnabled(state)
        self._map["dy"].setEnabled(state)

    def _disable_corrections(self, state):
        if state == 0:
            self._disable_det_efficiency(False)
            self._map["det_efficiency"].setChecked(False)
            self._map["flipping_ratio"].setChecked(False)
            self._disable_subtract_sample_background(False)
            self._map["subtract_background_from_sample"].setChecked(False)

    def _disable_subtract_sample_background(self, state):
        self._map["background_factor"].setEnabled(state)

    def _automatic_binning_clicked(self, state):
        self._map["two_theta_min"].setEnabled(not state)
        self._map["two_theta_max"].setEnabled(not state)
        self._map["two_theta_bin_size"].setEnabled(not state)
        self._map["omega_min"].setEnabled(not state)
        self._map["omega_max"].setEnabled(not state)
        self._map["omega_bin_size"].setEnabled(not state)
        self.sig_auto_binning_clicked.emit(state)

    def _get_wavelength(self, state):
        if state:
            self.sig_get_wavelength.emit()

    def _two_theta_min_changed(self):
        self.sig_two_theta_min_changed.emit()

    def _two_theta_max_changed(self):
        self.sig_two_theta_max_changed.emit()

    def _omega_min_changed(self):
        self.sig_omega_min_changed.emit()

    def _omega_max_changed(self):
        self.sig_omega_max_changed.emit()

    def _attach_signal_slots(self):
        self._map["wavelength"].valueChanged.connect(self.deactivate_get_wavelength)
        self._map["get_wavelength"].stateChanged.connect(self._get_wavelength)
        self._map["det_efficiency"].stateChanged.connect(self._disable_det_efficiency)
        self._map["subtract_background_from_sample"].stateChanged.connect(self._disable_subtract_sample_background)
        self._map["corrections"].clicked.connect(self._disable_corrections)
        self._map["sum_vana_sf_nsf"].stateChanged.connect(self._disable_sum_vanadium)
        self._map["ignore_vana_fields"].stateChanged.connect(self._disable_ignore_vana)
        self._map["use_dx_dy"].stateChanged.connect(self._disable_lattice)
        self._map["two_theta_min"].valueChanged.connect(self._two_theta_min_changed)
        self._map["two_theta_max"].valueChanged.connect(self._two_theta_max_changed)
        self._map["omega_min"].valueChanged.connect(self._omega_min_changed)
        self._map["omega_max"].valueChanged.connect(self._omega_max_changed)
        self._map["automatic_binning"].stateChanged.connect(self._automatic_binning_clicked)

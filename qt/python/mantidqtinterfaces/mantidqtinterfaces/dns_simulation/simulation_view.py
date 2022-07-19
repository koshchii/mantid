# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
"""
DNS View for simulation elastic DNS data
"""
from os.path import expanduser

from qtpy.QtCore import Signal
from qtpy.QtWidgets import QFileDialog, QTableWidgetItem

from mantidqt.utils.qt import load_ui

from mantidqtinterfaces.dns_powder_tof.data_structures.dns_view import DNSView


class MyTableWidgetItem(QTableWidgetItem):
    def __init__(self, number):
        QTableWidgetItem.__init__(self, number, QTableWidgetItem.UserType)
        self.__number = float(number)

    def __lt__(self, other):
        # pylint: disable=W0212
        return self.__number < other.__number


class DNSSimulationView(DNSView):
    """
        Widget that lets user select DNS data directories
    """
    NAME = "Simulation"

    def __init__(self, parent):
        super().__init__(parent)
        self._content = load_ui(__file__, 'simulation.ui', baseinstance=self)

        # input signals
        self._content.pB_load_cif.clicked.connect(self._open_filename_dialog)
        self._content.pB_clear_cif.clicked.connect(self._clearfilename)
        self._content.dSB_a.valueChanged.connect(self._unitcellchanged)
        self._content.dSB_c.valueChanged.connect(self._unitcellchanged)
        self._content.dSB_b.valueChanged.connect(self._unitcellchanged)
        self._content.dSB_alpha.valueChanged.connect(self._unitcellchanged)
        self._content.dSB_beta.valueChanged.connect(self._unitcellchanged)
        self._content.dSB_gamma.valueChanged.connect(self._unitcellchanged)
        self._content.lE_Spacegroup.textEdited.connect(self._unitcellchanged)
        self._content.cB_inplane.toggled.connect(self._inplane_unique_switched)
        self._content.cB_unique.toggled.connect(self._inplane_unique_switched)
        self._content.dSB_sample_rot.valueChanged.connect(
            self._sample_rotchanged)
        self._content.lE_hkl1.textEdited.connect(self._clear_d_tooltip)
        self._content.lE_hkl2.textEdited.connect(self._clear_d_tooltip)
        self._content.lE_hkl2_p.textEdited.connect(self._clear_d_tooltip)
        self._content.cB_fix_omega.toggled.connect(self._fixomegaoffset)
        self._content.dSB_wavelength.valueChanged.connect(
            self._wavelength_changed)
        self._content.pB.clicked.connect(self._calculate_clicked)

        self._map = {
            'hkl2': self._content.lE_hkl2,
            'cif_filename': self._content.lE_cif_filename,
            'omega_offset': self._content.dSB_omega_offset,
            'spacegroup': self._content.lE_Spacegroup,
            'wavelength': self._content.dSB_wavelength,
            'det_rot': self._content.dSB_det_rot,
            'det_number': self._content.SB_det_number,
            'fix_omega': self._content.cB_fix_omega,
            'hkl1': self._content.lE_hkl1,
            'beta': self._content.dSB_beta,
            'alpha': self._content.dSB_alpha,
            'unique': self._content.cB_unique,
            'hkl2_p': self._content.lE_hkl2_p,
            'a': self._content.dSB_a,
            'c': self._content.dSB_c,
            'b': self._content.dSB_b,
            'sample_rot': self._content.dSB_sample_rot,
            'inplane': self._content.cB_inplane,
            'gamma': self._content.dSB_gamma,
        }

    def add_tab(self, tab, name, pos=-1):
        self._content.tabs.insertTab(pos, tab, name)

    # custom signals for presenter
    sig_cif_set = Signal(str)
    sig_unitcell_changed = Signal()
    sig_inplane_unique_switched = Signal()
    sig_wavelength_changed = Signal()
    sig_calculate_clicked = Signal()
    sig_fixomegaoffset = Signal()
    sig_sample_rot_changed = Signal()

    # emitting custom signals for presenter
    def _calculate_clicked(self):
        self.sig_calculate_clicked.emit()

    def _fixomegaoffset(self, checked):
        self._map['omega_offset'].setEnabled(not checked)
        self._calculate_clicked()

    def _inplane_unique_switched(self):
        self.sig_inplane_unique_switched.emit()

    def _unitcellchanged(self):
        self.sig_unitcell_changed.emit()

    def _wavelength_changed(self):
        self.sig_wavelength_changed.emit()

    def _clear_d_tooltip(self):
        """called if hkl changed but no new calculation of d was done"""
        self._map['hkl1'].setToolTip('')
        self._map['hkl2'].setToolTip('')
        self._map['hkl2_p'].setToolTip('')

    def _clearfilename(self):
        """removes cif filename"""
        self._content.lE_cif_filename.setText('')

    def _open_filename_dialog(self):
        """open dialog to select ciffile and calls loadCif()"""
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name = QFileDialog.getOpenFileName(
            self,
            "QFileDialog.getOpenFileName()",
            expanduser("~"),
            "Cif Files (*.cif);;All Files (*)",
            options=options)
        if isinstance(
                file_name, tuple
        ):  # Qt4 vs QT5 hack QT5 returns string, filter as second argument
            file_name = file_name[0]
        if file_name:
            self._content.lE_cif_filename.setText(file_name)
            self.sig_cif_set.emit(file_name)

    def _sample_rotchanged(self, _value):
        self.sig_sample_rot_changed.emit()

    def set_d_tooltip(self, d_hkl1, d_hkl2, d_hkl2_p):
        self._map['hkl1'].setToolTip(f'd: {d_hkl1:.3f}')
        self._map['hkl2'].setToolTip(f'd: {d_hkl2:.3f}')
        self._map['hkl2_p'].setToolTip(f'd: {d_hkl2_p:.3f}')

    def set_hkl2_p(self, q2_p):
        self.set_single_state(self._map['hkl2_p'],
                              f"[{q2_p[0]:5.3f},{q2_p[1]:5.3f},"
                              f"{q2_p[2]:5.3f}]")

    def set_ki(self, ki):
        self._content.l_kival.setText(f"{ki:.3f}")

    def set_omega_offset(self, offset):
        self.set_single_state(self._map['omega_offset'], offset)

    def set_spacegroup(self, spacegroup):
        self.set_single_state(self._map['spacegroup'], spacegroup)

from qtpy.QtWidgets import QProxyStyle, QStyle

from mantidqtinterfaces.dns_single_crystal_elastic.data_structures.dns_dialog import DNSDialog


class SpinboxNorepeatStyle(QProxyStyle):  # turn off auto repeat on spinbox
    def styleHint(self, hint, option=None, widget=None, returnData=None):
        if hint == QStyle.SH_SpinBox_KeyPressAutoRepeatRate:
            return 10**5
        if hint == QStyle.SH_SpinBox_ClickAutoRepeatRate:
            return 10**5
        if hint == QStyle.SH_SpinBox_ClickAutoRepeatThreshold:
            return 10**5
        return super().styleHint(hint, option, widget, returnData)


class DNSOmegaOffsetDialog(DNSDialog):
    def __init__(self, parent, omega_offset):
        super().__init__(files=__file__, ui="/omega_offset.ui")
        self._content.dSB_omegaoffset.setStyle(SpinboxNorepeatStyle())
        self._initial_omega_offset = omega_offset
        self._content.dSB_omegaoffset.setValue(self._initial_omega_offset)
        self._content.dSB_omegaoffset.setFocus()
        self._content.dSB_omegaoffset.valueChanged.connect(self._omega_offset_changed)
        self._content.pB_discard.clicked.connect(self._discard)
        self._content.pB_apply.clicked.connect(self._apply)

        self.parent = parent

    def _omega_offset_changed(self, omega_offset):
        self.parent.sig_update_omega_offset.emit(omega_offset)

    def _apply(self):
        self.reject()

    def _discard(self):
        self._omega_offset_changed(self._initial_omega_offset)
        self.reject()


class DNSdxdyDialog(DNSDialog):
    def __init__(self, parent, dx, dy):
        super().__init__(files=__file__, ui="/dxdy.ui")
        self._initial_dx = dx
        self._initial_dy = dy
        self._content.dSB_dx.setStyle(SpinboxNorepeatStyle())
        self._content.dSB_dy.setStyle(SpinboxNorepeatStyle())
        self._content.dSB_dx.setValue(dx)
        self._content.dSB_dx.setFocus()
        self._content.dSB_dx.valueChanged.connect(self._dxdy_changed)
        self._content.dSB_dy.setValue(dy)
        self._content.dSB_dy.valueChanged.connect(self._dxdy_changed)
        self._content.pB_discard.clicked.connect(self._discard)
        self._content.pB_apply.clicked.connect(self._apply)
        self._parent = parent

    def _dxdy_changed(self, dx=None, dy=None):
        if dx is None or dy is None:
            dx = self._content.dSB_dx.value()
            dy = self._content.dSB_dy.value()
        self._parent.sig_update_dxdy.emit(dx, dy)

    def _apply(self):
        self.reject()

    def _discard(self):
        self._dxdy_changed(dx=self._initial_dx, dy=self._initial_dy)
        self.reject()

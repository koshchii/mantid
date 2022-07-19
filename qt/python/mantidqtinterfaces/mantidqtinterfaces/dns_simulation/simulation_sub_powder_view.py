# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
"""
DNS View for simulation elastic DNS data
"""
from mantidqt.utils.qt import load_ui

from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.backend_qt5agg import \
    NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

from qtpy import QtWidgets
from qtpy.QtCore import Signal

from mantidqtinterfaces.dns_powder_tof.data_structures.dns_view import DNSView


class DNSSimulationSubPowderView(DNSView):
    """
        Sub widget for simulation powder plots
    """

    def __init__(self, parent):
        super().__init__(parent)
        content = load_ui(__file__, 'simulation_powder_widget.ui',
                          baseinstance=self)
        self._ax = None
        self._map = {
            'powder_start': content.dSB_powder_start,
            'shift': content.dSB_shift,
            'powder_end': content.dSB_powder_end,
            'labels': content.cB_labels,

        }
        powd_layout = content.powd_plot_layout
        self._powd_static_canvas = FigureCanvas(
            Figure(figsize=(5, 3), dpi=200))
        self._powd_static_canvas.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                               QtWidgets.QSizePolicy.Expanding)
        self._powd_toolbar = NavigationToolbar(self._powd_static_canvas, self)
        powd_layout.addWidget(self._powd_toolbar)
        powd_layout.addWidget(self._powd_static_canvas)

        # Signals
        self._map['labels'].toggled.connect(self._powderplot_clicked)
        content.pB_ps.clicked.connect(self._powderplot_clicked)

    # custom signals for presenter
    sig_powderplot_clicked = Signal()

    def _powderplot_clicked(self):
        self.sig_powderplot_clicked.emit()

    def annotate_reflection(self, label, x, y):
        self._ax.annotate(label, (x, y), fontsize=10)

    def start_powderplot(self, x, y):
        self._powd_static_canvas.figure.clear()
        self._ax = self._powd_static_canvas.figure.subplots()
        self._ax.plot(x, y, zorder=1)
        self._ax.set_xlabel('2 theta', fontsize=14)
        self._ax.set_ylabel("Intensity (M*F2)", fontsize=14)

    def finish_powderplot(self):
        self._powd_toolbar.update()
        self._ax.figure.canvas.draw()

from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.backend_qt5agg import \
    NavigationToolbar2QT as NavigationToolbar
from matplotlib.cm import get_cmap
from matplotlib.figure import Figure

from qtpy import QtWidgets
from qtpy.QtCore import Signal

from mantidqt.utils.qt import load_ui

from mantidqtinterfaces.dns_powder_tof.data_structures.dns_view import DNSView


class DNSSimulationSubScView(DNSView):
    """
        Sub widget to plot sc simulation data
    """

    def __init__(self, parent):
        super().__init__(parent)
        content = load_ui(__file__, 'simulation_sc_widget.ui',
                          baseinstance=self)
        self._cid = None
        self._ax = None
        self._content = content
        self._map = {
            'sc_det_start': content.dSB_sc_det_start,
            'sc_sam_end': content.dSB_sc_sam_end,
            'sc_det_end': content.dSB_sc_det_end,
            'sc_sam_start': content.dSB_sc_sam_start,
            'sc_oofset_warning': content.l_warning_ooset,
            'sc_show_hkl': content.l_show_hkl,
            'plot_sc': content.pB_sc,
        }
        self._map['plot_sc'].clicked.connect(self._scplot_clicked)

        # m atplotlib layout
        sc_layout = self._content.sc_plot_layout
        self._sc_static_canvas = FigureCanvas(Figure(figsize=(5, 3), dpi=200))
        self._sc_static_canvas.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                             QtWidgets.QSizePolicy.Expanding)
        self._sc_toolbar = NavigationToolbar(self._sc_static_canvas, self)
        sc_layout.addWidget(self._sc_toolbar)
        sc_layout.addWidget(self._sc_static_canvas)

        self._sc_static_canvas.mpl_connect('axes_enter_event',
                                           self._mouseonplot)
        self._sc_static_canvas.mpl_connect('axes_leave_event',
                                           self._mouseoutplot)

    sig_mouse_pos_changed = Signal(float, float)
    sig_scplot_clicked = Signal()

    def set_off_warning(self, text):
        self._content.l_warning_ooset.setText(text)

    def _mouse_pos_changed(self, event):
        x = event.xdata
        y = event.ydata
        if x is not None and y is not None:
            self.sig_mouse_pos_changed.emit(x, y)

    def _scplot_clicked(self):
        self.sig_scplot_clicked.emit()

    def _mouseonplot(self, _event):
        self._cid = self._sc_static_canvas.mpl_connect('motion_notify_event',
                                                       self._mouse_pos_changed)

    def _mouseoutplot(self, _event):
        self._sc_static_canvas.mpl_disconnect(self._cid)

    def annotate_refl(self, label, x, y):
        self._ax.annotate(label, (x, y), fontsize=10, zorder=200)

    def scatter_plot(self, x, y, intensity, inten_max_min):
        cm = get_cmap('plasma')
        sc = self._ax.scatter(x, y, c=intensity,
                              vmin=inten_max_min[1],
                              vmax=inten_max_min[0],
                              s=300,
                              cmap=cm,
                              zorder=20)
        cb = self._sc_static_canvas.figure.colorbar(sc)
        cb.set_label('Intensity', fontsize=14)
        cb.ax.zorder = -1

    def start_sc_plot(self, line):
        self._sc_static_canvas.figure.clear()
        self._ax = self._sc_static_canvas.figure.subplots()
        self._ax.fill(line[:, 0], line[:, 1], zorder=1)

    def finish_sc_plot(self, xlabel, ylabel):
        self._ax.set_xlabel(xlabel, fontsize=14)
        self._ax.set_ylabel(ylabel, fontsize=14)
        self._sc_toolbar.update()
        self._ax.grid(color='grey', linestyle=':', linewidth=1)
        self._ax.figure.canvas.draw()

    def set_hkl_position_on_plot(self, hkl):
        self._content.l_show_hkl.setText(
            f'hkl = [{hkl[0]:5.2f} , {hkl[1]:5.2f} , {hkl[2]:5.2f}]')

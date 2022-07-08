# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +

"""
DNS View for simulation elastic DNS data.
"""

from mantidqt.utils.qt import load_ui

from qtpy import QtCore, QtGui, QtWidgets
from qtpy.QtCore import Signal
from qtpy.QtWidgets import QTableWidgetItem

from mantidqtinterfaces.dns_powder_tof.data_structures.dns_view import DNSView


class MyTableWidgetItem(QTableWidgetItem):
    def __init__(self, number):
        QTableWidgetItem.__init__(self, number, QTableWidgetItem.UserType)
        self.__number = float(number)

    def __lt__(self, other):
        # pylint: disable=W0212
        return self.__number < other.__number


class DNSSimulationSubTableView(DNSView):
    """
    Sub Widget to show Table of Reflections for Simulation
    """
    def __init__(self, parent):
        super().__init__(parent)
        self._content = load_ui(__file__, 'simulation_table_widget.ui',
                                baseinstance=self)

        self._content.table.setEditTriggers(
            QtWidgets.QTableWidget.NoEditTriggers)
        self._content.table.sortByColumn(3, 0)
        # sort after q unlesss user changes it
        self._map = {
        }
        self._cellinfo = None

    sig_table_item_clicked = Signal(float, float)

    def _get_det_rot(self, item):
        return float(self._content.table.item(item.row(),
                                              item.column() - 2).text())

    def _tableitemdclicked(self, item):
        if item.column() != 11:
            return
        dr = self._get_det_rot(item)
        sr = float(item.text())
        self.sig_table_item_clicked.emit(dr, sr)

    def set_mult_tooltip(self, tooltip):
        self._cellinfo.setToolTip(tooltip)

    def set_bg_color(self, identified):
        if identified:
            self._cellinfo.setBackground(QtGui.QColor(100, 200, 100))
        else:
            self._cellinfo.setBackground(QtGui.QColor(255, 255, 255))

    def create_tableitem(self, reflstr):
        self._cellinfo = MyTableWidgetItem(reflstr)
        # explicit conversion to int to fix deprecated bug, in current
        # python version used by mantid
        self._cellinfo.setTextAlignment(int(QtCore.Qt.AlignRight
                                            | QtCore.Qt.AlignVCenter))

    def set_tableitem(self, row, col):
        self._content.table.setItem(row, col, self._cellinfo)

    def start_table(self, rowcount, columncount):
        self._content.table.clearContents()
        self._content.table.setSortingEnabled(False)
        self._content.table.setColumnCount(columncount)
        self._content.table.setRowCount(rowcount)

    def finish_table(self):
        self._content.table.resizeColumnsToContents()
        self._content.table.itemDoubleClicked.connect(self._tableitemdclicked)
        self._content.table.setSortingEnabled(True)

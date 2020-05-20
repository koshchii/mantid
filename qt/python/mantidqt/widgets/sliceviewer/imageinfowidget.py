# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2020 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +
#  This file is part of the mantidqt package
#
#
from mantidqt.utils.qt import import_qt

from qtpy.QtWidgets import QVBoxLayout, QCheckBox, QWidget
from qtpy.QtCore import Qt

ImageInfoWidget_cpp = import_qt('.._common', 'mantidqt.widgets', 'ImageInfoWidget')


class ImageInfoWidget(QWidget):

    def __init__(self, workspace, parent):
        super(QWidget, self).__init__(parent)

        layout = QVBoxLayout()
        self.workspace = workspace
        self.track_cursor = QCheckBox("Track Cursor", self)
        self.track_cursor.setChecked(True)
        self.table_widget = ImageInfoWidget_cpp(workspace, self)
        layout.addWidget(self.track_cursor, 0, Qt.AlignRight)
        layout.addWidget(self.table_widget)
        self.setLayout(layout)

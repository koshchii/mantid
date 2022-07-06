# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +

"""
DNS View for simulation elastic DNS data.
"""

from mantidqtinterfaces.dns_simulation.simulation_sub_common_presenter import \
    DNSSimulationSubCommonPresenter

TABLE_HEAD = [
    'h', 'k', 'l', 'q', 'd', 'two_theta', 'fs', 'mult', 'diff', 'det_rot',
    'channel', 'sample_rot'
]
FORMAT_STRING = [
    ' {0:.0f} ', ' {0:.0f} ', ' {0:.0f} ', ' {0:.2f} ', ' {0:.2f} ',
    ' {0:.2f} ', ' {0:.0f} ', ' {0:.0f} ', ' {0:.2f} ', ' {0:.2f} ',
    ' {0:.0f} ', ' {0:.2f} '
]


class DNSSimulationSubTablePresenter(DNSSimulationSubCommonPresenter):
    """
    Sub Widget to show Table of Reflections for Simulation
    """
    def __init__(self, parent=None, view=None, model=None):
        super().__init__(parent, view, model)
        self.parent = parent
        self.view = view
        self.model = model
        self.view.sig_table_item_clicked.connect(self._table_item_double_clicked)
        self._sub_dict = None

    def process_request(self, sub_dict):
        self._sub_dict = sub_dict
        filtered_refls = sub_dict['filtered_refls']
        two_theta_limit = sub_dict['two_theta_limit']
        self._write_table(filtered_refls, two_theta_limit)

    def _table_item_double_clicked(self, det_rot, sample_rot):
        """
        Sets the omega offset based on identified reflection.
        """
        self.parent.parent.presenter.back_call_from_table_item_clicked(
            det_rot, sample_rot)

    def _write_table(self, refls, two_theta_limit):  #
        """
        Writes a list of reflections to the table.
        """
        self.view.start_table(len(refls), len(TABLE_HEAD))
        row = 0
        for refl in refls:
            for col, head in enumerate(TABLE_HEAD):
                refl_str = FORMAT_STRING[col].format(getattr(refl, head))
                self.view.create_tableitem(refl_str)
                self._add_mult_tooltip(refl, col)
                self._color_identified(refl, two_theta_limit)
                self.view.set_tableitem(row, col)
            row += 1
        self.view.finish_table()

    def _color_identified(self, refl, two_theta_limit):
        self.view.set_bg_color(refl.diff < two_theta_limit)

    def _add_mult_tooltip(self, refl, col):
        if TABLE_HEAD[col] == 'mult':
            self.view.set_mult_tooltip(str(refl.equivalents))

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

TABLEHEAD = [
    'h', 'k', 'l', 'q', 'd', 'tth', 'fs', 'mult', 'diff', 'det_rot',
    'channel', 'sample_rot'
]
FORMATSTRING = [
    ' {0:.0f} ', ' {0:.0f} ', ' {0:.0f} ', ' {0:.2f} ', ' {0:.2f} ',
    ' {0:.2f} ', ' {0:.0f} ', ' {0:.0f} ', ' {0:.2f} ', ' {0:.2f} ',
    ' {0:.0f} ', ' {0:.2f} '
]


class DNSSimulationSubTablePresenter(DNSSimulationSubCommonPresenter):
    """
        Sub Widget to show Table of Relfections for Simulation
    """

    def __init__(self, parent=None, view=None, model=None):
        super().__init__(parent, view, model)
        self.parent = parent
        self.view = view
        self.model = model
        self.view.sig_table_item_clicked.connect(self._tableitemdclicked)
        self._sub_dict = None

    def process_request(self, sub_dict):
        self._sub_dict = sub_dict
        filtered_refls = sub_dict['filtered_refls']
        tth_limit = sub_dict['tth_limit']
        self._writetable(filtered_refls, tth_limit)

    def _tableitemdclicked(self, det_rot, sample_rot):
        """ sets the omega offset based on identified reflection """
        self.parent.parent.presenter.back_call_from_tableitem_clicked(
            det_rot, sample_rot)

    def _writetable(self, refls, tthlimit):  #
        """writes a list of reflections to the table"""
        self.view.start_table(len(refls), len(TABLEHEAD))
        row = 0
        for refl in refls:
            for col, head in enumerate(TABLEHEAD):
                reflstr = FORMATSTRING[col].format(getattr(refl, head))
                self.view.create_tableitem(reflstr)
                self._add_mult_tooltip(refl, col)
                self._color_identified(refl, tthlimit)
                self.view.set_tableitem(row, col)
            row += 1
        self.view.finish_table()

    def _color_identified(self, refl, tthlimit):
        self.view.set_bg_color(refl.diff < tthlimit)

    def _add_mult_tooltip(self, refl, col):
        if TABLEHEAD[col] == 'mult':
            self.view.set_mult_tooltip(str(refl.equivalents))

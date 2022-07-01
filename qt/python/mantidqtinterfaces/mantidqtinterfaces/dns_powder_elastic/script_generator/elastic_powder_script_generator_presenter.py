# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
"""
DNS script generator for elastic powder data
"""

from mantidqtinterfaces.dns_powder_tof.script_generator. \
    common_script_generator_presenter import DNSScriptGeneratorPresenter


class DNSElasticPowderScriptGeneratorPresenter(DNSScriptGeneratorPresenter):
    def __init__(self, name=None, parent=None, view=None, model=None):
        super().__init__(parent=parent, name=name, view=view, model=model)
        self._plotlist = []

    def get_option_dict(self):
        if self.view is not None:
            self.own_dict.update(self.view.get_state())
        self.own_dict['script_path'] = self._script_path
        self.own_dict['script_number'] = self._script_number
        self.own_dict['script_text'] = self._script_text
        self.own_dict['plotlist'] = self._plotlist
        return self.own_dict

    def _finish_script_run(self):
        self._plotlist = self.model.get_plotlist()

# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +

"""
DNS View for simulation elastic DNS data.
"""


class DNSSimulationSubCommonPresenter:
    """
    Sub Widget to show Table of Reflections for Simulation.
    """
    def __init__(self, parent=None, view=None, model=None):
        # pylint: disable=unused-argument
        super().__init__()
        self.view = view
        self.model = model

    def process_request(self, sub_dict):
        pass

    def get_option_dict(self):
        """
        Return own options from view.
        """
        if self.view is not None:
            return self.view.get_state()
        return {}

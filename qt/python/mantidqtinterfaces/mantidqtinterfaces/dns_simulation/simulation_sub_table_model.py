# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +

"""
DNS View for simulation elastic DNS data.
"""


class DNSSimulationSubTableModel:
    """
    Sub Widget to show Table of Reflections for Simulation.
    """
    def __init__(self, parent):
        # pylint: disable=unused-argument
        super().__init__()

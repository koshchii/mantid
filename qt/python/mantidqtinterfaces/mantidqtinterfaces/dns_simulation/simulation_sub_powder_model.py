# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +

"""
DNS View for simulation elastic DNS data.
"""

import mantidqtinterfaces.dns_simulation.simulation_helpers as sim_help

from mantid.simpleapi import CreateWorkspace


class DNSSimulationSubPowderModel:
    """
    Sub Widget to show Table of Reflections for Simulation.
    """
    def __init__(self, parent=None):
        # pylint: disable=unused-argument
        super().__init__()

    @staticmethod
    def create_powder_profile(refls, start, end, shift):
        two_theta = sim_help.get_two_theta_range(start, end, shift)
        intensity = two_theta * 0
        for refl in refls:
            intensity += sim_help.get_intensity_prof(refl, shift, two_theta)
        CreateWorkspace(OutputWorkspace='mat_simulation',
                        DataX=sim_help.get_two_theta_bins(two_theta),
                        DataY=intensity,
                        NSpec=1,
                        UnitX='Degrees')
        return [two_theta, intensity]

    @staticmethod
    def get_annotation_list(refls, start, end, shift, intensity):
        two_theta_end = sim_help.get_two_theta_end(end, shift)
        start = sim_help.get_two_theta_start(start, shift)
        refls = sim_help.get_unique_refl(refls)
        annotate_list = [[], [], []]
        for refl in refls:
            if (refl.tth + shift <= two_theta_end
                    and round(refl.tth + shift, 2) not in annotate_list[0]):
                x_numb = int((refl.tth - start + shift) / 0.1)
                annotate_list[0].append(round(refl.tth + shift, 2))
                annotate_list[1].append(refl.hkl)
                annotate_list[2].append(intensity[x_numb])
        return annotate_list

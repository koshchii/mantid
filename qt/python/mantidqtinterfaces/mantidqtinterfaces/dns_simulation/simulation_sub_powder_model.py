# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
"""
DNS View for simulation elastic DNS data
"""
import mantidqtinterfaces.dns_simulation.simulation_helpers as sim_help

from mantid.simpleapi import CreateWorkspace


class DNSSimulationSubPowderModel:
    """
        Sub Widget to show Table of Relfections for Simulation
    """

    def __init__(self, parent=None):
        # pylint: disable=unused-argument
        super().__init__()

    @staticmethod
    def create_powder_profile(refls, start, end, shift):
        tth = sim_help.get_tth_range(start, end, shift)
        intensity = tth * 0
        for refl in refls:
            intensity += sim_help.get_intensity_prof(refl, shift, tth)
        CreateWorkspace(OutputWorkspace='mat_simulation',
                        DataX=sim_help.get_tth_bins(tth),
                        DataY=intensity,
                        NSpec=1,
                        UnitX='Degrees')
        return [tth, intensity]

    @staticmethod
    def get_annotation_list(refls, start, end, shift, intensity):
        tth_end = sim_help.get_tth_end(end, shift)
        start = sim_help.get_tth_start(start, shift)
        refls = sim_help.get_unique_refl(refls)
        annotate_list = [[], [], []]
        for refl in refls:
            if (refl.tth + shift <= tth_end
                    and round(refl.tth + shift, 2) not in annotate_list[0]):
                xnumb = int((refl.tth - start + shift) / 0.1)
                annotate_list[0].append(round(refl.tth + shift, 2))
                annotate_list[1].append(refl.hkl)
                annotate_list[2].append(intensity[xnumb])
        return annotate_list

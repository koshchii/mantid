# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +

"""
DNS widget for simulation elastic DNS data.
"""

import numpy as np
import mantidqtinterfaces.dns_simulation.simulation_helpers as sim_help


class DNSSimulationSubScModel:
    """
    Sub Widget to show Table of Reflections for Simulation.
    """
    def __init__(self, parent):
        # pylint: disable=unused-argument
        super().__init__()

    @staticmethod
    def create_dns_surface(oriented_lattice, q1, q2, wavelength, options):
        return sim_help.return_dns_surface_shape(oriented_lattice, q1, q2, wavelength,
                                                 options)

    @staticmethod
    def get_min_max_int(refls):
        if refls.any():
            return [refls[:, 2].min(), refls[:, 2].max()]
        return [0, 1]

    @staticmethod
    def get_hkl_on_plot(x, y, hkl1, hkl2):
        hkl1 = np.asarray(hkl1)
        hkl2_p = np.asarray(hkl2)
        if x is not None and y is not None:
            return hkl1 * x + hkl2_p * y
        return None

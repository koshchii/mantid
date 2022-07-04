# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
"""
Model for DNS simulation
"""

import mantidqtinterfaces.dns_simulation.simulation_helpers as sim_help
import numpy as np
from mantidqtinterfaces.dns_powder_tof.data_structures.dns_obs_model import \
    DNSObsModel
from mantidqtinterfaces.dns_powder_tof.data_structures.object_dict import \
    ObjectDict
from mantidqtinterfaces.dns_powder_elastic.helpers.converters import \
    convert_hkl_string_to_float
from mantid.geometry import (CrystalStructure, OrientedLattice,
                             ReflectionConditionFilter, ReflectionGenerator,
                             SpaceGroupFactory, UnitCell)
from mantid.simpleapi import CreateWorkspace, LoadCIF


class DNSSimulationModel(DNSObsModel):
    """
    Model for DNS simulation

    generates HKL list from CIF files or given lattice parameters and
    creates powder diffractogram and single crystal diffraction map
    """

    def __init__(self, parent):
        super().__init__(parent)
        self._sim_ws = None
        self._orilat = None
        self._refls = None
        self._generator = None
        self._cryst = None
        self._non_rot_lat = None

    @staticmethod
    def filter_refls(refls, inplane, unique):
        if inplane and not unique:
            return sim_help.get_inplane_refl(refls)
        if not inplane and unique:
            return sim_help.get_unique_refl(refls)
        if inplane and unique:
            return sim_help.get_unique_inplane_refl(refls)
        return refls

    def get_ds(self, hkl1_v, hkl2_v, hkl2_p_v):
        return self._generator.getDValues([
            sim_help.list_to_v3d(hkl1_v),
            sim_help.list_to_v3d(hkl2_v),
            sim_help.list_to_v3d(hkl2_p_v)
        ])

    def get_hkl2_p(self):
        q2_p = self._non_rot_lat.getvVector()
        q2_p = q2_p / np.linalg.norm(q2_p)
        return q2_p

    def _set_lattice(self, q1, q2):
        self._orilat.setUFromVectors(sim_help.list_to_v3d(q1),
                                     sim_help.list_to_v3d(q2))
        self._non_rot_lat.setUFromVectors(sim_help.list_to_v3d(q1),
                                          sim_help.list_to_v3d(q2))

    @staticmethod
    def _get_refl_filter(cif_filename):
        if not cif_filename:
            return ReflectionConditionFilter.SpaceGroup
        return ReflectionConditionFilter.StructureFactor

    def _get_filtered_hkls(self, wavelength, cif_filename):
        reffilter = self._get_refl_filter(cif_filename)
        maxq = sim_help.max_q(wavelength)
        hkls_unique = self._generator.getUniqueHKLsUsingFilter(
            1 / maxq, 100, reffilter)
        hkls = self._generator.getHKLsUsingFilter(1 / maxq, 100, reffilter)
        return [hkls, hkls_unique]

    def _get_q_d_val(self, hkls):
        dval = self._generator.getDValues(hkls)
        qval = [2 * np.pi / d for d in dval]
        return [dval, qval]

    def _get_ub(self, fix, omega_offset):
        ub = self._orilat.getUB()
        if fix:
            ub = sim_help.rotate_ub(omega_offset, ub)
            self._orilat.setUB(ub)
        return ub

    def get_refls_and_set_orientation(self, options):
        q1 = options['hkl1_v']
        q2 = options['hkl2_v']
        wvl = options['wavelength']
        cif_fn = options['cif_filename']
        det_rot = options['det_rot']
        det_nb = options['det_number']
        fix_omega = options['fix_omega']
        omega_offset = options['omega_offset']

        self._setcellfromparameters(options)
        self._set_lattice(q1, q2)
        self._generator = ReflectionGenerator(self._cryst)
        hkls, hkls_unique = self._get_filtered_hkls(wvl, cif_fn)
        dval, qval = self._get_q_d_val(hkls)
        f_squared = self._generator.getFsSquared(hkls)
        pg = self._cryst.getSpaceGroup().getPointGroup()
        ub = self._get_ub(fix_omega, omega_offset)
        identify_tth = sim_help.det_rot_nmb_to_tth(det_rot, det_nb)
        self._refls = []
        for i in range(len(hkls) - 1, 0, -1):
            # reverse order to get positive hkl
            refl = ObjectDict()
            refl.hkl = hkls[i]
            refl.unique = refl.hkl in hkls_unique
            refl.q = qval[i]
            refl.fs = f_squared[i]
            refl.d = dval[i]
            refl.tth = sim_help.q_to_tth(refl.q, wvl)
            refl.fs_lc = sim_help.lorentz_correction(refl.fs, refl.tth)
            if not cif_fn:
                # if dummy scatterer set I = 1
                refl.fs = 1
                refl.fs_lc = 1
            refl.h = refl.hkl[0]
            refl.k = refl.hkl[1]
            refl.l = refl.hkl[2]  # noqa: E741, E743
            refl.equivalents = pg.getEquivalents(refl.hkl)
            refl.mult = len(refl.equivalents)
            refl.diff = abs(identify_tth - refl.tth)
            refl.det_rot, refl.channel = sim_help.tth_to_rot_nmb(refl.tth)
            refl.inplane = sim_help.check_inplane(q1, q2, refl.hkl)

            refl.channel, refl.det_rot = sim_help.shift_channels_below_23(
                refl.channel, refl.det_rot)
            refl.omega = sim_help.hkl_to_omega(refl.hkl, ub, wvl,
                                               refl.tth)
            refl.sample_rot = sim_help.omega_to_samp_rot(
                refl.omega, refl.det_rot)
            self._refls.append(refl)
        return self._refls

    @staticmethod
    def _create_sim_ws():
        return CreateWorkspace(OutputWorkspace='__sim_ws',
                               DataX=[0],
                               DataY=[0],
                               NSpec=1,
                               UnitX='Degrees')

    def get_orilat(self):
        return self._orilat

    def load_cif(self, file_name):
        """uses mantid to load CIF and set crystalstructure and oriented
        lattice, calls set unitcell"""
        self._sim_ws = self._create_sim_ws()
        LoadCIF(self._sim_ws, file_name)
        self._cryst = self._sim_ws.sample().getCrystalStructure()
        unitcell = self._cryst.getUnitCell()
        self._orilat = OrientedLattice(unitcell)
        self._non_rot_lat = OrientedLattice(unitcell)
        load_dict = {
            'a': unitcell.a(),
            'b': unitcell.b(),
            'c': unitcell.c(),
            'alpha': unitcell.alpha(),
            'beta': unitcell.beta(),
            'gamma': unitcell.gamma(),
            'spacegroup': self._cryst.getSpaceGroup().getHMSymbol()
        }
        return load_dict

    def return_reflections_in_map(self, q1, q2, refls):
        refls = sim_help.return_qx_qx_inplane_refl(self._orilat, q1, q2, refls)
        return refls

    @staticmethod
    def _get_hm_spacegroup(spacegroup):
        if spacegroup.isdigit():  # if user gives SG number convert to HM
            spacegroup = SpaceGroupFactory.subscribedSpaceGroupSymbols(
                int(spacegroup))
            spacegroup = spacegroup[0]  # HM is not unique for number
        return spacegroup

    def _get_atom_str(self, cif_filename):
        if not cif_filename:
            return "Si 0 0 1 1.0 0.01"  # dummy
        return ';'.join(self._cryst.getScatterers())

    @staticmethod
    def _get_cell_str(unitcell):
        return f"{unitcell.a()} {unitcell.b()} {unitcell.c()} " \
               f"{unitcell.alpha()} {unitcell.beta()} {unitcell.gamma()}"

    def _setcellfromparameters(self, options):
        spg = options['spacegroup']
        a = options['a']
        b = options['b']
        c = options['c']
        alpha = options['alpha']
        beta = options['beta']
        gamma = options['gamma']
        if options['cifset']:
            return
        spacegroup = self._get_hm_spacegroup(spg)
        unitcell = UnitCell(a, b, c, alpha, beta, gamma, Unit=0)
        self._orilat = OrientedLattice(unitcell)
        self._non_rot_lat = OrientedLattice(unitcell)
        cellstr = self._get_cell_str(unitcell)
        atomstr = self._get_atom_str(options['cif_filename'])
        self._cryst = CrystalStructure(cellstr, spacegroup, atomstr)

    @staticmethod
    def get_ki(wavelength):
        return sim_help.ki_from_wavelength(wavelength)

    @staticmethod
    def validate_hkl(hkl1, hkl2):
        if hkl1[0] is None or hkl2[0] is None:
            error = 'Could not parse hkl, enter 3 comma seperated numbers.'
            return [False, error]
        if (np.cross(hkl1, hkl2) == 0).all():
            error = 'hkl1 cannot be paralell hkl2'
            return [False, error]
        return [True, '']

    @staticmethod
    def get_hkl_vector_dict(hkl1, hkl2):
        return {
            'hkl1_v': convert_hkl_string_to_float(hkl1),
            'hkl2_v': convert_hkl_string_to_float(hkl2)
        }

    @staticmethod
    def get_oof_from_ident(det_rot, sample_rot, id_sr, id_dr):
        return (sample_rot - det_rot) - (id_sr - id_dr)

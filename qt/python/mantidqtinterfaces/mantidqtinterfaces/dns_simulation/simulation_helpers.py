# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,

#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
"""
Helper functions only for DNS simulation calculations
"""

from numpy import arccos, arcsin, cos, log, pi, rad2deg, radians, sin, sqrt, \
    tan

import numpy as np
from mantid.kernel import V3D


def tth_to_q(tth, wavelength):
    return np.pi * 4 * np.sin(tth / 2.0) / wavelength


def tth_to_d(tth, wavelength):
    return wavelength / (2 * np.sin(tth / 2.0))


def ang_to_hkl(twotheta, omega, ub_inv, wavelength):
    """Return hkl for diffractometer angles and given inverse UB """
    theta = radians(twotheta) / 2.0
    omega = radians(omega) - theta
    uphi = np.array([-cos(omega), 0, -sin(omega)])
    hphi = 2 * sin(theta) / wavelength * uphi
    hkl = np.dot(ub_inv, hphi)
    return hkl


def check_inplane(q1, q2, hkl):
    """Return True if hkl is in the plane specified by q1 and q2"""
    q1q2 = np.column_stack((q1, q2))
    return np.linalg.det(np.column_stack((q1q2, hkl))) == 0


def det_rot_nmb_to_tth(det_rot, det_number):
    """Return 2theta for det_rot """
    return -det_rot + det_number * 5


def hkl_to_omega(hkl, ub, wavelength, tth):
    """Return omega for given hkl """
    hphi = np.transpose(np.dot(hkl, np.transpose(ub)))
    uphi = hphi * wavelength / 2.0 / sin(radians(tth / 2.0))
    sign = np.sign(arcsin(round(float(uphi[2]), 10)))
    omega = -sign * rad2deg(arccos((round(float(uphi[0]), 10)))) + tth / 2.0
    # rounding necesarry to prevent floats slightly larger 1 which will
    # kill arccos
    return omega


def ki_from_wavelength(wavelength):
    """sets ki for given wavelength"""
    return 2 * pi / wavelength


def list_to_v3d(liste):
    """make a mantid vector from a list"""
    v3dv = V3D(liste[0], liste[1], liste[2])
    return v3dv


def lorentz_correction(intensity, tth):
    """apply lorentz correction for simulation"""
    return intensity / sin(radians(tth))


def max_q(wavelength):
    """ Return maximum q at DNS for PA detector """
    maxq = 2 * sin(radians(135.5 / 2.0)) / wavelength
    return maxq


def omega_to_samp_rot(omega, det_rot):
    """Convert Omega to sample_rot """
    return omega + det_rot


def q_to_tth(q, wavelength):
    """Return two theta for given q"""
    return rad2deg(2 * arcsin(q * wavelength / (4 * pi)))


def rotate_ub(omegaoffset, ubm):
    """Rotate the UB matrix by given angle around z-axis """
    w1 = radians(omegaoffset)
    rm = np.asarray([[cos(w1), 0, -sin(w1)], [0, 1, 0], [sin(w1), 0, cos(w1)]])
    ubm = np.dot(rm, ubm)
    return ubm


def get_tth_path(tth_rang, omega_rang):
    """returns tth list of the surface boundary of dns map"""
    return np.concatenate(
        (tth_rang[0] * np.ones(len(omega_rang)), tth_rang,
         tth_rang[-1] * np.ones(len(omega_rang)), np.flip(tth_rang)))


def get_omega_path(tth_rang, omega_rang):
    """returns omega list of the surface boundary of dns map"""
    return np.concatenate(
        (np.flip(omega_rang), omega_rang[0] * np.ones(len(tth_rang)),
         omega_rang, omega_rang[-1] * np.ones(len(tth_rang))))


def return_dns_surface_shape(orilat, q1, q2, wavelength, options):
    """returns qx qy array of the surface boundary of dns sc measurement"""
    tth_start = options['sc_det_start']
    tth_end = options['sc_det_end']
    o_start = options['sc_sam_start']
    o_end = options['sc_sam_end']
    tth_rang = get_tth_range_sc(tth_start, tth_end)
    omega_rang = get_omega_range_sc(o_start, o_end, tth_start)
    ubm_inv = np.linalg.inv(orilat.getUB())
    line = []
    # create shape of dns measured surface in tth and omega
    tthr = get_tth_path(tth_rang, omega_rang)
    omegar = get_omega_path(tth_rang, omega_rang)
    for i, omega in enumerate(omegar):
        hkl = ang_to_hkl(tthr[i], omega, ubm_inv, wavelength)
        qx, qy = return_qxqy(orilat, q1, q2, hkl)
        line.append([qx, qy])
    return np.asarray(line)


def get_angle_q_vs_hkl(hkl, q, orilat):
    return radians(orilat.recAngle(hkl[0], hkl[1], hkl[2], q[0], q[1], q[2]))


def return_qxqy(orilat, q1, q2, hkl):
    """
    returns projection of hkl along q1 and q2
    """
    angle_q1_hkl = get_angle_q_vs_hkl(hkl, q1, orilat)
    angle_q2_hkl = get_angle_q_vs_hkl(hkl, q2, orilat)
    n_q1 = np.linalg.norm(orilat.qFromHKL(q1))
    n_q2 = np.linalg.norm(orilat.qFromHKL(q2))
    n_hkl = np.linalg.norm(orilat.qFromHKL(hkl))
    qx = cos(angle_q1_hkl) * n_hkl / n_q1
    qy = cos(angle_q2_hkl) * n_hkl / n_q2
    return [qx, qy]


def return_qx_qx_inplane_refl(orilat, q1, q2, refls):
    """ returns a qx,qx,Intensity,h,k,l list for reflections in plane """
    reflections = []
    for refl in refls:
        if refl.inplane:
            qx, qy = return_qxqy(orilat, q1, q2, refl.hkl)
            intensity = refl.fs_lc
            reflections.append([qx, qy, intensity, refl.h, refl.k, refl.l])
    return np.asarray(reflections)


def shift_channels_below_23(channel, det_rot):
    """
    shifting detector_rot to higher angles
    if reflection is not on last detector
    """
    while channel > 23:
        channel += -1
        det_rot += -5
    return [channel, det_rot]


def tth_to_rot_nmb(tth):
    """return det_rot and detector number for given two theta """
    det_rot = -5 - tth % 5
    channel = tth // 5 - 1
    return [det_rot, channel]


def get_omega_range_sc(o_start, o_end, tth_start):
    omega_start = o_start - tth_start
    omega_end = o_end - tth_start
    return get_1deg_angle_linspace(omega_start, omega_end)


def get_tth_end(end, shift):
    return -end + 23 * 5 + shift


def get_tth_start(start, shift):
    return -start + shift


def get_1deg_angle_linspace(start, end):
    return np.linspace(start, end, int(abs((end - start))) + 1)


def get_tth_range_sc(start, end):
    tth_start = get_tth_start(start, 0)
    tth_end = get_tth_end(end, 0)
    return get_1deg_angle_linspace(tth_start, tth_end)


def get_fwhm(shiftedtth):
    u = 0.1791  # thats what icsd has for neutrons
    v = -0.4503  # should be in options dialog at later step
    w = 0.4
    return sqrt(u + v * tan(radians(shiftedtth) / 2.0)
                + w * tan(radians(shiftedtth) / 2.0) ** 2)


def get_peak_width(fwhm):
    return fwhm / (2 * sqrt(2 * log(2)))


def get_peak(peak_width, tth, shiftedtth):
    return (1 / (peak_width * sqrt(2 * pi))
            * np.exp(-0.5 * ((tth - shiftedtth) / peak_width) ** 2))


def get_intensity_prof(refl, shift, tth):
    shiftedtth = refl.tth + shift
    fwhm = get_fwhm(shiftedtth)
    peak_width = get_peak_width(fwhm)
    peak = get_peak(peak_width, tth, shiftedtth)
    prof = refl.fs_lc * peak
    return prof


def get_tth_range(start, end, shift):
    tth_step = 0.1
    tth_end = get_tth_end(end, shift)
    tth_start = get_tth_start(start, shift)
    tth = np.linspace(tth_start, tth_end,
                      int((tth_end - tth_start) / tth_step) + 1)
    return tth


def get_tth_bins(tth, tth_step=0.1):
    bins = tth - tth_step / 2.0
    bins = np.append(bins, tth[-1] + tth_step / 2.0)
    return bins


def get_unique_refl(refls):
    return [refl for refl in refls if refl.unique]


def get_inplane_refl(refls):
    return [refl for refl in refls if refl.inplane]


def get_unique_inplane_refl(refls):
    newrefls = []
    for refl in refls:
        if refl.inplane:
            if not any(refl.hkl in nre.equivalents for nre in newrefls):
                newrefls.append(refl)
    return newrefls

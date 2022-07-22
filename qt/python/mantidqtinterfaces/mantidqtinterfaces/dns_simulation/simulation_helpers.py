# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,

#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +

"""
Helper functions only for DNS simulation calculations.
"""

from numpy import arccos, arcsin, cos, log, pi, rad2deg, radians, sin, sqrt, \
    tan
import numpy as np
from mantid.kernel import V3D


def two_theta_to_q(two_theta, wavelength):
    q = 4.0 * pi * sin(two_theta / 2.0) / wavelength
    return q


def two_theta_to_d(two_theta, wavelength):
    d = wavelength / (2.0 * sin(two_theta / 2.0))
    return d


def angle_to_hkl(two_theta, omega, ub_inv, wavelength):
    """
    Return hkl for diffractometer angles and given inverse UB.
    """
    theta = radians(two_theta) / 2.0
    omega = radians(omega) - theta
    u_phi = np.array([-cos(omega), 0, -sin(omega)])
    h_phi = 2.0 * sin(theta) / wavelength * u_phi
    hkl = np.dot(ub_inv, h_phi)
    return hkl


def check_inplane(q1, q2, hkl):
    """
    Return True if hkl is in the plane specified by q1 and q2.
    """
    q1q2 = np.column_stack((q1, q2))
    is_inplane = np.linalg.det(np.column_stack((q1q2, hkl))) == 0
    return is_inplane


def det_rot_number_to_two_theta(det_rot, det_number):
    """
    Return 2theta for det_rot.
    """
    two_theta = -det_rot + det_number * 5.0
    return two_theta


def hkl_to_omega(hkl, ub, wavelength, two_theta):
    """
    Return omega for given hkl.
    """
    h_phi = np.transpose(np.dot(hkl, np.transpose(ub)))
    u_phi = h_phi * wavelength / 2.0 / sin(radians(two_theta / 2.0))
    sign = np.sign(arcsin(round(float(u_phi[2]), 10)))
    omega = -sign * rad2deg(arccos((round(float(u_phi[0]), 10)))) + two_theta / 2.0
    # rounding is necessary to prevent floats slightly larger than 1, which will
    # kill arccos
    return omega


def ki_from_wavelength(wavelength):
    """
    Sets k_i for given wavelength.
    """
    k_i = 2.0 * pi / wavelength
    return k_i


def list_to_v3d(input_list):
    """
    Make a mantid vector from a list.
    """
    v3dv = V3D(input_list[0], input_list[1], input_list[2])
    return v3dv


def lorentz_correction(intensity, two_theta):
    """
    Apply lorentz correction for simulation.
    """
    correction = intensity / sin(radians(two_theta))
    return correction


def max_q(wavelength):
    """
    Return maximum q at DNS for PA detector.
    """
    max_two_theta = 135.5
    maximum_q = 4.0 * pi * sin(radians(max_two_theta/ 2.0)) / wavelength
    return maximum_q


def omega_to_sample_rot(omega, det_rot):
    """
    Convert Omega to sample_rot.
    """
    sample_rot = omega + det_rot
    return sample_rot


def q_to_two_theta(q, wavelength):
    """
    Return two theta for given q.
    """
    two_theta = rad2deg(2.0 * arcsin(q * wavelength / (4.0 * pi)))
    return two_theta


def rotate_ub(omega_offset, ubm):
    """
    Rotate the UB matrix by given angle around z-axis.
    """
    w1 = radians(omega_offset)
    rm = np.asarray([[cos(w1), 0, -sin(w1)], [0, 1, 0], [sin(w1), 0, cos(w1)]])
    ubm = np.dot(rm, ubm)
    return ubm


def get_two_theta_path(two_theta_range, omega_range):
    """
    Returns two_theta list of the surface boundary of dns map.
    """
    two_theta_path = np.concatenate(
        (two_theta_range[0] * np.ones(len(omega_range)), two_theta_range,
         two_theta_range[-1] * np.ones(len(omega_range)), np.flip(two_theta_range)))
    return two_theta_path


def get_omega_path(two_theta_range, omega_range):
    """
    Returns omega list of the surface boundary of dns map.
    """
    omega_path = np.concatenate(
        (np.flip(omega_range), omega_range[0] * np.ones(len(two_theta_range)),
         omega_range, omega_range[-1] * np.ones(len(two_theta_range))))
    return omega_path


def return_dns_surface_shape(oriented_lattice, q1, q2, wavelength, options):
    """
    Returns qx, qy array of the surface boundary of dns single
    crystal measurement.
    """
    two_theta_start = options['sc_det_start']
    two_theta_end = options['sc_det_end']
    o_start = options['sc_sam_start']
    o_end = options['sc_sam_end']
    two_theta_range = get_two_theta_range_sc(two_theta_start, two_theta_end)
    omega_range = get_omega_range_single_crystal(o_start, o_end, two_theta_start)
    ubm_inv = np.linalg.inv(oriented_lattice.getUB())
    line = []
    # create shape of dns measured surface in two_theta and omega
    two_theta_r = get_two_theta_path(two_theta_range, omega_range)
    omega_r = get_omega_path(two_theta_range, omega_range)
    for i, omega in enumerate(omega_r):
        hkl = angle_to_hkl(two_theta_r[i], omega, ubm_inv, wavelength)
        qx, qy = return_qxqy(oriented_lattice, q1, q2, hkl)
        line.append([qx, qy])
    return np.asarray(line)


def get_angle_q_vs_hkl(hkl, q, oriented_lattice):
    angle = radians(oriented_lattice.recAngle(hkl[0], hkl[1], hkl[2], q[0], q[1], q[2]))
    return angle


def return_qxqy(oriented_lattice, q1, q2, hkl):
    """
    Returns projection of hkl along q1 and q2
    """
    angle_q1_hkl = get_angle_q_vs_hkl(hkl, q1, oriented_lattice)
    angle_q2_hkl = get_angle_q_vs_hkl(hkl, q2, oriented_lattice)
    n_q1 = np.linalg.norm(oriented_lattice.qFromHKL(q1))
    n_q2 = np.linalg.norm(oriented_lattice.qFromHKL(q2))
    n_hkl = np.linalg.norm(oriented_lattice.qFromHKL(hkl))
    qx = cos(angle_q1_hkl) * n_hkl / n_q1
    qy = cos(angle_q2_hkl) * n_hkl / n_q2
    return [qx, qy]


def return_qx_qx_inplane_refl(oriented_lattice, q1, q2, refls):
    """
    Returns a qx, qx, Intensity, h, k, l list for reflections in plane.
    """
    reflections = []
    for refl in refls:
        if refl.inplane:
            qx, qy = return_qxqy(oriented_lattice, q1, q2, refl.hkl)
            intensity = refl.fs_lc
            reflections.append([qx, qy, intensity, refl.h, refl.k, refl.l])
    return np.asarray(reflections)


def shift_channels_below_23(channel, det_rot):
    """
    Shifting detector_rot to higher angles
    if reflection is not on last detector.
    """
    while channel > 23:
        channel += -1
        det_rot += -5
    return [channel, det_rot]


def two_theta_to_rot_number(two_theta):
    """
    Return det_rot and detector number for given two theta.
    """
    det_rot = -5 - two_theta % 5
    channel = two_theta // 5 - 1
    return [det_rot, channel]


def get_omega_range_single_crystal(o_start, o_end, two_theta_start):
    omega_start = o_start - two_theta_start
    omega_end = o_end - two_theta_start
    omega_range_sc = get_1deg_angle_linspace(omega_start, omega_end)
    return omega_range_sc


def get_two_theta_end(end, shift):
    two_theta_end = -end + 23 * 5 + shift
    return two_theta_end


def get_two_theta_start(start, shift):
    two_theta_start = -start + shift
    return two_theta_start


def get_1deg_angle_linspace(start, end):
    angle_linspace = np.linspace(start, end, int(abs((end - start))) + 1)
    return angle_linspace


def get_two_theta_range_sc(start, end):
    two_theta_start = get_two_theta_start(start, 0)
    two_theta_end = get_two_theta_end(end, 0)
    two_theta_range = get_1deg_angle_linspace(two_theta_start, two_theta_end)
    return two_theta_range


def get_fwhm(shifted_two_theta):
    u = 0.1791  # that's what icsd has for neutrons
    v = -0.4503  # should be in options dialog at later step
    w = 0.4
    fwhm = sqrt(u + v * tan(radians(shifted_two_theta) / 2.0)
                + w * tan(radians(shifted_two_theta) / 2.0) ** 2)
    return fwhm


def get_peak_width(fwhm):
    peak_width = fwhm / (2.0 * sqrt(2.0 * log(2.0)))
    return peak_width


def get_peak(peak_width, two_theta, shifted_two_theta):
    peak = (1.0 / (peak_width * sqrt(2.0 * pi))
            * np.exp(-0.5 * ((two_theta - shifted_two_theta) / peak_width) ** 2))
    return peak


def get_intensity_prof(refl, shift, two_theta):
    shifted_two_theta = refl.tth + shift
    fwhm = get_fwhm(shifted_two_theta)
    peak_width = get_peak_width(fwhm)
    peak = get_peak(peak_width, two_theta, shifted_two_theta)
    prof = refl.fs_lc * peak
    return prof


def get_two_theta_range(start, end, shift):
    two_theta_step = 0.1
    two_theta_end = get_two_theta_end(end, shift)
    two_theta_start = get_two_theta_start(start, shift)
    two_theta = np.linspace(two_theta_start, two_theta_end,
                            int((two_theta_end - two_theta_start) / two_theta_step) + 1)
    return two_theta


def get_two_theta_bins(two_theta, two_theta_step=0.1):
    bins = two_theta - two_theta_step / 2.0
    bins = np.append(bins, two_theta[-1] + two_theta_step / 2.0)
    return bins


def get_unique_refl(refls):
    unique_refl = [refl for refl in refls if refl.unique]
    return unique_refl


def get_inplane_refl(refls):
    inplane_refl = [refl for refl in refls if refl.inplane]
    return inplane_refl


def get_unique_inplane_refl(refls):
    new_refls = []
    for refl in refls:
        if refl.inplane:
            if not any(refl.hkl in nre.equivalents for nre in new_refls):
                new_refls.append(refl)
    return new_refls

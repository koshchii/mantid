# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,

#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +

"""
DNS single crystal elastic plot tab model of DNS reduction GUI.
"""

import numpy as np

import mantidqtinterfaces.dns_single_crystal_elastic.plot.elastic_single_crystal_helpers as helper
from mantidqtinterfaces.dns_powder_tof.data_structures.dns_obs_model import DNSObsModel
from mantidqtinterfaces.dns_single_crystal_elastic.data_structures.dns_single_crystal_map import DNSScMap
from mantidqtinterfaces.dns_powder_tof.data_structures.object_dict import ObjectDict
from mantidqtinterfaces.dns_single_crystal_elastic.plot.elastic_single_crystal_helpers import get_hkl_intensity_from_cursor


class DNSElasticSCPlotModel(DNSObsModel):
    """
    Model for DNS plot calculations. It converts data from (omega, 2theta)
    space into (q_x, q_y) and (n_x, n_y). Also, creates hull of DNS data
    to filter points.
    """

    def __init__(self, parent):
        super().__init__(parent)
        self._single_crystal_map = None

        self._data = ObjectDict()
        self._data.x = None
        self._data.y = None
        self._data.z = None
        self._data.z_min = None
        self._data.z_max = None
        self._data.pz_min = None
        self._data.triang = None
        self._data.z_triang = None

    def create_single_crystal_map(self, data_array, options, initial_values=None):
        two_theta = data_array["two_theta_array"]
        omega = data_array["omega_array"]
        z_mesh = data_array["intensity"]
        error = data_array["error"]
        parameter = {
            "wavelength": options["wavelength"],
            "dx": options["dx"],
            "dy": options["dy"],
            "hkl1": options["hkl1"],
            "hkl2": options["hkl2"],
            "omega_offset": options["omega_offset"],
        }
        if initial_values is not None:
            parameter.update(initial_values)
        self._single_crystal_map = DNSScMap(two_theta=two_theta, omega=omega, z_mesh=z_mesh, error_mesh=error, parameter=parameter)
        return self._single_crystal_map

    def has_data(self):
        return self._data.x is not None

    def get_projections(self, xlim, ylim):
        limits = np.append(xlim, ylim)
        x, y, z = helper.filter_flattened_meshes(self._data.x, self._data.y, self._data.z, limits)
        x_projection = helper.get_projection(x, z)
        y_projection = helper.get_projection(y, z)
        return x_projection, y_projection

    def get_interpolated_quadmesh(self, interpolate, axis_type):
        new_mesh_name = axis_type + "_mesh_interpolated"
        old_mesh_name = axis_type + "_mesh"
        if interpolate:
            self._single_crystal_map.interpolate_quad_mesh(interpolate)
            x, y, z = getattr(self._single_crystal_map, new_mesh_name)
        else:
            x, y, z = getattr(self._single_crystal_map, old_mesh_name)
        self._data.x = x
        self._data.y = y
        self._data.z = z
        return x, y, z

    def get_interpolated_triangulation(self, interpolate, axis_type, switch):
        mesh_name = axis_type + "_mesh"
        self._single_crystal_map.triangulate(mesh_name=mesh_name, switch=switch)
        self._single_crystal_map.mask_triangles(mesh_name=mesh_name)
        triangulator_refiner, z_refiner = self._single_crystal_map.interpolate_triangulation(interpolate)
        self._data.triang = triangulator_refiner
        self._data.z_triang = z_refiner
        # this is important to get the limits
        x, y, z = getattr(self._single_crystal_map, mesh_name)
        self._data.x = x
        self._data.y = y
        self._data.z = z
        return triangulator_refiner, z_refiner

    def get_xy_dy_ratio(self):
        return self._single_crystal_map.dx / self._single_crystal_map.dy

    def get_aspect_ratio(self, axis_type):
        if axis_type["fix_aspect"]:
            if axis_type["type"] == "hkl":
                ratio = self.get_xy_dy_ratio()
                return ratio
            return 1
        return "auto"

    def get_axis_labels(self, axis_type, crystal_axes, switch=False):
        if crystal_axes:
            return self._single_crystal_map.get_crystal_axis_names()
        hkl1 = self._single_crystal_map.hkl1
        hkl2 = self._single_crystal_map.hkl2
        axis_labels = {
            "two_theta_and_omega": ["2\u03B8 (deg)", "\u03C9 (deg)"],
            "qxqy": [r"$q_{x} \ (\AA^{-1})$", r"$q_{y} \ (\AA^{-1})$"],
            "hkl": [f"[{hkl1}] (r.l.u.)", f"[{hkl2}] (r.l.u.)"],
        }
        labels = axis_labels[axis_type]
        if switch:
            labels.reverse()
        return labels

    def get_changing_hkl_components(self):
        return self._single_crystal_map.get_changing_hkl_components()

    def get_format_coord(self, axis_type):
        # adds z and hkl label to cursor position
        # this is somehow bad since it backlinks the model function
        # from the view
        def format_coord(x, y):
            h, k, l, z, error = get_hkl_intensity_from_cursor(self._single_crystal_map, axis_type, x, y)
            return f"x = {x: 2.4f}, y = {y: 2.4f}, hkl = ({h: 2.2f}, {k: 2.2f}, {l: 2.2f}),\n Intensity = {z: 6.4f} ± {error:6.4f}"

        return format_coord

    @staticmethod
    def get_m_limits(*args):
        return [helper.string_range_to_float(arg) for arg in args]

    @staticmethod
    def get_mz_limit(arg):
        return helper.string_range_to_float(arg)

    def get_data_z_min_max(self, xlim=None, ylim=None):
        return helper.get_z_min_max(self._data.z, xlim, ylim, self._data.x, self._data.y)

    def get_data_xy_lim(self, switch):
        limits = [[min(self._data.x.flatten()), max(self._data.x.flatten())], [min(self._data.y.flatten()), max(self._data.y.flatten())]]
        if switch:
            limits.reverse()
        return limits

    @staticmethod
    def switch_axis(x, y, z, switch):
        if switch:  # switch x and y axes
            nx = np.transpose(y)
            ny = np.transpose(x)
            nz = np.transpose(z)
            return nx, ny, nz
        return x, y, z

    def get_omega_offset(self):
        return self._single_crystal_map["omega_offset"]

    def get_dx_dy(self):
        return self._single_crystal_map["dx"], self._single_crystal_map["dy"]

    def prepare_data_for_saving(self):
        x = self._data.x.flatten()
        y = self._data.y.flatten()
        z = self._data.z.flatten()
        data_combined = np.array(list(zip(x, y, z)))
        return data_combined

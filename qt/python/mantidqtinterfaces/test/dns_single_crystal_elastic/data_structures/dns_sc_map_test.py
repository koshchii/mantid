# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
import unittest
from unittest import mock
from unittest.mock import patch
import numpy as np
from mantidqtinterfaces.dns_sc_elastic.data_structures.dns_sc_map import \
    DNSScMap
from mantidqtinterfaces.dns_powder_tof.helpers.helpers_for_testing import \
    get_fake_elastic_sc_dataset, \
    get_fake_elastic_sc_options
from mantidqtinterfaces.dns_powder_tof.data_structures.object_dict \
    import ObjectDict
from mantidqtinterfaces.dns_sc_elastic.data_structures.dns_sc_map \
    import (_get_mesh,
            _is_rectangular_mesh,
            _correct_rect_grid,
            _correct_omegaoffset,
            _get_unique,
            _get_q_mesh,
            _get_hkl_mesh,
            _get_interpolated)


class DNSScMapTest(unittest.TestCase):
    # pylint: disable=protected-access, too-many-public-methods, too-many-arguments

    @classmethod
    def setUpClass(cls):
        data_array = get_fake_elastic_sc_dataset()
        options = get_fake_elastic_sc_options()
        ttheta = data_array['ttheta']
        omega = data_array['omega']
        z_mesh = data_array['intensity']
        error = data_array['error']
        parameter = {'wavelength': options['wavelength'],
                     'dx': options['dx'],
                     'dy': options['dy'],
                     'hkl1': options['hkl1'],
                     'hkl2': options['hkl2'],
                     'omega_offset': options['omega_offset']
                     }
        cls.map = DNSScMap(ttheta=ttheta,
                           omega=omega,
                           z_mesh=z_mesh,
                           error_mesh=error,
                           parameter=parameter)

    def setUp(self):
        self.map.hkl1 = '1,2,3'
        self.map.hkl2 = '2,3,4'

    def test___init__(self):
        self.assertIsInstance(self.map, DNSScMap)
        self.assertIsInstance(self.map, ObjectDict)
        self.assertTrue(hasattr(self.map, 'omega_intp'))
        self.assertEqual(self.map.rectangular_grid, True)
        self.assertTrue((self.map.ttheta == [0, 1, 2]).all())
        self.assertTrue((self.map.omega == [4, 5]).all())
        self.assertEqual(self.map.omegaoffset, 0)
        self.assertEqual(self.map.dx, 1)
        self.assertEqual(self.map.dy, 2)
        self.assertEqual(self.map.hkl1, '1,2,3')
        self.assertEqual(self.map.hkl2, '2,3,4')
        self.assertEqual(self.map.wavelength, 4.74)
        self.assertIsInstance(self.map.tthomega_mesh, list)
        self.assertIsInstance(self.map.hkl_mesh, list)
        self.assertIsInstance(self.map.qxqy_mesh, list)
        self.assertEqual(len(self.map.tthomega_mesh), 3)
        self.assertEqual(len(self.map.hkl_mesh), 3)
        self.assertEqual(len(self.map.qxqy_mesh), 3)
        self.assertTrue(
            (self.map.omega_mesh == np.array([[4, 5], [4, 5], [4, 5]])).all())
        self.assertTrue(
            (self.map.ttheta_mesh == np.array([[0, 0], [1, 1], [2, 2]])).all())
        self.assertTrue(np.allclose(self.map.qx_mesh,
                                    [[0., 0.],
                                     [-0.00141237, -0.00181517],
                                     [-0.00242151, -0.00322754]]))
        self.assertTrue(np.allclose(self.map.qy_mesh,
                                    [[0., 0.],
                                     [-0.02309205, -0.02306389],
                                     [-0.04620524, -0.04615594]]))
        self.assertTrue(np.allclose(self.map.hklx_mesh,
                                    [[0., 0.],
                                     [-0.00022479, -0.00028889],
                                     [-0.0003854, -0.00051368]]))
        self.assertTrue(np.allclose(self.map.hkly_mesh,
                                    [[0., 0.],
                                     [-0.00735043, -0.00734146],
                                     [-0.01470759, -0.01469189]]))
        self.assertTrue(
            (self.map.z_mesh == np.array(
                [[8.0, 11.0], [9.0, 12.0], [10.0, 13.0]])).all())
        self.assertTrue(
            (self.map.error_mesh == np.array(
                [[14.0, 17.0], [15.0, 18.0], [16.0, 19.0]])).all())

    def test__get_mesh(self):
        data_array = get_fake_elastic_sc_dataset()
        ttheta = data_array['ttheta']
        omega = data_array['omega']
        z_mesh = data_array['intensity']
        z_mesh[0, 1] = np.nan
        testv = _get_mesh(omega, ttheta, z_mesh)
        self.assertTrue((testv[0] == np.array([4, 4, 5, 4, 5])).all())
        self.assertTrue((testv[1] == np.array([0, 1, 1, 2, 2])).all())
        self.assertTrue((testv[2] == np.array([8, 9, 12, 10, 13])).all())

    def test__is_rectangular_mesh(self):
        data_array = get_fake_elastic_sc_dataset()
        ttheta = data_array['ttheta']
        omega = data_array['omega']
        z_mesh = data_array['intensity']
        testv = _is_rectangular_mesh(omega, ttheta, z_mesh)
        self.assertTrue(testv)
        omega = [1]
        testv = _is_rectangular_mesh(omega, ttheta, z_mesh)
        self.assertFalse(testv)

    def test__correct_rect_grid(self):
        data_array = get_fake_elastic_sc_dataset()
        ttheta = data_array['ttheta']
        omega = data_array['omega']
        z_mesh = data_array['intensity']
        omega_mesh, ttheta_mesh = 0, 1
        testv = _correct_rect_grid(z_mesh, omega_mesh, ttheta_mesh, omega,
                                   ttheta)
        self.assertEqual(testv[0].shape, (3, 2))
        self.assertEqual(testv[1].shape, (3, 2))
        self.assertEqual(testv[2].shape, (3, 2))
        self.assertTrue(testv[3])
        omega = [1]
        testv = _correct_rect_grid(z_mesh, omega_mesh, ttheta_mesh, omega,
                                   ttheta)
        self.assertTrue((testv[0] == z_mesh).all())
        self.assertEqual(testv[1], 0)
        self.assertEqual(testv[2], 1)
        self.assertFalse(testv[3])

    def test__correct_omegaoffset(self):
        data_array = get_fake_elastic_sc_dataset()
        omega = data_array['omega']
        testv = _correct_omegaoffset(omega, 3)
        self.assertEqual(testv[0], 1)
        self.assertEqual(testv[1], 2)

    def test_get_unique(self):
        testv = _get_unique([1, 1, 2], [3, 3, 3])
        self.assertTrue((testv[0] == [1, 2]).all())
        self.assertTrue((testv[1] == [3]).all())

    @patch('mantidqtinterfaces.dns_sc_elastic.data_structures.dns_sc_map.'
           'angle_to_q')
    def test__get_q_mesh(self, mock_angle):
        testv = _get_q_mesh(1, 2, 3)
        mock_angle.assert_called_once_with(ttheta=2,
                                           omega=1,
                                           wavelength=3)
        self.assertEqual(testv, mock_angle.return_value)

    def test__get_hkl_mesh(self):
        testv = _get_hkl_mesh(2, 3, 4, 5)
        self.assertAlmostEqual(testv[0], 1.2732395447351628)
        self.assertAlmostEqual(testv[1], 2.3873241463784303)

    def test__get_interpolated(self):
        testar = np.array([1, 2, 3])
        testv = _get_interpolated(testar, 2)
        car = [1., 1.25, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3.]
        self.assertTrue((testv == car).all())

    def test_create_np_array(self):
        self.map.omega_intp = None
        self.map.rectangular_grid = True
        testv = self.map.create_np_array()
        self.assertIsNone(testv[1])
        car = np.array([[0, 4, 0, 0, 0, 0, 8, 14],
                        [0, 5, 0, 0, 0, 0, 11, 17],
                        [1, 4, -1.41237e-03, -2.30920e-02, -2.24785e-04,
                         -7.35042e-03, 9, 15],
                        [1, 5, -1.81516e-03, -2.30638e-02, -2.88892e-04,
                         -7.34146e-03, 12, 18],
                        [2, 4, -2.42151e-03, -4.62052e-02, -3.85395e-04,
                         -1.47075e-02, 10, 16],
                        [2, 5, -3.22753e-03, -4.61559e-02, -5.13678e-04,
                         -1.46918e-02, 13, 19]])
        self.assertTrue(testv[0].shape == (6, 8))
        self.assertTrue(np.allclose(testv[0], car))
        self.map.omega_intp = True
        self.map.rectangular_grid = True
        self.map.omega_mesh_intp = self.map.omega_mesh
        self.map.ttheta_mesh_intp = self.map.ttheta_mesh
        self.map.qx_mesh_intp = self.map.qx_mesh
        self.map.qy_mesh_intp = self.map.qy_mesh
        self.map.hklx_mesh_intp = self.map.hklx_mesh
        self.map.hkly_mesh_intp = self.map.hkly_mesh
        self.map.z_mesh_intp = self.map.z_mesh
        testv = self.map.create_np_array()
        self.assertTrue(testv[1].shape == (6, 7))
        self.assertTrue(np.allclose(testv[1], car[:, 0:7]))

    @patch('mantidqtinterfaces.dns_sc_elastic.data_structures.dns_sc_map.'
           'file_helper.create_dir_from_filename')
    @patch('mantidqtinterfaces.dns_sc_elastic.data_structures.dns_sc_map.'
           'DNSScMap.create_np_array')
    @patch('mantidqtinterfaces.dns_sc_elastic.data_structures.dns_sc_map.'
           'np.savetxt')
    def test_save_ascii(self, mock_save, mock_create_array, mock_cdir):
        header = ' 2theta, omega, qx,      qy ,      hklx,' \
                 '     hkly,         Intensity,          Error'
        self.map.rectangular_grid = False
        mock_create_array.return_value = 1, 2
        self.map.save_ascii('123.txt')
        mock_cdir.assert_called_once_with('123_no_interp.txt')
        mock_save.assert_called_once_with(
            '123_no_interp.txt',
            1,
            fmt='%7.3f %7.3f %8.5f %8.5f %9.5f %9.5f %15.5f %15.5f',
            delimiter=' ',
            newline='\n',
            header=header
        )
        mock_save.reset_mock()
        self.map.rectangular_grid = True
        self.map.save_ascii('123.txt')
        self.assertEqual(mock_save.call_count, 2)
        mock_save.assert_called_with(
            '123_interp.txt',
            2,
            fmt='%7.3f %7.3f %8.5f %8.5f %9.5f %9.5f %15.5f',
            delimiter=' ',
            newline='\n',
            header=header[0:-16]
        )

    @patch('mantidqtinterfaces.dns_sc_elastic.data_structures.dns_sc_map.'
           'scipy')
    def test__get_z_mesh_interp(self, mock_scipy):
        self.map.omega_intp = 3
        self.map.ttheta_intp = 4
        testv = self.map._get_z_mesh_interp()
        mock_scipy.interpolate.interp2d.assert_called_once()
        cargs = mock_scipy.interpolate.interp2d.call_args_list[0][0]
        self.assertTrue(np.allclose(cargs[0], np.asarray([4, 5])))
        self.assertTrue(np.allclose(cargs[1], np.asarray([0, 1, 2])))
        zara = np.array([[8.0, 11.0], [9.0, 12.0], [10.0, 13.0]])
        self.assertTrue(np.allclose(cargs[2], zara))
        mock_scipy.interpolate.interp2d.return_value.assert_called_once_with(
            3, 4)
        self.assertEqual(
            testv, mock_scipy.interpolate.interp2d.return_value.return_value)

    @patch('mantidqtinterfaces.dns_sc_elastic.data_structures.dns_sc_map.'
           'np.meshgrid')
    @patch('mantidqtinterfaces.dns_sc_elastic.data_structures.dns_sc_map.'
           'DNSScMap._get_z_mesh_interp')
    @patch('mantidqtinterfaces.dns_sc_elastic.data_structures.dns_sc_map.'
           '_get_hkl_mesh')
    @patch('mantidqtinterfaces.dns_sc_elastic.data_structures.dns_sc_map.'
           '_get_q_mesh')
    @patch('mantidqtinterfaces.dns_sc_elastic.data_structures.dns_sc_map.'
           '_get_interpolated')
    def test_interpolate_quad_mesh(self, mock_get_interpolated,
                                   mock_get_q_mesh, mock_get_hkl_mesh,
                                   mock_get_z_mesh_interp, mock_meshgrid):
        self.map.rectangular_grid = False
        mock_get_q_mesh.return_value = 1, 2
        mock_get_hkl_mesh.return_value = 3, 4
        mock_get_z_mesh_interp.return_value = 8
        mock_meshgrid.return_value = 4, 5
        self.map.interpolate_quad_mesh()
        mock_get_interpolated.assert_not_called()
        self.map.rectangular_grid = True
        self.map.interpolate_quad_mesh(interp=0)
        mock_get_interpolated.assert_not_called()
        mock_get_interpolated.reset_mock()
        self.map.interpolate_quad_mesh()
        carg = mock_get_interpolated.call_args_list
        self.assertEqual(len(carg), 2)
        self.assertTrue(np.allclose(carg[0][0][0], np.asarray([0, 1, 2])))
        self.assertEqual(carg[0][0][1], 3)
        self.assertTrue(np.allclose(carg[1][0][0], np.asarray([4, 5])))
        self.assertEqual(carg[1][0][1], 3)
        self.assertEqual(self.map.ttheta_intp,
                         mock_get_interpolated.return_value)
        self.assertEqual(self.map.omega_intp,
                         mock_get_interpolated.return_value)
        mock_meshgrid.assert_called_with(mock_get_interpolated.return_value,
                                         mock_get_interpolated.return_value)

        mock_get_q_mesh.assert_called_once_with(4, 5, 4.74)
        mock_get_hkl_mesh.assert_called_once_with(1, 2, 1, 2)
        mock_get_z_mesh_interp.assert_called_once()
        self.assertEqual(self.map.z_mesh_intp, 8)
        self.assertEqual(self.map.omega_mesh_intp, 4)
        self.assertEqual(self.map.ttheta_mesh_intp, 5)
        self.assertEqual(self.map.qx_mesh_intp, 1)
        self.assertEqual(self.map.qy_mesh_intp, 2)
        self.assertEqual(self.map.hklx_mesh_intp, 3)
        self.assertEqual(self.map.hkly_mesh_intp, 4)
        self.assertEqual(self.map.tthomega_mesh_intp, [5, 4, 8])
        self.assertEqual(self.map.hkl_mesh_intp, [3, 4, 8])
        self.assertEqual(self.map.qxqy_mesh_intp, [1, 2, 8])

    @patch('mantidqtinterfaces.dns_sc_elastic.data_structures.dns_sc_map.'
           'tri')
    def test_triangulate(self, mock_tri):
        self.map.hkl_mesh = [np.asarray([0, 1]), np.asarray([2, 3]),
                             np.asarray([4, 5])]
        testv = self.map.triangulate('hkl_mesh', switch=False)
        mock_tri.Triangulation.assert_called_once()
        self.assertTrue(
            (mock_tri.Triangulation.call_args_list[0][0][0]
             == np.asarray([0, 1])).all())
        self.assertTrue(
            (mock_tri.Triangulation.call_args_list[0][0][1]
             == np.asarray([2, 3])).all())
        self.assertEqual(testv, mock_tri.Triangulation.return_value)
        mock_tri.reset_mock()
        testv = self.map.triangulate('hkl_mesh', switch=True)
        self.assertTrue(
            (mock_tri.Triangulation.call_args_list[0][0][0]
             == np.asarray([2, 3])).all())

    @patch('mantidqtinterfaces.dns_sc_elastic.data_structures.dns_sc_map.'
           'UniformTriRefiner')
    @patch('mantidqtinterfaces.dns_sc_elastic.data_structures.dns_sc_map.'
           'LinearTriInterpolator')
    def test_interpolate_triangulation(self, mock_lin_intp, mock_uni):
        refimock = mock.Mock()
        mock_uni.return_value.refine_field.return_value = 4, refimock
        zmock = mock.Mock()
        self.map.triang = None
        self.map.z_mesh = zmock
        self.assertIsNone(self.map.interpolate_triangulation())
        self.map.triang = 1
        testv = self.map.interpolate_triangulation()
        mock_lin_intp.assert_called_once_with(1, zmock.flatten.return_value)
        mock_uni.return_value.refine_field.assert_called_once_with(
            zmock, subdiv=0, triinterpolator=mock_lin_intp.return_value)
        self.assertEqual(testv[0], 1)
        self.assertEqual(testv[1], zmock.flatten.return_value)
        testv = self.map.interpolate_triangulation(interp=1)
        self.assertEqual(testv[0], 4)
        self.assertEqual(testv[1], refimock.flatten.return_value)

    @patch('mantidqtinterfaces.dns_sc_elastic.data_structures.dns_sc_map.path')
    def test_get_dns_map_border(self, mock_path):
        testv = self.map.get_dns_map_border('qxqy')
        self.assertEqual(testv, mock_path.Path.return_value)
        testarray = np.array([[0., 0.],
                              [0., 0.],
                              [0., 0.],
                              [-0.00141237, -0.02309205],
                              [-0.00242151, -0.04620524],
                              [-0.00242151, -0.04620524],
                              [-0.00322754, -0.04615594],
                              [-0.00322754, -0.04615594],
                              [-0.00181517, -0.02306389],
                              [0., 0.]])
        self.assertTrue(np.allclose(
            mock_path.Path.call_args_list[0][0][0], testarray))
        mock_path.reset_mock()
        self.map.get_dns_map_border('hkl')
        testarray = np.array([[0., 0.],
                              [0., 0.],
                              [0., 0.],
                              [-0.00022479, -0.00735043],
                              [-0.0003854, -0.01470759],
                              [-0.0003854, -0.01470759],
                              [-0.00051368, -0.01469189],
                              [-0.00051368, -0.01469189],
                              [-0.00028889, -0.00734146],
                              [0., 0.]])
        self.assertTrue(np.allclose(
            mock_path.Path.call_args_list[0][0][0], testarray))

    def test_mask_triangles(self):
        mock_triang = mock.Mock()
        self.map.triang = mock_triang
        self.map.triang.triangles = np.array([[5, 4, 2],
                                              [2, 4, 0],
                                              [5, 2, 3],
                                              [3, 2, 0]])
        testv = self.map.mask_triangles('tthomega')
        self.assertEqual(testv, mock_triang)
        testv = self.map.mask_triangles('qxqy_mesh')
        mock_triang.set_mask.assert_called_once()
        carg = mock_triang.set_mask.call_args_list[0][0][0]
        self.assertTrue((carg == np.array([False, True, False, False])).all())

    def test_return_changing_indexes(self):
        testv = self.map.return_changing_indexes()
        self.assertEqual(testv, [1, 2])  # 3 changing indexes, does not work
        self.map.hkl1 = '0, 1, 0'
        self.map.hkl2 = '1, 1, 0'
        testv = self.map.return_changing_indexes()
        self.assertEqual(testv, [0, 1])
        self.map.hkl1 = '0, 1, 2'
        self.map.hkl2 = '0, 1, 1'
        testv = self.map.return_changing_indexes()
        self.assertEqual(testv, [1, 2])

    def test_get_crystal_axis_names(self):
        testv = self.map.get_crystal_axis_names()
        self.assertEqual(testv, ('k (r.l.u)', 'l (r.l.u)'))

    def test_get_changing_hkl_components(self):
        testv = self.map.get_changing_hkl_components()
        self.assertEqual(testv, (2.0, 3.0, 3.0, 4.0))


if __name__ == '__main__':
    unittest.main()

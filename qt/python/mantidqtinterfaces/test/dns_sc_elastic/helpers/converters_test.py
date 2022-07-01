# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +

import unittest

from mantidqtinterfaces.dns_sc_elastic.helpers.converters import \
    d_spacing_from_lattice


class DNSconvertersTest(unittest.TestCase):
    def setUp(self):
        pass

    def test_d_spacing_from_lattice(self):
        testv = d_spacing_from_lattice(2, 3, 4, 78, 86, 85, [1, 2, 3])
        self.assertAlmostEqual(testv, 0.9979628311312633)


if __name__ == '__main__':
    unittest.main()

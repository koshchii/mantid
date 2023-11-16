# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
import unittest
from unittest import mock

from mantidqtinterfaces.dns_powder_tof.data_structures.dns_obs_model import \
    DNSObsModel
from mantidqtinterfaces.dns_powder_tof.options.common_options_model import \
    DNSCommonOptionsModel
from mantidqtinterfaces.dns_sc_elastic.options.elastic_sc_options_model \
    import DNSElasticSCOptionsModel
from mantidqtinterfaces.dns_powder_tof.helpers.helpers_for_testing import \
    get_fake_elastic_sc_options


class DNSElasticSCOptionsModelTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        parent = mock.Mock()
        cls.model = DNSElasticSCOptionsModel(parent)

    def test___init__(self):
        self.assertIsInstance(self.model, DNSElasticSCOptionsModel)
        self.assertIsInstance(self.model, DNSCommonOptionsModel)
        self.assertIsInstance(self.model, DNSObsModel)

    def test_get_dx_dy(self):
        options = get_fake_elastic_sc_options()
        testv = self.model.get_dx_dy(options)
        self.assertAlmostEqual(testv[0], 0.9979628311312633)
        self.assertAlmostEqual(testv[1], 0.6427233372178484)


if __name__ == '__main__':
    unittest.main()

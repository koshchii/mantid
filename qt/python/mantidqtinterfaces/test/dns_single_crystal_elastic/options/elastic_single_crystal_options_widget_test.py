# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2024 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +

import unittest
from unittest import mock

from mantidqt.gui_helper import get_qapplication

from mantidqtinterfaces.dns_powder_tof.data_structures.dns_widget import DNSWidget
from mantidqtinterfaces.dns_single_crystal_elastic.options.elastic_single_crystal_options_model import DNSElasticSCOptionsModel
from mantidqtinterfaces.dns_single_crystal_elastic.options.elastic_single_crystal_options_presenter import DNSElasticSCOptionsPresenter
from mantidqtinterfaces.dns_single_crystal_elastic.options.elastic_single_crystal_options_view import DNSElasticSCOptionsView
from mantidqtinterfaces.dns_single_crystal_elastic.options.elastic_single_crystal_options_widget import DNSElasticSCOptionsWidget

app, within_mantid = get_qapplication()


class DNSElasticSCOptionsWidgetTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        parent = mock.Mock()
        parent.view = None
        cls.widget = DNSElasticSCOptionsWidget("elastic_single_crystal_options", parent)

    def test___init__(self):
        self.assertIsInstance(self.widget, DNSElasticSCOptionsWidget)
        self.assertIsInstance(self.widget, DNSWidget)
        self.assertIsInstance(self.widget.view, DNSElasticSCOptionsView)
        self.assertIsInstance(self.widget.model, DNSElasticSCOptionsModel)
        self.assertIsInstance(self.widget.presenter, DNSElasticSCOptionsPresenter)


if __name__ == "__main__":
    unittest.main()

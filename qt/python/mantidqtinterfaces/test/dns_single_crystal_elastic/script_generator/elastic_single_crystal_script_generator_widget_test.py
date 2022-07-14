# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
import unittest
from unittest import mock

from mantidqt.gui_helper import get_qapplication

# yapf: disable
from mantidqtinterfaces.dns_powder_tof.data_structures.dns_widget \
    import DNSWidget
from mantidqtinterfaces.dns_powder_tof.script_generator.\
    common_script_generator_view import DNSScriptGeneratorView
from mantidqtinterfaces.dns_sc_elastic.script_generator.\
    elastic_sc_script_generator_model import DNSElasticSCScriptGeneratorModel
from mantidqtinterfaces.dns_sc_elastic.script_generator.\
    elastic_sc_script_generator_presenter \
    import DNSElasticSCScriptGeneratorPresenter  # yapf: disable
from mantidqtinterfaces.dns_sc_elastic.script_generator.\
    elastic_sc_script_generator_widget import \
    DNSElasticSCScriptGeneratorWidget

# yapf: enable

app, within_mantid = get_qapplication()


class DNSElasticSCScriptGeneratorWidgetTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        parent = mock.Mock()
        parent.view = None
        cls.widget = DNSElasticSCScriptGeneratorWidget(
            'elastic_sc_script_generator', parent)

    def test___init__(self):
        self.assertIsInstance(self.widget, DNSElasticSCScriptGeneratorWidget)
        self.assertIsInstance(self.widget, DNSWidget)
        self.assertIsInstance(self.widget.view, DNSScriptGeneratorView)
        self.assertIsInstance(self.widget.model,
                              DNSElasticSCScriptGeneratorModel)
        self.assertIsInstance(self.widget.presenter,
                              DNSElasticSCScriptGeneratorPresenter)


if __name__ == '__main__':
    unittest.main()

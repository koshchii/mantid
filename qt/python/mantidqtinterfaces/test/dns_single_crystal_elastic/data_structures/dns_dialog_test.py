# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +

import unittest
from unittest import mock
from unittest.mock import patch

from mantidqt.gui_helper import get_qapplication
from qtpy.QtWidgets import QDialog

from mantidqtinterfaces.dns_single_crystal_elastic.data_structures.dns_dialog \
    import DNSDialog

app, within_mantid = get_qapplication()


class DNSDialogTest(unittest.TestCase):

    @classmethod
    @patch('mantidqtinterfaces.dns_single_crystal_elastic.data_structures.'
           'dns_dialog.loadUi')
    @patch('mantidqtinterfaces.dns_single_crystal_elastic.data_structures.'
           'dns_dialog.Path')
    def setUpClass(cls, mock_path, mock_loadui):
        # pylint: disable=arguments-differ
        cls.path_parent = mock.Mock()
        cls.path_parent.absolute = mock.Mock()
        cls.path_parent.absolute.return_value = '456'
        mock_loadui.return_value = '123'
        cls.mock_loadui = mock_loadui
        mock_path.return_value = cls.path_parent
        cls.dialog = DNSDialog(files='123', ui='456')

    def test___init__(self):
        self.assertIsInstance(self.dialog, DNSDialog)
        self.assertIsInstance(self.dialog, QDialog)
        self.assertTrue(hasattr(self.dialog, '_content'))
        self.mock_loadui.assert_called_once()


if __name__ == '__main__':
    unittest.main()

# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +
from mantid.api import ITableWorkspace
from Muon.GUI.Common.ADSHandler.ADS_calls import check_if_workspace_exist, retrieve_ws
from Muon.GUI.Common.contexts.corrections_context import CorrectionsContext
from Muon.GUI.Common.contexts.muon_data_context import MuonDataContext
from Muon.GUI.Common.utilities.run_string_utils import run_list_to_string

DEAD_TIME_FROM_FILE = "FromFile"
DEAD_TIME_FROM_ADS = "FromADS"
DEAD_TIME_TABLE_KEY = "dead-time"


class CorrectionsModel:
    """
    The CorrectionsModel calculates Dead Time corrections and Background corrections.
    """

    def __init__(self, data_context: MuonDataContext, corrections_context: CorrectionsContext):
        """Initialize the CorrectionsModel with empty data."""
        self._data_context = data_context
        self._corrections_context = corrections_context

    @property
    def number_of_run_strings(self) -> int:
        """Returns the number of runs currently loaded into the context."""
        return len(self._data_context.current_runs)

    def run_number_strings(self) -> list:
        """Returns a list of run number strings from the context. Its a string because of co-add mode."""
        return [run_list_to_string(run_list) for run_list in self._data_context.current_runs]

    def current_runs(self) -> list:
        """Returns the currently selected run number in a list. If co-add is selected, there are multiple runs."""
        for run_list in self._data_context.current_runs:
            if run_list_to_string(run_list) == self._corrections_context.current_run_string:
                return run_list
        return None

    def set_current_run_string(self, run_string: str) -> None:
        """Sets the currently selected run string shown in the view."""
        self._corrections_context.current_run_string = run_string if run_string != "" else None

    def dead_times_average(self) -> float:
        """Returns the average dead time for the currently selected run."""
        dead_times = self._current_dead_times()
        return sum(dead_times) / len(dead_times) if len(dead_times) > 0 else 0.0

    def dead_times_range(self) -> tuple:
        """Returns the minimum and maximum dead time for the currently selected run."""
        dead_times = self._current_dead_times()
        return min(dead_times, default=0.0), max(dead_times, default=0.0)

    def _current_dead_times(self) -> list:
        """Returns a list of dead times for the currently displayed run and dead time mode."""
        table_name = self._current_dead_time_table_name()
        table = retrieve_ws(table_name) if table_name else None
        return table.toDict()[DEAD_TIME_TABLE_KEY] if table is not None else []

    def _current_dead_time_table_name(self) -> str:
        """Returns the name of the dead time table for the currently display run and dead time mode."""
        if self.is_dead_time_source_from_data_file():
            table_name = self._get_default_dead_time_table_for_run(self.current_runs())
        elif self.is_dead_time_source_from_workspace():
            table_name = self._corrections_context.dead_time_table_name
        elif self.is_dead_time_source_from_none():
            table_name = None
        return table_name

    def _get_default_dead_time_table_for_run(self, runs: list) -> str:
        """Returns the default table to use for dead time corrections for a specific run."""
        if runs is not None:
            run_data = self._data_context.get_loaded_data_for_run(runs)
            return run_data["DataDeadTimeTable"] if run_data is not None else None
        else:
            return None

    def set_dead_time_source_to_from_file(self) -> None:
        """Sets the dead time source to be 'FromFile'."""
        self._corrections_context.dead_time_source = DEAD_TIME_FROM_FILE
        self._corrections_context.dead_time_table_name = None

    def set_dead_time_source_to_from_ads(self, table_name: str) -> None:
        """Sets the dead time source to be 'FromADS'."""
        if table_name == "None":
            self.set_dead_time_source_to_none()
        else:
            self._corrections_context.dead_time_source = DEAD_TIME_FROM_ADS
            self._corrections_context.dead_time_table_name = table_name

    def set_dead_time_source_to_none(self) -> None:
        """Sets the dead time source to be 'None'."""
        self._corrections_context.dead_time_source = None
        self._corrections_context.dead_time_table_name = None

    def is_dead_time_source_from_data_file(self) -> bool:
        """Returns true if the dead time should be retrieved from a data file."""
        return self._corrections_context.dead_time_source == DEAD_TIME_FROM_FILE

    def is_dead_time_source_from_workspace(self) -> bool:
        """Returns true if the dead time should be retrieved from a workspace."""
        return self._corrections_context.dead_time_source == DEAD_TIME_FROM_ADS

    def is_dead_time_source_from_none(self) -> bool:
        """Returns true if the dead time should not be retrieved from any."""
        return self._corrections_context.dead_time_source is None

    def validate_selected_dead_time_workspace(self, table_name: str) -> str:
        """Validates the selected dead time workspace. Returns a string containing an error message if its invalid."""
        if check_if_workspace_exist(table_name):
            table = retrieve_ws(table_name)
            return self._validate_dead_time_table(table)
        else:
            return f"Workspace '{table_name}' does not exist in the ADS."

    def _validate_dead_time_table(self, table: ITableWorkspace) -> str:
        """Validates that the dead time table provided has the expected format."""
        if not isinstance(table, ITableWorkspace):
            return "The dead time table selected is not a Table Workspace."
        column_names = table.getColumnNames()
        number_of_rows = table.rowCount()
        if len(column_names) != 2:
            return f"Expected 2 columns, found {str(max(0, len(column_names)))} columns."
        if column_names[0] != "spectrum" or column_names[1] != DEAD_TIME_TABLE_KEY:
            return f"Columns have incorrect names."
        if number_of_rows != self._data_context.current_workspace.getNumberHistograms():
            return f"The number of histograms does not match the number of rows in dead time table ({number_of_rows})."
        return ""

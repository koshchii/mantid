# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2022 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +
import os
import systemtesting
import shutil

import mantid.simpleapi as mantid
from mantid import config

from isis_powder.osiris import Osiris

DIRS = config["datasearch.directories"].split(";")

# Setup various path details

inst_name = "OSIRIS"
# Relative to system data folder
working_folder_name = "ISIS_Powder"

# Relative to working folder
input_folder_name = "input"
output_folder_name = "output"

# Relative to input folder
calibration_folder_name = os.path.join("calibration", inst_name.lower())
calibration_map_rel_path = os.path.join("yaml_files", "osiris_system_test_mapping.yaml")

# Generate paths for the tests
# This implies DIRS[0] is the system test data folder
working_dir = os.path.join(DIRS[0], working_folder_name)

input_dir = os.path.join(working_dir, input_folder_name)
output_dir = os.path.join(working_dir, output_folder_name)

calibration_map_path = os.path.join(input_dir, calibration_map_rel_path)
calibration_dir = os.path.join(input_dir, calibration_folder_name)


def run_diffraction_focusing(
    sample_runs,
    user_name,
    output_file,
    merge_drange=False,
    subtract_empty_can=False,
    vanadium_normalisation=False,
):
    sample_runs = sample_runs  # Choose full drange set in the cycle 1_1

    osiris_inst_obj = setup_inst_object(user_name)

    # Run diffraction focusing
    osiris_inst_obj.run_diffraction_focusing(
        run_number=sample_runs,
        merge_drange=merge_drange,
        subtract_empty_can=subtract_empty_can,
        vanadium_normalisation=vanadium_normalisation,
    )

    focussed_rel_path = os.path.join("1_1", user_name, output_file)
    focused_path = os.path.join(output_dir, focussed_rel_path)

    foccussed_ws = mantid.Load(Filename=focused_path)

    return foccussed_ws


def setup_mantid_paths():
    config["datasearch.directories"] += ";" + input_dir


def setup_inst_object(user_name, with_container=False):
    inst_obj = Osiris(
        user_name=user_name,
        calibration_mapping_file=calibration_map_path,
        calibration_directory=calibration_dir,
        output_directory=output_dir,
    )

    return inst_obj


def _try_delete(path):
    try:
        # Use this instead of os.remove as we could be passed a non-empty dir
        if os.path.isdir(path):
            shutil.rmtree(path)
        else:
            os.remove(path)
    except OSError:
        print("Could not delete output file at: ", path)


class OSIRISDiffractionFocusingWithSubtractionTest(systemtesting.MantidSystemTest):
    calibration_results = None
    existing_config = config["datasearch.directories"]

    def requiredFiles(self):
        return self._gen_required_files()

    def runTest(self):
        setup_mantid_paths()
        self.calibration_results = run_diffraction_focusing(
            "119977",
            "Test",
            "OSI119977_d_spacing.nxs",
            subtract_empty_can=True,
        )

    def validate(self):
        foccussed_ws = self.calibration_results
        return (foccussed_ws.name(), "OSI119977_d_spacing.nxs")

    def cleanup(self):
        try:
            _try_delete(output_dir)
        finally:
            mantid.mtd.clear()
            config["datasearch.directories"] = self.existing_config

    def _gen_required_files(self):
        required_run_files = [
            "OSI82717.nxs",  # empty can
            "OSIRIS00119963.nxs",  # van
            "OSIRIS00119977.nxs",
        ]  # sample
        input_files = [os.path.join(input_dir, file) for file in required_run_files]
        input_files.append(calibration_map_path)
        return input_files


class OSIRISDiffractionFocusingWithMergingTest(systemtesting.MantidSystemTest):
    existing_config = config["datasearch.directories"]

    def requiredFiles(self):
        return self._gen_required_files()

    def runTest(self):
        setup_mantid_paths()
        self.calibration_results = run_diffraction_focusing(
            "119977-119978",
            "Test_Merge",
            "OSI119977-119978_d_spacing.nxs",
            merge_drange=True,
        )

    def validate(self):
        foccussed_ws = self.calibration_results
        return (foccussed_ws.name(), "OSI119977-119978_d_spacing.nxs")

    def cleanup(self):
        try:
            _try_delete(output_dir)
        finally:
            mantid.mtd.clear()
            config["datasearch.directories"] = self.existing_config

    def _gen_required_files(self):
        required_run_files = [
            "OSI82717.nxs",
            "OSI82718.nxs",  # empty can
            "OSIRIS00119963.nxs",
            "OSIRIS00119964.nxs",  # van
            "OSIRIS00119977.nxs",
            "OSIRIS00119978.nxs",
        ]  # sample
        input_files = [os.path.join(input_dir, file) for file in required_run_files]
        input_files.append(calibration_map_path)
        return input_files

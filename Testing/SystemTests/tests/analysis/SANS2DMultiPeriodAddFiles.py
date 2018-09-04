#pylint: disable=no-init
from __future__ import (absolute_import, division, print_function)
import stresstesting
from mantid.simpleapi import *
from mantid import config
import ISISCommandInterface as ii
import sans.command_interface.ISISCommandInterface as ii2


class SANS2DMultiPeriodAddFiles(stresstesting.MantidStressTest):

    def requiredMemoryMB(self):
        """Requires 2.5Gb"""
        return 2500

    def runTest(self):
        pass
        ii.SANS2D()
        ii.Set1D()
        ii.Detector("rear-detector")
        ii.MaskFile('MASKSANS2Doptions.091A')
        ii.Gravity(True)
        ii.add_runs( ('5512', '5512') ,'SANS2D', 'nxs', lowMem=True)

        #one period of a multi-period Nexus file
        ii.AssignSample('5512-add.nxs', period=7)

        ii.WavRangeReduction(2, 4, ii.DefaultTrans)
        paths = [os.path.join(config['defaultsave.directory'],'SANS2D00005512-add.nxs'),
                 os.path.join(config['defaultsave.directory'],'SANS2D00005512.log')]
        for path in paths:
            if os.path.exists(path):
                os.remove(path)

    def validate(self):
    # Need to disable checking of the Spectra-Detector map because it isn't
    # fully saved out to the nexus file (it's limited to the spectra that
    # are actually present in the saved workspace).
        self.disableChecking.append('SpectraMap')
        self.disableChecking.append('Instrument')
        self.disableChecking.append('Axes')

        return '5512p7rear_1D_2.0_4.0Phi-45.0_45.0','SANS2DMultiPeriodAddFiles.nxs'


class LARMORMultiPeriodAddEventFiles(stresstesting.MantidStressTest):
    def requiredMemoryMB(self):
        """Requires 2.5Gb"""
        return 2500

    def runTest(self):
        ii.LARMOR()
        ii.Set1D()
        ii.Detector("DetectorBench")
        ii.MaskFile('USER_LARMOR_151B_LarmorTeam_80tubes_BenchRot1p4_M4_r3699.txt')
        ii.Gravity(True)
        ii.add_runs( ('13065', '13065') ,'LARMOR', 'nxs', lowMem=True)

        ii.AssignSample('13065-add.nxs')
        ii.WavRangeReduction(2, 4, ii.DefaultTrans)

        # Clean up
        to_clean = ["13065_sans_nxs",
                    "13065p1rear_1D_2.0_4.0_incident_monitor",
                    "13065p2rear_1D_2.0_4.0_incident_monitor",
                    "13065p3rear_1D_2.0_4.0_incident_monitor",
                    "13065p4rear_1D_2.0_4.0_incident_monitor",
                    "80tubeCalibration_1-05-2015_r3157-3160"]
        for workspace in to_clean:
            DeleteWorkspace(workspace)

        paths = [os.path.join(config['defaultsave.directory'],'LARMOR00013065-add.nxs'),
                 os.path.join(config['defaultsave.directory'],'SANS2D00013065.log')]  # noqa
        for path in paths:
            if os.path.exists(path):
                os.remove(path)

    def validate(self):
    # Need to disable checking of the Spectra-Detector map because it isn't
    # fully saved out to the nexus file (it's limited to the spectra that
    # are actually present in the saved workspace).
        self.disableChecking.append('SpectraMap')
        self.disableChecking.append('Instrument')
        self.disableChecking.append('Axes')

        return "13065p1rear_1D_2.0_4.0" , "LARMORMultiPeriodAddEventFiles.nxs"


class SANS2DMultiPeriodAddFiles_V2(stresstesting.MantidStressTest):

    def requiredMemoryMB(self):
        """Requires 2.5Gb"""
        return 2500

    def runTest(self):
        ii2.UseCompatibilityMode()
        ii2.SANS2D()
        ii2.Set1D()
        ii2.Detector("rear-detector")
        ii2.MaskFile('MASKSANS2Doptions.091A')
        ii2.Gravity(True)
        ii2.AddRuns(('5512', '5512'), 'SANS2D', 'nxs', lowMem=True)

        # one period of a multi-period Nexus file
        ii2.AssignSample('5512-add.nxs', period=7)

        ii2.WavRangeReduction(2, 4, ii2.DefaultTrans)
        paths = [os.path.join(config['defaultsave.directory'], 'SANS2D00005512-add.nxs'),
                 os.path.join(config['defaultsave.directory'], 'SANS2D00005512.log')]
        for path in paths:
            if os.path.exists(path):
                os.remove(path)

    def validate(self):
        # Need to disable checking of the Spectra-Detector map because it isn't
        # fully saved out to the nexus file (it's limited to the spectra that
        # are actually present in the saved workspace).
        self.disableChecking.append('SpectraMap')
        self.disableChecking.append('Instrument')
        self.disableChecking.append('Axes')

        return '5512p7rear_1D_2.0_4.0Phi-45.0_45.0', 'SANS2DMultiPeriodAddFiles.nxs'


class LARMORMultiPeriodAddEventFilesTest_V2(stresstesting.MantidStressTest):
    def requiredMemoryMB(self):
        """Requires 2.5Gb"""
        return 2500

    def runTest(self):
        ii2.UseCompatibilityMode()
        ii2.LARMOR()
        ii2.Set1D()
        ii2.Detector("DetectorBench")
        ii2.MaskFile('USER_LARMOR_151B_LarmorTeam_80tubes_BenchRot1p4_M4_r3699.txt')
        ii2.Gravity(True)
        ii2.AddRuns(('13065', '13065'), 'LARMOR', 'nxs', lowMem=True)

        ii2.AssignSample('13065-add.nxs')
        ii2.WavRangeReduction(2, 4, ii2.DefaultTrans)

        # Clean up
        for element in AnalysisDataService.getObjectNames():
            if AnalysisDataService.doesExist(element) and element != "13065p1rear_1D_2.0_4.0":
                AnalysisDataService.remove(element)

        paths = [os.path.join(config['defaultsave.directory'], 'LARMOR00013065-add.nxs'),
                 os.path.join(config['defaultsave.directory'], 'SANS2D00013065.log')]  # noqa
        for path in paths:
            if os.path.exists(path):
                os.remove(path)

    def validate(self):
        # Need to disable checking of the Spectra-Detector map because it isn't
        # fully saved out to the nexus file (it's limited to the spectra that
        # are actually present in the saved workspace).
        self.disableChecking.append('SpectraMap')
        self.disableChecking.append('Instrument')
        self.disableChecking.append('Axes')

        return "13065p1rear_1D_2.0_4.0", "LARMORMultiPeriodAddEventFiles.nxs"

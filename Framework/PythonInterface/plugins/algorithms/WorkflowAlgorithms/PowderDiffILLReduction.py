from __future__ import (absolute_import, division, print_function)

import os
import numpy as np
from mantid.kernel import StringListValidator, Direction, FloatArrayProperty, \
    FloatArrayOrderedPairsValidator, VisibleWhenProperty, PropertyCriterion
from mantid.api import PythonAlgorithm, MultipleFileProperty, FileProperty, \
    FileAction, Progress, MatrixWorkspaceProperty
from mantid.simpleapi import *


class PowderDiffILLReduction(PythonAlgorithm):

    _calibration_file = None
    _roc_file = None
    _normalise_option = None
    _region_of_interest = []
    _observable = None
    _sort_x_axis = None
    _unit = None
    _out_name = None
    _progress = None
    _crop_negative = None
    _zero_counting_option = None
    _zero_cells = []
    _rebin_width = 0

    def _hide(self, name):
        return '__' + self._out_name + '_' + name

    def _hide_run(selfs, runnumber):
        return '__' + runnumber

    def category(self):
        return "ILL\\Diffraction;Diffraction\\Reduction"

    def summary(self):
        return 'Performs powder diffraction data reduction for ILL instrument D20.'

    def name(self):
        return "PowderDiffILLReduction"

    def validateInputs(self):
        issues = dict()
        return issues

    def PyInit(self):

        self.declareProperty(MultipleFileProperty('Run', extensions=['nxs']),
                             doc='File path of run(s).')

        self.declareProperty(FileProperty('CalibrationFile', '',
                                          action=FileAction.OptionalLoad, extensions=['nxs']),
                             doc='File containing the detector efficiencies.')

        self.declareProperty(FileProperty('ROCCorrectionFile', '',
                                          action=FileAction.OptionalLoad, extensions=['nxs']),
                             doc='File containing the radial oscillating collimator (ROC) corrections.')

        self.declareProperty(name='NormaliseTo',
                             defaultValue='None',
                             validator=StringListValidator(['None', 'Time', 'Monitor', 'ROI']),
                             doc='Normalise to time, monitor or ROI counts.')

        thetaRangeValidator = FloatArrayOrderedPairsValidator()

        self.declareProperty(FloatArrayProperty(name='ROI', values=[0, 153.6], validator=thetaRangeValidator),
                             doc='Regions of interest for normalisation [in scattering angle in degrees].')

        normaliseToROI = VisibleWhenProperty('NormaliseTo', PropertyCriterion.IsEqualTo, 'ROI')
        self.setPropertySettings('ROI', normaliseToROI)

        self.declareProperty(name='Observable',
                             defaultValue='sample.temperature',
                             doc='Scanning observable, a Sample Log entry.')

        self.declareProperty(name='SortObservableAxis',
                             defaultValue=False,
                             doc='Whether or not to sort the scanning observable axis.')

        self.declareProperty(name='RebinToWidth', defaultValue=0,
                             doc='Bin width to rebin the observable axis. Default is to not rebin.')

        self.declareProperty(name='CropNegative2Theta', defaultValue=True,
                             doc='Whether or not to crop out the bins corresponding to negative scattering angle.')

        self.declareProperty(name='ZeroCountingCells', defaultValue='Interpolate',
                             validator=StringListValidator(['Crop','Interpolate','Leave']),
                             doc='Crop out the zero counting cells or interpolate the counts from the neighbours.')

        self.declareProperty(name='Unit',
                             defaultValue='ScatteringAngle',
                             validator=StringListValidator(['ScatteringAngle', 'MomentumTransfer', 'dSpacing']),
                             doc='The unit of the reduced diffractogram.')

        self.declareProperty(MatrixWorkspaceProperty('OutputWorkspace', '',
                                                     direction=Direction.Output),
                             doc='Output workspace containing the reduced data.')

    def PyExec(self):

        self._configure()
        temp_ws = self._hide('temp')
        joined_ws = self._hide('joined')
        mon_ws = self._hide('mon')
        to_group = self._load()
        GroupWorkspaces(InputWorkspaces=to_group, OutputWorkspace=temp_ws)

        if self._normalise_option == 'Time':
            for ws in to_group:
                # normalise to time here, before joining, since the duration is in sample logs
                duration = mtd[ws].getRun().getLogData('duration').value
                Scale(InputWorkspace=ws,OutputWorkspace=ws,Factor=1./duration)

        try:
            ConjoinXRuns(InputWorkspaces=temp_ws, SampleLogAsXAxis=self._observable, OutputWorkspace=joined_ws)
        except RuntimeError:
            raise ValueError('Invalid scanning observable')

        DeleteWorkspace(temp_ws)

        ExtractMonitors(InputWorkspace=joined_ws, DetectorWorkspace=joined_ws, MonitorWorkspace=mon_ws)

        if self._normalise_option == 'Monitor':
            Divide(LHSWorkspace=joined_ws, RHSWorkspace=mon_ws, OutputWorkspace=joined_ws)
        elif self._normalise_option == 'ROI':
            self._normalise_to_roi(joined_ws)

        DeleteWorkspace(mon_ws)

        if self._calibration_file:
            calib_ws = self._hide('calib')
            LoadNexusProcessed(Filename=self._calibration_file, OutputWorkspace=calib_ws)
            Multiply(LHSWorkspace=joined_ws, RHSWorkspace=calib_ws, OutputWorkspace=joined_ws)
            DeleteWorkspace(calib_ws)

        if self._roc_file:
            roc_ws = self._hide('roc')
            LoadNexusProcessed(Filename=self._roc_file, OutputWorkspace=roc_ws)
            Multiply(LHSWorkspace=joined_ws, RHSWorkspace=roc_ws, OutputWorkspace=joined_ws)
            DeleteWorkspace(roc_ws)

        if self._sort_x_axis:
            SortXAxis(InputWorkspace=joined_ws, OutputWorkspace=joined_ws)

        theta_ws = self._hide('theta')
        ConvertSpectrumAxis(InputWorkspace=joined_ws, OutputWorkspace=theta_ws, Target='SignedTheta', OrderAxis=False)
        theta_axis = mtd[theta_ws].getAxis(1).extractValues()
        DeleteWorkspace(theta_ws)
        first_positive_theta = int(np.where(theta_axis > 0)[0][0])

        if self._crop_negative:
            self.log().information('First positive 2theta at workspace index: ' + str(first_positive_theta))
            CropWorkspace(InputWorkspace=joined_ws, OutputWorkspace=joined_ws, StartWorkspaceIndex=first_positive_theta)

        self._find_zero_cells(joined_ws)

        if self._zero_counting_option == 'Crop':
            self._crop_zero_cells(joined_ws, self._zero_cells)
        elif self._zero_counting_option == 'Interpolate':
            self._interpolate_zero_cells(joined_ws, theta_axis)

        target = 'SignedTheta'
        if self._unit == 'MomentumTransfer':
            target = 'ElasticQ'
        elif self._unit == 'dSpacing':
            target = 'ElasticDSpacing'

        ConvertSpectrumAxis(InputWorkspace=joined_ws, OutputWorkspace=joined_ws, Target=target)

        if self._rebin_width != 0:
            ConvertToHistogram(InputWorkspace=joined_ws, OutputWorkspace=joined_ws)
            x_axis = mtd[joined_ws].getAxis(0).extractValues()
            bin_widths = np.diff(x_axis)
            for spectrum in range(mtd[joined_ws].getNumberHistograms()):
                mtd[joined_ws].setY(spectrum, mtd[joined_ws].readY(spectrum) * bin_widths)
                mtd[joined_ws].setE(spectrum, mtd[joined_ws].readE(spectrum) * bin_widths)
            Rebin(InputWorkspace=joined_ws, OutputWorkspace=joined_ws, Params=self._rebin_width)
            Scale(InputWorkspace=joined_ws, Factor=1 / self._rebin_width)
            ConvertToPointData(InputWorkspace=joined_ws, OutputWorkspace=joined_ws)

        Transpose(InputWorkspace=joined_ws, OutputWorkspace=joined_ws)
        RenameWorkspace(InputWorkspace=joined_ws, OutputWorkspace=self._out_name)
        self.setProperty('OutputWorkspace', self._out_name)

    def _load(self):
        """
            Loads the list of runs
            If sum is requested, MergeRuns is called
            @return : the list of the loaded ws names
        """
        runs = self.getPropertyValue('Run')
        to_group = []
        self._progress = Progress(self, start=0.0, end=1.0, nreports=runs.count(',') + runs.count('+') + 1)

        for runs_list in runs.split(','):
            runs_sum = runs_list.split('+')
            if len(runs_sum) == 1:
                runnumber = os.path.basename(runs_sum[0]).split('.')[0]
                run = self._hide_run(runnumber)
                LoadILLDiffraction(Filename=runs_sum[0], OutputWorkspace=run)
                self._progress.report('Loaded run #' + runnumber)
                to_group.append(run)
            else:
                for i, run in enumerate(runs_sum):
                    runnumber = os.path.basename(run).split('.')[0]
                    if i == 0:
                        first = self._hide_run(runnumber + '_multiple')
                        LoadILLDiffraction(Filename=run, OutputWorkspace=first)
                        self._progress.report('Loaded run #' + runnumber)
                        to_group.append(first)
                    else:
                        run = self._hide_run(runnumber)
                        LoadILLDiffraction(Filename=run, OutputWorkspace=run)
                        self._progress.report('Loaded run #' + runnumber)
                        MergeRuns(InputWorkspaces=[first, run], OutputWorkspace=first)
                        DeleteWorkspace(Workspace=run)

        return to_group

    def _configure(self):
        """
            Configures the input properties
        """
        self._out_name = self.getPropertyValue('OutputWorkspace')
        self._observable = self.getPropertyValue('Observable')
        self._sort_x_axis = self.getProperty('SortObservableAxis').value
        self._normalise_option = self.getPropertyValue('NormaliseTo')
        self._calibration_file = self.getPropertyValue('CalibrationFile')
        self._roc_file = self.getPropertyValue('ROCCorrectionFile')
        self._unit = self.getPropertyValue('Unit')
        self._crop_negative = self.getProperty('CropNegative2Theta').value
        self._zero_counting_option = self.getPropertyValue('ZeroCountingCells')
        self._rebin_width = self.getProperty('RebinToWidth').value
        if self._normalise_option == 'ROI':
            self._region_of_interest = self.getProperty('ROI').value

    def _find_zero_cells(self, ws):
        """
            Finds the cells counting zeros
            @param ws: the input workspace
        """
        self._zero_cells = []
        size = mtd[ws].blocksize()
        for spectrum in range(mtd[ws].getNumberHistograms()):
            counts = mtd[ws].readY(spectrum)
            if np.count_nonzero(counts) < size/5:
                self._zero_cells.append(spectrum)
        self._zero_cells.sort()
        self.log().information('Found zero counting cells at indices: ' + str(self._zero_cells))

    def _crop_zero_cells(self, ws, wsIndexList):
        """
            Crops out the spectra corresponding to zero counting pixels
            @param ws: the input workspace
            @param wsIndexList: list of workspace indices to crop out
        """
        MaskDetectors(Workspace=ws, WorkspaceIndexList=wsIndexList)
        ExtractUnmaskedSpectra(InputWorkspace=ws, OutputWorkspace=ws)

    def _interpolate_zero_cells(self, ws, theta_axis):
        """
            Interpolates the counts of zero counting cells linearly from the
            nearest non-zero neighbour cells
            @param ws: the input workspace
            @param theta_axis: the unordered signed 2theta axis
        """
        unable_to_interpolate = []
        for cell in self._zero_cells:
            prev_cell = cell - 1
            next_cell = cell + 1

            while prev_cell in self._zero_cells:
                prev_cell-=1
            while next_cell in self._zero_cells:
                next_cell+=1

            if prev_cell == -1:
                self.log().notice('Unable to interpolate for cell #'+str(cell)+
                                  ': no non-zero neighbour cell was found on the left side. Bin will be cropped.')
                unable_to_interpolate.append(cell)
            if next_cell == mtd[ws].getNumberHistograms():
                self.log().notice('Unable to interpolate for cell #'+str(cell)+
                                  ': no non-zero neighbour cell was found on the right side. Bin will be cropped.')
                unable_to_interpolate.append(cell)

            if prev_cell >= 0 and next_cell < mtd[ws].getNumberHistograms():
                theta_prev = theta_axis[prev_cell]
                theta = theta_axis[cell]
                theta_next = theta_axis[next_cell]
                counts_prev = mtd[ws].readY(prev_cell)
                errors_prev = mtd[ws].readE(prev_cell)
                counts_next = mtd[ws].readY(next_cell)
                errors_next = mtd[ws].readE(next_cell)
                coefficient = (theta - theta_prev) / (theta_next - theta_prev)
                counts = counts_prev + coefficient * (counts_next - counts_prev)
                errors = errors_prev + coefficient * (errors_next - errors_prev)
                mtd[ws].setY(cell,counts)
                mtd[ws].setE(cell,errors)

        self._crop_zero_cells(ws, unable_to_interpolate)

    def _normalise_to_roi(self, ws):
        """
            Normalises counts to the sum of counts in the region-of-interest
            @param ws : input workspace with raw spectrum axis
        """
        roi_ws = self._hide('roi')
        theta_ws = self._hide('theta_ROI')
        ConvertSpectrumAxis(InputWorkspace=ws, OutputWorkspace=theta_ws, Target='SignedTheta')
        roi_pattern = self._parse_roi(theta_ws)
        SumSpectra(InputWorkspace=ws, OutputWorkspace=roi_ws, ListOfWorkspaceIndices=roi_pattern)
        SumSpectra(InputWorkspace=roi_ws, OutputWorkspace=roi_ws)
        Divide(LHSWorkspace=ws, RHSWorkspace=roi_ws, OutputWorkspace=ws)
        DeleteWorkspace(roi_ws)
        DeleteWorkspace(theta_ws)

    def _parse_roi(self, ws):
        """
            Parses the regions of interest string from 2theta ranges to workspace indices
            @param ws : input workspace with 2theta as spectrum axis
            Returns: roi as workspace indices, e.g. 7-20,100-123
        """
        result = ''
        axis = mtd[ws].getAxis(1).extractValues()
        index = 0
        while index < len(self._region_of_interest):
            start = self._region_of_interest[index]
            end = self._region_of_interest[index+1]
            start_index = np.argwhere(axis > start)
            end_index = np.argwhere(axis < end)
            result += str(start_index[0][0])+'-'+str(end_index[-1][0])
            result += ','
            index += 2
        self.log().information('ROI summing pattern is '+result[:-1])
        return result[:-1]

# Register the algorithm with Mantid
AlgorithmFactory.subscribe(PowderDiffILLReduction)

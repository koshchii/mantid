// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include "FitTab.h"
#include "Common/InterfaceUtils.h"
#include "Common/SettingsHelper.h"
#include "FitPlotView.h"

#include "MantidAPI/MultiDomainFunction.h"
#include "MantidQtWidgets/Common/WorkspaceUtils.h"

#include <QString>

using namespace Mantid::API;

namespace MantidQt {
namespace CustomInterfaces {
namespace Inelastic {

FitTab::FitTab(QWidget *parent, std::string const &tabName)
    : InelasticTab(parent), m_uiForm(new Ui::FitTab), m_dataPresenter(), m_fittingPresenter(), m_plotPresenter(),
      m_outOptionsPresenter(), m_fittingAlgorithm() {
  m_uiForm->setupUi(parent);
  parent->setWindowTitle(QString::fromStdString(tabName));
}

void FitTab::setup() {
  connect(m_uiForm->pbRun, SIGNAL(clicked()), this, SLOT(runTab()));
  updateResultOptions();
}

void FitTab::setupOutputOptionsPresenter(bool const editResults) {
  auto model = std::make_unique<FitOutputOptionsModel>();
  m_outOptionsPresenter =
      std::make_unique<FitOutputOptionsPresenter>(this, m_uiForm->ovOutputOptionsView, std::move(model));
  m_outOptionsPresenter->setEditResultVisible(editResults);
}

void FitTab::setupPlotView(std::optional<std::pair<double, double>> const &xPlotBounds) {
  auto model = std::make_unique<FitPlotModel>();
  m_plotPresenter = std::make_unique<FitPlotPresenter>(this, m_uiForm->dockArea->m_fitPlotView, std::move(model));
  m_plotPresenter->setFittingData(m_dataPresenter->getFittingData());
  m_plotPresenter->setFitOutput(m_fittingPresenter->getFitOutput());
  if (xPlotBounds) {
    m_plotPresenter->setXBounds(*xPlotBounds);
  }
  m_plotPresenter->updatePlots();
}

void FitTab::setModelFitFunction() {
  auto func = m_fittingPresenter->fitFunction();
  m_plotPresenter->setFitFunction(func);
  m_fittingPresenter->setFitFunction(func);
}

void FitTab::handleTableStartXChanged(double startX, WorkspaceID workspaceID, WorkspaceIndex spectrum) {
  if (m_plotPresenter->isCurrentlySelected(workspaceID, spectrum)) {
    m_plotPresenter->setStartX(startX);
    m_plotPresenter->updateGuess();
  }
}

void FitTab::handleTableEndXChanged(double endX, WorkspaceID workspaceID, WorkspaceIndex spectrum) {
  if (m_plotPresenter->isCurrentlySelected(workspaceID, spectrum)) {
    m_plotPresenter->setEndX(endX);
    m_plotPresenter->updateGuess();
  }
}

void FitTab::handleFunctionListChanged(const std::map<std::string, std::string> &functionStrings) {
  m_fittingPresenter->updateFunctionListInBrowser(functionStrings);
}

void FitTab::handleStartXChanged(double startX) {
  m_plotPresenter->setStartX(startX);
  m_dataPresenter->setStartX(startX, m_plotPresenter->getActiveWorkspaceID());
  updateParameterEstimationData();
  m_plotPresenter->updateGuess();
  m_dataPresenter->updateTableFromModel();
}

void FitTab::handleEndXChanged(double endX) {
  m_plotPresenter->setEndX(endX);
  m_dataPresenter->setEndX(endX, m_plotPresenter->getActiveWorkspaceID());
  updateParameterEstimationData();
  m_plotPresenter->updateGuess();
  m_dataPresenter->updateTableFromModel();
}

void FitTab::updateFitOutput(bool error) {
  disconnect(m_batchAlgoRunner, SIGNAL(batchComplete(bool)), this, SLOT(updateFitOutput(bool)));

  if (error) {
    m_fittingPresenter->cleanFailedRun(m_fittingAlgorithm);
    m_fittingAlgorithm.reset();
  } else {
    m_fittingPresenter->addOutput(m_fittingAlgorithm);
  }
}

void FitTab::updateSingleFitOutput(bool error) {
  disconnect(m_batchAlgoRunner, SIGNAL(batchComplete(bool)), this, SLOT(updateSingleFitOutput(bool)));

  if (error) {
    m_fittingPresenter->cleanFailedSingleRun(m_fittingAlgorithm, m_plotPresenter->getActiveWorkspaceID());
    m_fittingAlgorithm.reset();
  } else {
    m_fittingPresenter->addSingleFitOutput(m_fittingAlgorithm, m_plotPresenter->getActiveWorkspaceID(),
                                           m_plotPresenter->getActiveWorkspaceIndex());
  }
}

/**
 * Performs necessary state changes when the fit algorithm was run
 * and completed within this interface.
 */
void FitTab::fitAlgorithmComplete(bool error) {
  m_plotPresenter->setFitSingleSpectrumIsFitting(false);
  enableFitButtons(true);
  enableOutputOptions(!error);
  m_fittingPresenter->setErrorsEnabled(!error);
  if (!error) {
    m_fittingPresenter->updateFitBrowserParameterValuesFromAlg(m_fittingAlgorithm,
                                                               m_dataPresenter->getNumberOfDomains());
    m_fittingPresenter->updateFitStatus(m_fittingAlgorithm, m_dataPresenter->getNumberOfDomains());
    setModelFitFunction();
  }
  m_plotPresenter->updatePlots();
  disconnect(m_batchAlgoRunner, SIGNAL(batchComplete(bool)), this, SLOT(fitAlgorithmComplete(bool)));
}

/**
 * Plots the spectra corresponding to the selected parameters
 */
void FitTab::handlePlotSelectedSpectra() {
  enableFitButtons(false);
  for (auto const &spectrum : m_outOptionsPresenter->getSpectraToPlot())
    m_plotter->plotSpectra(spectrum.first, std::to_string(spectrum.second), SettingsHelper::externalPlotErrorBars());
  m_outOptionsPresenter->clearSpectraToPlot();
  enableFitButtons(true);
  m_outOptionsPresenter->setPlotting(false);
}

void FitTab::handleSingleFitClicked(WorkspaceID workspaceID, WorkspaceIndex spectrum) {
  if (validate()) {
    m_plotPresenter->setFitSingleSpectrumIsFitting(true);
    enableFitButtons(false);
    enableOutputOptions(false);
    m_fittingPresenter->setFittingMode(FittingMode::SIMULTANEOUS);
    runSingleFit(m_fittingPresenter->getSingleFit(workspaceID, spectrum));
  }
}

bool FitTab::validate() {
  UserInputValidator validator;
  m_dataPresenter->validate(validator);
  m_fittingPresenter->validate(validator);

  const auto error = validator.generateErrorMessage().toStdString();
  if (!error.empty()) {
    displayWarning(error);
  }
  return error.empty();
}

/**
 * Called when the 'Run' button is called in the InelasticTab.
 */
void FitTab::run() {
  enableFitButtons(false);
  enableOutputOptions(false);
  m_fittingPresenter->updateFittingModeFromBrowser();
  runFitAlgorithm(m_fittingPresenter->getFittingAlgorithm(m_fittingPresenter->getFittingMode()));
}

/**
 * Enables or disables the 'Run', 'Fit Single Spectrum' and other related
 * buttons
 * @param enable :: true to enable buttons
 */
void FitTab::enableFitButtons(bool enable) {
  m_uiForm->pbRun->setText(enable ? "Run" : "Running...");
  m_uiForm->pbRun->setEnabled(enable);
  m_plotPresenter->setFitSingleSpectrumEnabled(enable);
  m_fittingPresenter->setFitEnabled(enable);
}

/**
 * Enables or disables the output options. It also sets the current result and
 * PDF workspaces to be plotted
 * @param enable :: true to enable buttons
 */
void FitTab::enableOutputOptions(bool enable) {
  if (enable) {
    m_outOptionsPresenter->setResultWorkspace(m_fittingPresenter->getResultWorkspace());
    setPDFWorkspace(m_fittingPresenter->getOutputBasename() + "_PDFs");
    m_outOptionsPresenter->setPlotTypes("Result Group");
  } else
    m_outOptionsPresenter->setMultiWorkspaceOptionsVisible(enable);

  m_outOptionsPresenter->setPlotEnabled(enable && m_outOptionsPresenter->isSelectedGroupPlottable());
  m_outOptionsPresenter->setEditResultEnabled(enable);
  m_outOptionsPresenter->setSaveEnabled(enable);
}

/**
 * Sets the active PDF workspace within the output options if one exists for the
 * current run
 * @param workspaceName :: the name of the PDF workspace if it exists
 */
void FitTab::setPDFWorkspace(std::string const &workspaceName) {
  auto const fabMinimizer = m_fittingPresenter->minimizer() == "FABADA";
  auto const enablePDFOptions = WorkspaceUtils::doesExistInADS(workspaceName) && fabMinimizer;

  if (enablePDFOptions) {
    m_outOptionsPresenter->setPDFWorkspace(WorkspaceUtils::getADSWorkspace<WorkspaceGroup>(workspaceName));
    m_outOptionsPresenter->setPlotWorkspaces();
  } else
    m_outOptionsPresenter->removePDFWorkspace();
  m_outOptionsPresenter->setMultiWorkspaceOptionsVisible(enablePDFOptions);
}

void FitTab::updateParameterEstimationData() {
  m_fittingPresenter->updateParameterEstimationData(
      m_dataPresenter->getDataForParameterEstimation(m_fittingPresenter->getEstimationDataSelector()));
  m_fittingPresenter->estimateFunctionParameters(m_plotPresenter->getActiveWorkspaceID(),
                                                 m_plotPresenter->getActiveWorkspaceIndex());
}

/*
 * Runs the specified fit algorithm and calls the algorithmComplete
 * method of this fit analysis tab once completed.
 *
 * @param fitAlgorithm      The fit algorithm to run.
 */
void FitTab::runFitAlgorithm(IAlgorithm_sptr fitAlgorithm) {
  connect(m_batchAlgoRunner, SIGNAL(batchComplete(bool)), this, SLOT(updateFitOutput(bool)));
  setupFit(std::move(fitAlgorithm));
  m_batchAlgoRunner->executeBatchAsync();
}

void FitTab::runSingleFit(IAlgorithm_sptr fitAlgorithm) {
  connect(m_batchAlgoRunner, SIGNAL(batchComplete(bool)), this, SLOT(updateSingleFitOutput(bool)));
  setupFit(std::move(fitAlgorithm));
  m_batchAlgoRunner->executeBatchAsync();
}

void FitTab::setupFit(IAlgorithm_sptr fitAlgorithm) {
  auto properties = m_fittingPresenter->fitProperties();
  m_fittingAlgorithm = fitAlgorithm;
  m_batchAlgoRunner->addAlgorithm(fitAlgorithm, std::move(properties));
  connect(m_batchAlgoRunner, SIGNAL(batchComplete(bool)), this, SLOT(fitAlgorithmComplete(bool)));
}

void FitTab::updateDataReferences() {
  m_fittingPresenter->updateFunctionBrowserData(static_cast<int>(m_dataPresenter->getNumberOfDomains()),
                                                m_dataPresenter->getDatasets(), m_dataPresenter->getQValuesForData(),
                                                m_dataPresenter->getResolutionsForFit());
  setModelFitFunction();
}

/**
 * Updates whether the options for plotting and saving fit results are
 * enabled/disabled.
 */
void FitTab::updateResultOptions() {
  const bool isFit = m_fittingPresenter->isPreviouslyFit(m_plotPresenter->getActiveWorkspaceID(),
                                                         m_plotPresenter->getActiveWorkspaceIndex());
  if (isFit)
    m_outOptionsPresenter->setResultWorkspace(m_fittingPresenter->getResultWorkspace());
  m_outOptionsPresenter->setPlotEnabled(isFit);
  m_outOptionsPresenter->setEditResultEnabled(isFit);
  m_outOptionsPresenter->setSaveEnabled(isFit);
}

std::string FitTab::tabName() const { return m_parentWidget->windowTitle().toStdString(); }

void FitTab::handleDataChanged() {
  updateDataReferences();
  m_fittingPresenter->removeFittingData();
  m_plotPresenter->updateAvailableSpectra();
  m_plotPresenter->updatePlots();
  m_plotPresenter->updateGuessAvailability();
  updateParameterEstimationData();
  updateResultOptions();
}

void FitTab::handleDataAdded(IAddWorkspaceDialog const *dialog) {
  if (m_dataPresenter->addWorkspaceFromDialog(dialog)) {
    m_fittingPresenter->addDefaultParameters();
  }
  updateDataReferences();
  m_plotPresenter->appendLastDataToSelection(m_dataPresenter->createDisplayNames());
  updateParameterEstimationData();
}

void FitTab::handleDataRemoved() {
  m_fittingPresenter->removeDefaultParameters();
  updateDataReferences();
  m_plotPresenter->updateDataSelection(m_dataPresenter->createDisplayNames());
  updateParameterEstimationData();
}

void FitTab::handlePlotSpectrumChanged() {
  auto const index = m_plotPresenter->getSelectedDomainIndex();
  m_fittingPresenter->setCurrentDataset(index);
}

void FitTab::handleFwhmChanged(double fwhm) {
  m_fittingPresenter->setFWHM(m_plotPresenter->getActiveWorkspaceID(), fwhm);
  m_fittingPresenter->updateFitBrowserParameterValues();
  m_plotPresenter->updateGuess();
}

void FitTab::handleBackgroundChanged(double value) {
  m_fittingPresenter->setBackground(m_plotPresenter->getActiveWorkspaceID(), value);
  setModelFitFunction();
  m_plotPresenter->updateGuess();
}

void FitTab::handleFunctionChanged() {
  setModelFitFunction();
  m_fittingPresenter->removeFittingData();
  m_plotPresenter->updatePlots();
  m_plotPresenter->updateFit();
  m_fittingPresenter->updateFitTypeString();
}

} // namespace Inelastic
} // namespace CustomInterfaces
} // namespace MantidQt

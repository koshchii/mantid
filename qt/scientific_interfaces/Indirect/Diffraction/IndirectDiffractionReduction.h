// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "DllConfig.h"

#include "Common/IndirectInterface.h"
#include "Common/OutputPlotOptionsPresenter.h"

#include "ui_IndirectDiffractionReduction.h"

#include "MantidQtWidgets/Common/BatchAlgorithmRunner.h"

namespace MantidQt {
namespace CustomInterfaces {

class DetectorGroupingOptions;

class IndirectDiffractionReduction : public IndirectInterface {
  Q_OBJECT

public:
  /// Default Constructor
  explicit IndirectDiffractionReduction(QWidget *parent = nullptr);
  ~IndirectDiffractionReduction() override;

  /// The name of the interface as registered into the factory
  static std::string name() { return "Diffraction"; }
  /// This interface's categories.
  static QString categoryInfo() { return "Indirect"; }

public slots:
  void instrumentSelected(const QString &instrumentName, const QString &analyserName, const QString &reflectionName);
  void run();
  void saveReductions();
  void runFilesChanged();
  void runFilesFinding();
  void runFilesFound();
  void algorithmComplete(bool error);
  void validateSpectrumMin(int value);
  void validateSpectrumMax(int value);

private:
  std::string documentationPage() const override;

  void initLayout() override;
  void initLocalPython() override;

  void loadSettings();
  void saveSettings();

  Mantid::API::IAlgorithm_sptr saveGSSAlgorithm(const std::string &filename);
  Mantid::API::IAlgorithm_sptr saveASCIIAlgorithm(const std::string &filename, const std::string &inputWsName);
  Mantid::API::IAlgorithm_sptr saveNexusProcessedAlgorithm(const std::string &filename, const std::string &inputWsName);
  Mantid::API::IAlgorithm_sptr saveAlgorithm(const std::string &saveAlgName, const std::string &filename,
                                             const std::string &inputWsName = "", const int &version = -1);
  Mantid::API::IAlgorithm_sptr convertUnitsAlgorithm(const std::string &inputWsName, const std::string &outputWsName,
                                                     const std::string &target);

  bool validateRebin();
  bool validateVanCal();
  bool validateCalOnly();

  Mantid::API::MatrixWorkspace_sptr loadInstrument(const std::string &instrumentName,
                                                   const std::string &reflection = "");

  void runGenericReduction(const QString &instName, const QString &mode);
  void connectRunButtonValidation(const MantidQt::API::FileFinderWidget *file_field);
  void runOSIRISdiffonlyReduction();

  void setRunIsRunning(bool running);
  void setButtonsEnabled(bool enabled);
  void setRunEnabled(bool enabled);
  void setSaveEnabled(bool enabled);

private:
  /// The settings dialog
  Ui::IndirectDiffractionReduction m_uiForm; /// The form generated using Qt Designer
  QDoubleValidator *m_valDbl;
  QString m_settingsGroup; /// The settings group
  MantidQt::API::BatchAlgorithmRunner *m_batchAlgoRunner;
  std::vector<std::string> m_plotWorkspaces;

  std::unique_ptr<OutputPlotOptionsPresenter> m_plotOptionsPresenter;
  DetectorGroupingOptions *m_groupingWidget;
};
} // namespace CustomInterfaces
} // namespace MantidQt

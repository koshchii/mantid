// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2024 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "DllConfig.h"
#include "MantidAPI/MatrixWorkspace_fwd.h"

#include <QStringList>

#include <memory>
#include <string>
#include <tuple>

namespace MantidQt {
namespace CustomInterfaces {

class IndirectPlotOptionsView;
class ISqwPresenter;

class MANTIDQT_INELASTIC_DLL ISqwView {

public:
  virtual void subscribePresenter(ISqwPresenter *presenter) = 0;

  virtual IndirectPlotOptionsView *getPlotOptions() const = 0;
  virtual void setFBSuffixes(QStringList const &suffix) = 0;
  virtual void setWSSuffixes(QStringList const &suffix) = 0;
  virtual std::tuple<double, double> getQRangeFromPlot() const = 0;
  virtual std::tuple<double, double> getERangeFromPlot() const = 0;
  virtual std::string getDataName() const = 0;
  virtual void plotRqwContour(Mantid::API::MatrixWorkspace_sptr rqwWorkspace) = 0;
  virtual void setDefaultQAndEnergy() = 0;
  virtual void setSaveEnabled(bool const enabled) = 0;
  virtual bool validate() = 0;
  virtual void showMessageBox(std::string const &message) const = 0;
  virtual void updateRunButton(bool const enabled, std::string const &enableOutputButtons = "unchanged",
                               std::string const &message = "Run", std::string const &tooltip = "") = 0;
};
} // namespace CustomInterfaces
} // namespace MantidQt
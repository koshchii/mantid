// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//     NScD Oak Ridge National Laboratory, European Spallation Source
//     & Institut Laue - Langevin
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidAlgorithms/MaskBinsIf.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidAPI/NumericAxis.h"
#include "MantidAPI/SpectraAxis.h"
#include "MantidHistogramData/HistogramIterator.h"
#include "MantidKernel/MultiThreaded.h"

#include <muParser.h>

namespace {

/**
 * @brief Makes a new instance of the muparser
 * @param y : count
 * @param e : error
 * @param x : bin center
 * @param dx : bin center error
 * @param s : spectrum axis value
 * @param i : workspace index
 * @param criterion : expression
 * @return muparser
 */
mu::Parser makeParser(double &y, double &e, double &x, double &dx, double &s,
                      double &i, const std::string &criterion) {
  mu::Parser muParser;
  muParser.DefineVar("y", &y);
  muParser.DefineVar("e", &e);
  muParser.DefineVar("x", &x);
  muParser.DefineVar("dx", &dx);
  muParser.DefineVar("s", &s);
  muParser.DefineVar("i", &i);
  muParser.SetExpr(criterion);
  return muParser;
}
}

namespace Mantid {
namespace Algorithms {

using namespace API;
using namespace Kernel;

// Register the algorithm into the AlgorithmFactory
DECLARE_ALGORITHM(MaskBinsIf)

//----------------------------------------------------------------------------------------------
/** Initialize the algorithm's properties.
 */
void MaskBinsIf::init() {
  declareProperty(Kernel::make_unique<WorkspaceProperty<MatrixWorkspace>>(
                      "InputWorkspace", "", Direction::Input),
                  "An input workspace.");
  declareProperty("Criterion", "",
                  "Masking criterion as a muparser expression; y: bin count, "
                  "e: bin error, x: bin center, dx: bin center error, s: "
                  "spectrum axis value, i: workspace index.");
  declareProperty(Kernel::make_unique<WorkspaceProperty<MatrixWorkspace>>(
                      "OutputWorkspace", "", Direction::Output),
                  "An output workspace.");
}

//----------------------------------------------------------------------------------------------
/** Validate the inputs.
 */
std::map<std::string, std::string> MaskBinsIf::validateInputs() {
  std::map<std::string, std::string> issues;
  double y = 0., e = 0., x = 0., dx = 0., s = 0., i = 0.;
  mu::Parser parser =
      makeParser(y, e, x, dx, s, i, getPropertyValue("Criterion"));
  try {
    parser.Eval();
  } catch (mu::Parser::exception_type &e) {
    issues["Criterion"] = "Invalid expression given: " + e.GetMsg();
  }
  return issues;
}

//----------------------------------------------------------------------------------------------
/** Execute the algorithm.
 */
void MaskBinsIf::exec() {
  const std::string criterion = getPropertyValue("Criterion");
  MatrixWorkspace_const_sptr inputWorkspace = getProperty("InputWorkspace");
  MatrixWorkspace_sptr outputWorkspace = getProperty("OutputWorkspace");
  if (inputWorkspace != outputWorkspace) {
    outputWorkspace = inputWorkspace->clone();
  }
  const auto spectrumAxis = outputWorkspace->getAxis(1);
  const auto numeric = dynamic_cast<NumericAxis *>(spectrumAxis);
  const auto spectrum = dynamic_cast<SpectraAxis *>(spectrumAxis);
  const bool spectrumOrNumeric = numeric || spectrum;
  PARALLEL_FOR_IF(Mantid::Kernel::threadSafe(*outputWorkspace))
  for (size_t index = 0; index < outputWorkspace->getNumberHistograms();
       ++index) {
    PARALLEL_START_INTERUPT_REGION
    double y, e, x, dx, s, i;
    dx = 0.;
    s = 0.;
    i = index;
    mu::Parser parser = makeParser(y, e, x, dx, s, i, criterion);
    if (spectrumOrNumeric) {
      s = spectrumAxis->getValue(index);
    }
    const auto &spectrum = outputWorkspace->histogram(index);
    const bool hasDx = outputWorkspace->hasDx(index);
    for (auto it = spectrum.begin(); it != spectrum.end(); ++it) {
      const auto bin = std::distance(spectrum.begin(), it);
      y = it->counts();
      x = it->center();
      e = it->countStandardDeviation();
      if (hasDx) {
        dx = it->centerError();
      }
      if (parser.Eval() != 0.) {
        outputWorkspace->maskBin(index, bin);
      }
    }
    PARALLEL_END_INTERUPT_REGION
  }
  PARALLEL_CHECK_INTERUPT_REGION
  setProperty("OutputWorkspace", outputWorkspace);
}

} // namespace Algorithms
} // namespace Mantid

#include "MantidDataHandling/JoinISISPolarizationEfficiencies.h"

#include "MantidAPI/FileProperty.h"
#include "MantidAPI/TextAxis.h"
#include "MantidAPI/WorkspaceFactory.h"
#include "MantidDataObjects/Workspace2D.h"
#include "MantidDataObjects/WorkspaceCreation.h"
#include "MantidHistogramData/Histogram.h"
#include "MantidHistogramData/Interpolate.h"
#include "MantidHistogramData/LinearGenerator.h"
#include "MantidKernel/make_unique.h"
#include <limits>

namespace Prop {

const std::string PP_WORKSPACE("Pp");
const std::string AP_WORKSPACE("Ap");
const std::string RHO_WORKSPACE("Rho");
const std::string ALPHA_WORKSPACE("Alpha");

const std::string P1_WORKSPACE("P1");
const std::string P2_WORKSPACE("P2");
const std::string F1_WORKSPACE("F1");
const std::string F2_WORKSPACE("F2");

const std::string OUTPUT_WORKSPACE("OutputWorkspace");

} // namespace Prop

using namespace Mantid::API;
using namespace Mantid::DataObjects;
using namespace Mantid::HistogramData;
using namespace Mantid::Kernel;

namespace Mantid {
namespace DataHandling {

// Register the algorithm into the AlgorithmFactory
DECLARE_ALGORITHM(JoinISISPolarizationEfficiencies)

//----------------------------------------------------------------------------------------------

/// Algorithms name for identification. @see Algorithm::name
const std::string JoinISISPolarizationEfficiencies::name() const {
  return "JoinISISPolarizationEfficiencies";
}

/// Algorithm's version for identification. @see Algorithm::version
int JoinISISPolarizationEfficiencies::version() const { return 1; }

/// Algorithm's category for identification. @see Algorithm::category
const std::string JoinISISPolarizationEfficiencies::category() const {
  return "DataHandling;ISIS\\Reflectometry";
}

/// Algorithm's summary for use in the GUI and help. @see Algorithm::summary
const std::string JoinISISPolarizationEfficiencies::summary() const {
  return "Joins workspaces containing ISIS reflectometry polarization "
         "efficiency factors into a single workspace ready to be used with "
         "PolarizationEfficiencyCor.";
}

//----------------------------------------------------------------------------------------------
/** Initialize the algorithm's properties.
 */
void JoinISISPolarizationEfficiencies::init() {

  declareProperty(
      Kernel::make_unique<WorkspaceProperty<MatrixWorkspace>>(
          Prop::PP_WORKSPACE, "", Kernel::Direction::Input,
          PropertyMode::Optional),
      "A matrix workspaces containing the Pp polarization efficiency.");

  declareProperty(
      Kernel::make_unique<WorkspaceProperty<MatrixWorkspace>>(
          Prop::AP_WORKSPACE, "", Kernel::Direction::Input,
          PropertyMode::Optional),
      "A matrix workspaces containing the Ap polarization efficiency.");

  declareProperty(
      Kernel::make_unique<WorkspaceProperty<MatrixWorkspace>>(
          Prop::RHO_WORKSPACE, "", Kernel::Direction::Input,
          PropertyMode::Optional),
      "A matrix workspaces containing the Rho polarization efficiency.");

  declareProperty(
      Kernel::make_unique<WorkspaceProperty<MatrixWorkspace>>(
          Prop::ALPHA_WORKSPACE, "", Kernel::Direction::Input,
          PropertyMode::Optional),
      "A matrix workspaces containing the Alpha polarization efficiency.");

  declareProperty(
      Kernel::make_unique<WorkspaceProperty<MatrixWorkspace>>(
          Prop::P1_WORKSPACE, "", Kernel::Direction::Input,
          PropertyMode::Optional),
      "A matrix workspaces containing the P1 polarization efficiency.");

  declareProperty(
      Kernel::make_unique<WorkspaceProperty<MatrixWorkspace>>(
          Prop::P2_WORKSPACE, "", Kernel::Direction::Input,
          PropertyMode::Optional),
      "A matrix workspaces containing the P2 polarization efficiency.");

  declareProperty(
      Kernel::make_unique<WorkspaceProperty<MatrixWorkspace>>(
          Prop::F1_WORKSPACE, "", Kernel::Direction::Input,
          PropertyMode::Optional),
      "A matrix workspaces containing the F1 polarization efficiency.");

  declareProperty(
      Kernel::make_unique<WorkspaceProperty<MatrixWorkspace>>(
          Prop::F2_WORKSPACE, "", Kernel::Direction::Input,
          PropertyMode::Optional),
      "A matrix workspaces containing the F2 polarization efficiency.");

  declareProperty(Kernel::make_unique<WorkspaceProperty<API::MatrixWorkspace>>(
                      Prop::OUTPUT_WORKSPACE, "", Direction::Output),
                  "An output workspace containing the efficiencies.");
}

//----------------------------------------------------------------------------------------------
/** Execute the algorithm.
 */
void JoinISISPolarizationEfficiencies::exec() {
  auto const propsFredrikze =
      getNonDefaultProperties({Prop::PP_WORKSPACE, Prop::AP_WORKSPACE,
                               Prop::RHO_WORKSPACE, Prop::ALPHA_WORKSPACE});
  auto const propsWildes =
      getNonDefaultProperties({Prop::P1_WORKSPACE, Prop::P2_WORKSPACE,
                               Prop::F1_WORKSPACE, Prop::F2_WORKSPACE});

  if (propsFredrikze.empty() && propsWildes.empty()) {
    throw std::invalid_argument(
        "At least one of the efficiency file names must be set.");
  }

  if (!propsFredrikze.empty() && !propsWildes.empty()) {
    throw std::invalid_argument(
        "Efficiencies belonging to different methods cannot mix.");
  }

  MatrixWorkspace_sptr efficiencies;
  if (!propsFredrikze.empty()) {
    efficiencies = createEfficiencies(propsFredrikze);
  } else {
    efficiencies = createEfficiencies(propsWildes);
  }

  setProperty("OutputWorkspace", efficiencies);
}

/// Get names of non-default properties out of a list of names
/// @param labels :: Names of properties to check.
std::vector<std::string>
JoinISISPolarizationEfficiencies::getNonDefaultProperties(
    std::vector<std::string> const &labels) const {
  std::vector<std::string> outputLabels;
  for (auto const &label : labels) {
    if (!isDefault(label)) {
      outputLabels.push_back(label);
    }
  }
  return outputLabels;
}

/// Load efficientcies from files and put them into a single workspace.
/// @param props :: Names of properties containg names of files to load.
MatrixWorkspace_sptr JoinISISPolarizationEfficiencies::createEfficiencies(
    std::vector<std::string> const &props) {
  std::vector<MatrixWorkspace_sptr> workspaces;
  for (auto const &propName : props) {
    MatrixWorkspace_sptr ws = getProperty(propName);
    if (ws->getNumberHistograms() != 1) {
      throw std::runtime_error(
          "Loaded workspace must contain a single histogram. Found " +
          std::to_string(ws->getNumberHistograms()));
    }
    workspaces.push_back(ws);
  }

  return createEfficiencies(props, workspaces);
}

/// Create the efficiency workspace by combining single spectra workspaces into
/// one.
/// @param labels :: Axis labels for each workspace.
/// @param workspaces :: Workspaces to put together.
MatrixWorkspace_sptr JoinISISPolarizationEfficiencies::createEfficiencies(
    std::vector<std::string> const &labels,
    std::vector<MatrixWorkspace_sptr> const &workspaces) {
  auto interpolatedWorkspaces = interpolateWorkspaces(workspaces);

  auto const &inWS = interpolatedWorkspaces.front();
  MatrixWorkspace_sptr outWS = WorkspaceFactory::Instance().create(
      inWS, labels.size(), inWS->x(0).size(), inWS->blocksize());
  auto axis1 = new TextAxis(labels.size());
  outWS->replaceAxis(1, axis1);
  outWS->getAxis(0)->setUnit("Wavelength");

  for (size_t i = 0; i < interpolatedWorkspaces.size(); ++i) {
    auto &ws = interpolatedWorkspaces[i];
    outWS->mutableX(i) = ws->x(0);
    outWS->mutableY(i) = ws->y(0);
    outWS->mutableE(i) = ws->e(0);
    axis1->setLabel(i, labels[i]);
  }

  return outWS;
}

/// Interpolate the workspaces so that all have the same blocksize.
/// @param workspaces :: The workspaces to interpolate.
/// @return A list of interpolated workspaces.
std::vector<MatrixWorkspace_sptr>
JoinISISPolarizationEfficiencies::interpolateWorkspaces(
    std::vector<MatrixWorkspace_sptr> const &workspaces) {
  size_t minSize(std::numeric_limits<size_t>::max());
  size_t maxSize(0);
  bool thereAreHistograms = false;
  bool allAreHistograms = true;

  // Find out if the workspaces need to be interpolated.
  for (auto const &ws : workspaces) {
    auto size = ws->blocksize();
    if (size < minSize) {
      minSize = size;
    }
    if (size > maxSize) {
      maxSize = size;
    }
    thereAreHistograms = thereAreHistograms || ws->isHistogramData();
    allAreHistograms = allAreHistograms && ws->isHistogramData();
  }

  if (thereAreHistograms != allAreHistograms) {
    throw std::invalid_argument("Cannot mix histograms and point data.");
  }

  // All same size, same type - nothing to do
  if (minSize == maxSize) {
    return workspaces;
  }

  // Interpolate those that need interpolating
  std::vector<MatrixWorkspace_sptr> interpolatedWorkspaces;
  for (auto const &ws : workspaces) {
    if (ws->blocksize() < maxSize) {
      if (allAreHistograms) {
        interpolatedWorkspaces.push_back(
            interpolateHistogramWorkspace(ws, maxSize));
      } else {
        interpolatedWorkspaces.push_back(
            interpolatePointDataWorkspace(ws, maxSize));
      }
    } else {
      interpolatedWorkspaces.push_back(ws);
    }
  }

  return interpolatedWorkspaces;
}

MatrixWorkspace_sptr
JoinISISPolarizationEfficiencies::interpolatePointDataWorkspace(
    MatrixWorkspace_sptr ws, size_t const maxSize) {
  auto const &x = ws->x(0);
  auto const startX = x.front();
  auto const endX = x.back();
  Counts yVals(maxSize, 0.0);
  auto const dX = (endX - startX) / double(maxSize - 1);
  Points xVals(maxSize, LinearGenerator(startX, dX));
  auto newHisto = Histogram(xVals, yVals);
  interpolateLinearInplace(ws->histogram(0), newHisto);
  auto interpolatedWS = boost::make_shared<Workspace2D>();
  interpolatedWS->initialize(1, newHisto);
  assert(interpolatedWS->y(0).size() == maxSize);
  return interpolatedWS;
}

MatrixWorkspace_sptr
JoinISISPolarizationEfficiencies::interpolateHistogramWorkspace(
    MatrixWorkspace_sptr ws, size_t const maxSize) {
  ws->setDistribution(true);
  auto const &x = ws->x(0);
  auto const dX = (x.back() - x.front()) / double(maxSize);
  std::vector<double> params(2 * maxSize + 1);
  for (size_t i = 0; i < maxSize; ++i) {
    params[2 * i] = x.front() + dX * double(i);
    params[2 * i + 1] = dX;
  }
  params.back() = x.back();
  auto alg = createChildAlgorithm("InterpolatingRebin");
  alg->setProperty("InputWorkspace", ws);
  alg->setProperty("Params", params);
  alg->setProperty("OutputWorkspace", "dummy");
  alg->execute();
  MatrixWorkspace_sptr interpolatedWS = alg->getProperty("OutputWorkspace");
  assert(interpolatedWS->y(0).size() == maxSize);
  assert(interpolatedWS->x(0).size() == maxSize + 1);
  return interpolatedWS;
}

} // namespace DataHandling
} // namespace Mantid

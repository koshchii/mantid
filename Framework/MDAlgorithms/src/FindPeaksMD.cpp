// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidMDAlgorithms/FindPeaksMD.h"
#include "MantidAPI/Run.h"
#include "MantidDataObjects/MDEventFactory.h"
#include "MantidDataObjects/MDHistoWorkspace.h"
#include "MantidGeometry/Crystal/EdgePixel.h"
#include "MantidGeometry/Instrument/Goniometer.h"
#include "MantidGeometry/Objects/InstrumentRayTracer.h"
#include "MantidKernel/BoundedValidator.h"
#include "MantidKernel/EnabledWhenProperty.h"
#include "MantidKernel/ListValidator.h"
#include "MantidKernel/VMD.h"

#include <map>
#include <vector>

using namespace Mantid::Kernel;
using namespace Mantid::API;
using namespace Mantid::DataObjects;
using namespace Mantid::DataObjects;

namespace Mantid {
namespace MDAlgorithms {
namespace {
// ---------- Template deduction of the event type
// --------------------------------
// See boost::type_traits documentation

/// Type trait to indicate that a general type is not a full MDEvent
template <typename MDE, size_t nd> struct IsFullEvent : boost::false_type {};

/// Specialization of type trait to indicate that a MDEvent is a full event
template <size_t nd> struct IsFullEvent<MDEvent<nd>, nd> : boost::true_type {};

/**
 * Specialization if isFullEvent for DataObjects
 * to return true
 */
template <typename MDE, size_t nd>
bool isFullMDEvent(const boost::true_type & /*unused*/) {
  return true;
}

/**
 * Specialization if isFullEvent for DataObjects
 * to return false
 */
template <typename MDE, size_t nd>
bool isFullMDEvent(const boost::false_type & /*unused*/) {
  return false;
}

/**
 * Returns true if the templated type is a full MDEvent
 */
template <typename MDE, size_t nd> bool isFullMDEvent() {
  return isFullMDEvent<MDE, nd>(IsFullEvent<MDE, nd>());
}

/**
 * Add the detectors from the given box as contributing detectors to the peak
 * @param peak :: The peak that relates to the box
 * @param box :: A reference to the box containing the peak
 */
template <typename MDE, size_t nd>
void addDetectors(DataObjects::Peak &peak, MDBoxBase<MDE, nd> &box,
                  const boost::true_type & /*unused*/) {
  if (box.getNumChildren() > 0) {
    std::cerr << "Box has children\n";
    addDetectors(peak, box, boost::true_type());
  }
  auto *mdBox = dynamic_cast<MDBox<MDE, nd> *>(&box);
  if (!mdBox) {
    throw std::invalid_argument("FindPeaksMD::addDetectors - Unexpected Box "
                                "type, cannot retrieve events");
  }
  const auto &events = mdBox->getConstEvents();
  auto itend = events.end();
  for (auto it = events.begin(); it != itend; ++it) {
    peak.addContributingDetID(it->getDetectorID());
  }
}

/// Add detectors based on lean events. Always throws as they do not know their
/// IDs
template <typename MDE, size_t nd>
void addDetectors(DataObjects::Peak & /*unused*/,
                  MDBoxBase<MDE, nd> & /*unused*/,
                  const boost::false_type & /*unused*/) {
  throw std::runtime_error("FindPeaksMD - Workspace contains lean events, "
                           "cannot include detector information");
}

/**
 * Add the detectors from the given box as contributing detectors to the peak
 * @param peak :: The peak that relates to the box
 * @param box :: A reference to the box containing the peak
 */
template <typename MDE, size_t nd>
void addDetectors(DataObjects::Peak &peak, MDBoxBase<MDE, nd> &box) {
  // Compile time deduction of the correct function call
  addDetectors(peak, box, IsFullEvent<MDE, nd>());
}
} // namespace

// Register the algorithm into the AlgorithmFactory
DECLARE_ALGORITHM(FindPeaksMD)

const std::string FindPeaksMD::volumeNormalization = "VolumeNormalization";
const std::string FindPeaksMD::numberOfEventsNormalization =
    "NumberOfEventsNormalization";

//----------------------------------------------------------------------------------------------
/** Constructor
 */
FindPeaksMD::FindPeaksMD()
    : peakWS(), peakRadiusSquared(), DensityThresholdFactor(0.0), m_maxPeaks(0),
      m_addDetectors(true), m_densityScaleFactor(1e-6), prog(nullptr), inst(),
      m_runNumber(-1), dimType(), m_goniometer() {}

//----------------------------------------------------------------------------------------------
/** Initialize the algorithm's properties.
 */
void FindPeaksMD::init() {
  declareProperty(std::make_unique<WorkspaceProperty<IMDWorkspace>>(
                      "InputWorkspace", "", Direction::Input),
                  "An input MDEventWorkspace or MDHistoWorkspace with at least "
                  "3 dimensions.");

  declareProperty(
      std::make_unique<PropertyWithValue<double>>("PeakDistanceThreshold", 0.1,
                                                  Direction::Input),
      "Threshold distance for rejecting peaks that are found to be too close "
      "from each other.\n"
      "This should be some multiple of the radius of a peak. Default: 0.1.");

  declareProperty(std::make_unique<PropertyWithValue<int64_t>>(
                      "MaxPeaks", 500, Direction::Input),
                  "Maximum number of peaks to find. Default: 500.");

  std::vector<std::string> strategy = {volumeNormalization,
                                       numberOfEventsNormalization};
  declareProperty(
      "PeakFindingStrategy", volumeNormalization,
      std::make_shared<StringListValidator>(strategy),
      "Strategy for finding peaks in an MD workspace."
      "1. VolumeNormalization: This is the default strategy. It will sort "
      "all boxes in the workspace by deacresing order of signal density "
      "(total weighted event sum divided by box volume).\n"
      "2.NumberOfEventsNormalization: This option is only valid for "
      "MDEventWorkspaces. It will use the total weighted event sum divided"
      "by the number of events. This can improve peak finding for "
      "histogram-based"
      "raw data which has been converted to an EventWorkspace. The threshold "
      "for"
      "peak finding can be controlled by the SingalThresholdFactor property "
      "which should"
      "be larger than 1. Note that this approach does not work for event-based "
      "raw data.\n");

  declareProperty(std::make_unique<PropertyWithValue<double>>(
                      "DensityThresholdFactor", 10.0, Direction::Input),
                  "The overall signal density of the workspace will be "
                  "multiplied by this factor \n"
                  "to get a threshold signal density below which boxes are NOT "
                  "considered to be peaks. See the help.\n"
                  "Default: 10.0");

  setPropertySettings("DensityThresholdFactor",
                      std::make_unique<EnabledWhenProperty>(
                          "PeakFindingStrategy",
                          Mantid::Kernel::ePropertyCriterion::IS_EQUAL_TO,
                          volumeNormalization));

  declareProperty(std::make_unique<PropertyWithValue<double>>(
                      "SignalThresholdFactor", 1.5, Direction::Input),
                  "The overal signal value (not density!) normalized by the "
                  "number of events is compared to the specified signal "
                  "threshold. Boxes which are below this threshold are NOT "
                  "considered to be peaks."
                  "This property is enabled when the PeakFindingStrategy has "
                  "been set to NumberOfEventsNormalization.\n"
                  "The value of boxes which contain peaks will be above 1. See "
                  "the below for more information.\n"
                  "Default: 1.50");

  setPropertySettings("SignalThresholdFactor",
                      std::make_unique<EnabledWhenProperty>(
                          "PeakFindingStrategy",
                          Mantid::Kernel::ePropertyCriterion::IS_EQUAL_TO,
                          numberOfEventsNormalization));

  declareProperty("CalculateGoniometerForCW", false,
                  "This will calculate the goniometer rotation (around y-axis "
                  "only) for a constant wavelength. This only works for Q "
                  "sample workspaces.");

  auto nonNegativeDbl = std::make_shared<BoundedValidator<double>>();
  nonNegativeDbl->setLower(0);
  declareProperty("Wavelength", DBL_MAX, nonNegativeDbl,
                  "Wavelength to use when calculating goniometer angle. If not"
                  "set will use the wavelength parameter on the instrument.");
  setPropertySettings("Wavelength",
                      std::make_unique<EnabledWhenProperty>(
                          "CalculateGoniometerForCW",
                          Mantid::Kernel::ePropertyCriterion::IS_NOT_DEFAULT));

  declareProperty("InnerGoniometer", false,
                  "Whether the goniometer to be calculated is the most inner "
                  "(phi) or most outer (omega)");
  setPropertySettings("InnerGoniometer",
                      std::make_unique<EnabledWhenProperty>(
                          "CalculateGoniometerForCW",
                          Mantid::Kernel::ePropertyCriterion::IS_NOT_DEFAULT));

  declareProperty("FlipX", false,
                  "Used when calculating goniometer angle if the q_lab x value "
                  "should be negative, hence the detector of the other side "
                  "(right) of the beam");
  setPropertySettings("FlipX",
                      std::make_unique<EnabledWhenProperty>(
                          "CalculateGoniometerForCW",
                          Mantid::Kernel::ePropertyCriterion::IS_NOT_DEFAULT));
  std::vector<std::string> peakTypes = {"Automatic", "Peak", "LeanElasticPeak"};
  declareProperty("OutputPeakType", "Automatic",
                  std::make_shared<StringListValidator>(peakTypes),
                  "Type of Peak in OutputWorkspace");
  declareProperty(std::make_unique<WorkspaceProperty<IPeaksWorkspace>>(
                      "OutputWorkspace", "", Direction::Output),
                  "An output PeaksWorkspace with the peaks' found positions.");

  declareProperty("AppendPeaks", false,
                  "If checked, then append the peaks in the output workspace "
                  "if it exists. \n"
                  "If unchecked, the output workspace is replaced (Default).");

  auto nonNegativeInt = std::make_shared<BoundedValidator<int>>();
  nonNegativeInt->setLower(0);
  declareProperty("EdgePixels", 0, nonNegativeInt,
                  "Remove peaks that are at pixels this close to edge. ");
}

//----------------------------------------------------------------------------------------------
/** Extract needed data from the workspace's experiment info */
void FindPeaksMD::readExperimentInfo(const ExperimentInfo_sptr &ei) {
  // Instrument associated with workspace
  inst = ei->getInstrument();
  // Find the run number
  m_runNumber = ei->getRunNumber();

  // Find the goniometer rotation matrix
  m_goniometer =
      Mantid::Kernel::Matrix<double>(3, 3, true); // Default IDENTITY matrix
  try {
    m_goniometer = ei->mutableRun().getGoniometerMatrix();
  } catch (std::exception &e) {
    g_log.warning() << "Error finding goniometer matrix. It will not be set in "
                       "the peaks found.\n";
    g_log.warning() << e.what() << '\n';
  }
}

void FindPeaksMD::checkWorkspaceDims(const IMDWorkspace_sptr &ws) {
  // Check that the workspace dimensions are in Q-sample-frame or Q-lab-frame.
  std::string dim0 = ws->getDimension(0)->getName();
  if (dim0 == "H") {
    dimType = HKL;
    throw std::runtime_error(
        "Cannot find peaks in a workspace that is already in HKL space.");
  } else if (dim0 == "Q_lab_x") {
    dimType = QLAB;
  } else if (dim0 == "Q_sample_x")
    dimType = QSAMPLE;
  else
    throw std::runtime_error(
        "Unexpected dimensions: need either Q_lab_x or Q_sample_x.");
}

void FindPeaksMD::determineOutputPeakType(const std::string peakType,
                                          const uint16_t nexp) {
  // This method will be expanded later to check a property on the
  // input workspace which can specify a default peak type for that
  // instrument.
  m_leanElasticPeak = false;
  if (peakType == "Automatic") {
    if (nexp == 0)
      m_leanElasticPeak = true;
  } else if (peakType == "LeanElasticPeak") {
    m_leanElasticPeak = true;
  } else { // Peak
    if (nexp == 0)
      throw std::runtime_error("Cannot create Peak output with 0 expInfo");
  }
}

//----------------------------------------------------------------------------------------------
/** Create and add a Peak to the output workspace
 *
 * @param Q :: Q_lab or Q_sample, depending on workspace
 * @param binCount :: bin count to give to the peak.
 * @param tracer :: Ray tracer to use for detector finding
 */
void FindPeaksMD::addPeak(const V3D &Q, const double binCount,
                          const Geometry::InstrumentRayTracer &tracer) {
  try {
    auto p = this->createPeak(Q, binCount, tracer);
    if (m_edge > 0) {
      if (edgePixel(inst, p->getBankName(), p->getCol(), p->getRow(), m_edge))
        return;
    }
    if (p->getDetectorID() != -1)
      peakWS->addPeak(*p);
  } catch (std::exception &e) {
    g_log.notice() << "Error creating peak at " << Q << " because of '"
                   << e.what() << "'. Peak will be skipped.\n";
  }
}

//----------------------------------------------------------------------------------------------
/** Create and add a LeanElasticPeak to the output workspace
 *
 * @param Q :: Q_lab or Q_sample, depending on workspace
 * @param binCount :: bin count to give to the peak.
 * @param useGoniometer :: if to include set goniometer from the input workspace
 * in the peak
 */
void FindPeaksMD::addLeanElasticPeak(const V3D &Q, const double binCount,
                                     const bool useGoniometer) {
  auto p = this->createLeanElasticPeak(Q, binCount, useGoniometer);
  peakWS->addPeak(*p);
}

/**
 * Creates a Peak object from Q & bin count
 * */
std::shared_ptr<DataObjects::Peak>
FindPeaksMD::createPeak(const Mantid::Kernel::V3D &Q, const double binCount,
                        const Geometry::InstrumentRayTracer &tracer) {
  std::shared_ptr<DataObjects::Peak> p;
  if (dimType == QLAB) {
    // Build using the Q-lab-frame constructor
    p = std::make_shared<Peak>(inst, Q);
    // Save gonio matrix for later
    p->setGoniometerMatrix(m_goniometer);
  } else if (dimType == QSAMPLE) {
    // Build using the Q-sample-frame constructor
    bool calcGoniometer = getProperty("CalculateGoniometerForCW");

    if (calcGoniometer) {
      // Calculate Q lab from Q sample and wavelength
      double wavelength = getProperty("Wavelength");
      if (wavelength == DBL_MAX) {
        if (inst->hasParameter("wavelength")) {
          wavelength = inst->getNumberParameter("wavelength").at(0);
        } else {
          throw std::runtime_error("Could not get wavelength, neither "
                                   "Wavelength algorithm property "
                                   "set nor instrument wavelength parameter");
        }
      }

      Geometry::Goniometer goniometer(m_goniometer);
      goniometer.calcFromQSampleAndWavelength(
          Q, wavelength, getProperty("FlipX"), getProperty("InnerGoniometer"));
      std::vector<double> angles = goniometer.getEulerAngles("YZY");
      g_log.information()
          << "Found goniometer rotation to be in YZY convention [" << angles[0]
          << ", " << angles[1] << ". " << angles[2]
          << "] degrees for Q sample = " << Q << "\n";
      p = std::make_shared<Peak>(inst, Q, goniometer.getR());

    } else {
      p = std::make_shared<Peak>(inst, Q, m_goniometer);
    }
  } else {
    throw std::invalid_argument(
        "Cannot Integrate peaks unless the dimension is QLAB or QSAMPLE");
  }

  try { // Look for a detector
    p->findDetector(tracer);
  } catch (...) { /* Ignore errors in ray-tracer */
  }

  p->setBinCount(binCount);
  // Save the run number found before.
  p->setRunNumber(m_runNumber);
  return p;
}

/**
 * Creates a Peak object from Q & bin count
 * */
std::shared_ptr<DataObjects::LeanElasticPeak>
FindPeaksMD::createLeanElasticPeak(const Mantid::Kernel::V3D &Q,
                                   const double binCount,
                                   const bool useGoniometer) {
  std::shared_ptr<DataObjects::LeanElasticPeak> p;
  if (dimType == QSAMPLE) {
    if (useGoniometer)
      p = std::make_shared<LeanElasticPeak>(Q, m_goniometer);
    else
      p = std::make_shared<LeanElasticPeak>(Q);
  } else {
    throw std::invalid_argument(
        "Cannot find peaks unless the dimension is QSAMPLE");
  }

  p->setBinCount(binCount);
  // Save the run number found before.
  p->setRunNumber(m_runNumber);
  return p;
}
//----------------------------------------------------------------------------------------------
/** Integrate the peaks of the workspace using parameters saved in the
 * algorithm class
 * @param ws ::  MDEventWorkspace to integrate
 */
template <typename MDE, size_t nd>
void FindPeaksMD::findPeaks(typename MDEventWorkspace<MDE, nd>::sptr ws) {
  if (nd < 3)
    throw std::invalid_argument("Workspace must have at least 3 dimensions.");

  if (isFullMDEvent<MDE, nd>()) {
    m_addDetectors = true;
  } else {
    m_addDetectors = false;
    g_log.warning("Workspace contains only lean events. Resultant "
                  "PeaksWorkspaces will not contain full detector "
                  "information.");
  }

  progress(0.01, "Refreshing Centroids");

  // TODO: This might be slow, progress report?
  // Make sure all centroids are fresh
  // ws->getBox()->refreshCentroid();

  // Calculate a threshold below which a box is too diffuse to be considered a
  // peak.
  signal_t threshold =
      m_useNumberOfEventsNormalization
          ? ws->getBox()->getSignalByNEvents() * m_signalThresholdFactor
          : ws->getBox()->getSignalNormalized() * DensityThresholdFactor;
  threshold *= m_densityScaleFactor;

  if (!std::isfinite(threshold)) {
    g_log.warning() << "Infinite or NaN overall density found. Your input data "
                       "may be invalid. Using a 0 threshold instead.\n";
    threshold = 0;
  }
  g_log.information() << "Threshold signal density: " << threshold << '\n';

  using boxPtr = API::IMDNode *;
  // We will fill this vector with pointers to all the boxes (up to a given
  // depth)
  typename std::vector<API::IMDNode *> boxes;

  // Get all the MDboxes
  progress(0.10, "Getting Boxes");
  ws->getBox()->getBoxes(boxes, 1000, true);

  // This pair is the <density, ptr to the box>
  using dens_box = std::pair<double, API::IMDNode *>;

  // Map that will sort the boxes by increasing density. The key = density;
  // value = box *.
  typename std::multimap<double, API::IMDNode *> sortedBoxes;

  // --------------- Sort and Filter by Density -----------------------------
  progress(0.20, "Sorting Boxes by Density");
  auto it1 = boxes.begin();
  auto it1_end = boxes.end();
  for (; it1 != it1_end; it1++) {
    auto box = *it1;
    double value = m_useNumberOfEventsNormalization
                       ? box->getSignalByNEvents()
                       : box->getSignalNormalized();
    value *= m_densityScaleFactor;
    // Skip any boxes with too small a signal value.
    if (value > threshold)
      sortedBoxes.insert(dens_box(value, box));
  }

  // --------------- Find Peak Boxes -----------------------------
  // List of chosen possible peak boxes.
  std::vector<API::IMDNode *> peakBoxes;

  prog = std::make_unique<Progress>(this, 0.30, 0.95, m_maxPeaks);

  // used for selecting method for calculating BinCount
  bool isMDEvent(ws->id().find("MDEventWorkspace") != std::string::npos);

  int64_t numBoxesFound = 0;
  // Now we go (backwards) through the map
  // e.g. from highest density down to lowest density.
  typename std::multimap<double, boxPtr>::reverse_iterator it2;
  auto it2_end = sortedBoxes.rend();
  for (it2 = sortedBoxes.rbegin(); it2 != it2_end; ++it2) {
    signal_t density = it2->first;
    boxPtr box = it2->second;
#ifndef MDBOX_TRACK_CENTROID
    coord_t boxCenter[nd];
    box->calculateCentroid(boxCenter);
#else
    const coord_t *boxCenter = box->getCentroid();
#endif

    // Compare to all boxes already picked.
    bool badBox = false;
    for (auto &peakBoxe : peakBoxes) {

#ifndef MDBOX_TRACK_CENTROID
      coord_t otherCenter[nd];
      (*it3)->calculateCentroid(otherCenter);
#else
      const coord_t *otherCenter = peakBoxe->getCentroid();
#endif

      // Distance between this box and a box we already put in.
      coord_t distSquared = 0.0;
      for (size_t d = 0; d < nd; d++) {
        coord_t dist = otherCenter[d] - boxCenter[d];
        distSquared += (dist * dist);
      }

      // Reject this box if it is too close to another previously found box.
      if (distSquared < peakRadiusSquared) {
        badBox = true;
        break;
      }
    }

    // The box was not rejected for another reason.
    if (!badBox) {
      if (numBoxesFound++ >= m_maxPeaks) {
        g_log.notice() << "Number of peaks found exceeded the limit of "
                       << m_maxPeaks << ". Stopping peak finding.\n";
        break;
      }

      peakBoxes.emplace_back(box);
      g_log.debug() << "Found box at ";
      for (size_t d = 0; d < nd; d++)
        g_log.debug() << (d > 0 ? "," : "") << boxCenter[d];
      g_log.debug() << "; Density = " << density << '\n';
      // Report progres for each box found.
      prog->report("Finding Peaks");
    }
  }

  prog->resetNumSteps(numBoxesFound, 0.95, 1.0);

  uint16_t nexp = ws->getNumExperimentInfo();

  if (nexp == 0) {
    // --- Convert the "boxes" to peaks ----
    for (auto box : peakBoxes) {
      // The center of the box = Q in the lab frame
#ifndef MDBOX_TRACK_CENTROID
      coord_t boxCenter[nd];
      box->calculateCentroid(boxCenter);
#else
      const coord_t *boxCenter = box->getCentroid();
#endif

      // Q of the centroid of the box
      V3D Q(boxCenter[0], boxCenter[1], boxCenter[2]);

      // The "bin count" used will be the box density.
      double binCount = box->getSignalNormalized() * m_densityScaleFactor;
      if (isMDEvent)
        binCount = static_cast<double>(box->getNPoints());

      // Create the peak
      addLeanElasticPeak(Q, binCount);

      // Report progres for each box found.
      prog->report("Adding Peaks");

    } // for each box found
  } else {
    for (uint16_t iexp = 0; iexp < ws->getNumExperimentInfo(); iexp++) {
      ExperimentInfo_sptr ei = ws->getExperimentInfo(iexp);
      this->readExperimentInfo(ei);

      Geometry::InstrumentRayTracer tracer(inst);
      // Copy the instrument, sample, run to the peaks workspace.
      peakWS->copyExperimentInfoFrom(ei.get());

      // --- Convert the "boxes" to peaks ----
      for (auto box : peakBoxes) {
        //  If no events from this experimental contribute to the box then skip
        if (nexp > 1) {
          auto *mdbox = dynamic_cast<MDBox<MDE, nd> *>(box);
          const std::vector<MDE> &events = mdbox->getEvents();
          if (std::none_of(events.cbegin(), events.cend(),
                           [&iexp, &nexp](MDE event) {
                             return event.getRunIndex() == iexp ||
                                    event.getRunIndex() >= nexp;
                           }))
            continue;
        }

        // If multiple goniometers than use the average one from the
        // events in the box, that matches this runIndex, this assumes
        // the events are added in same order as the goniometers
        if (ei->run().getNumGoniometers() > 1) {
          const std::vector<MDE> &events =
              dynamic_cast<MDBox<MDE, nd> *>(box)->getEvents();
          double sum = 0;
          double count = 0;
          for (const auto &event : events) {
            if (event.getRunIndex() == iexp) {
              sum += event.getGoniometerIndex();
              count++;
            }
          }
          m_goniometer =
              ei->mutableRun().getGoniometerMatrix(lrint(sum / count));
        }
        // The center of the box = Q in the lab frame

#ifndef MDBOX_TRACK_CENTROID
        coord_t boxCenter[nd];
        box->calculateCentroid(boxCenter);
#else
        const coord_t *boxCenter = box->getCentroid();
#endif

        // Q of the centroid of the box
        V3D Q(boxCenter[0], boxCenter[1], boxCenter[2]);

        // The "bin count" used will be the box density or the number of events
        // in the box
        double binCount = box->getSignalNormalized() * m_densityScaleFactor;
        if (isMDEvent)
          binCount = static_cast<double>(box->getNPoints());

        if (m_leanElasticPeak) {
          addLeanElasticPeak(Q, binCount, true);
        } else {
          try {
            auto p = this->createPeak(Q, binCount, tracer);
            if (m_addDetectors) {
              auto mdBox = dynamic_cast<MDBoxBase<MDE, nd> *>(box);
              if (!mdBox) {
                throw std::runtime_error("Failed to cast box to MDBoxBase");
              }
              addDetectors(*p, *mdBox);
            }
            if (p->getDetectorID() != -1) {
              if (m_edge > 0) {
                if (!edgePixel(inst, p->getBankName(), p->getCol(), p->getRow(),
                               m_edge))
                  peakWS->addPeak(*p);
                ;
              } else {
                peakWS->addPeak(*p);
              }
              g_log.information() << "Add new peak with Q-center = " << Q[0]
                                  << ", " << Q[1] << ", " << Q[2] << "\n";
            }
          } catch (std::exception &e) {
            g_log.notice() << "Error creating peak at " << Q << " because of '"
                           << e.what() << "'. Peak will be skipped.\n";
          }
        }

        // Report progress for each box found.
        prog->report("Adding Peaks");

      } // for each box found
    }
  }
  g_log.notice() << "Number of peaks found: " << peakWS->getNumberPeaks()
                 << '\n';
}

//----------------------------------------------------------------------------------------------
/** Find peaks in the given MDHistoWorkspace
 *
 * @param ws :: MDHistoWorkspace
 */
void FindPeaksMD::findPeaksHisto(
    const Mantid::DataObjects::MDHistoWorkspace_sptr &ws) {
  size_t nd = ws->getNumDims();
  if (nd < 3)
    throw std::invalid_argument("Workspace must have at least 3 dimensions.");

  g_log.warning("Workspace is an MDHistoWorkspace. Resultant PeaksWorkspaces "
                "will not contain full detector information.");

  // This pair is the <density, box index>
  using dens_box = std::pair<double, size_t>;

  // Map that will sort the boxes by increasing density. The key = density;
  // value = box index.
  std::multimap<double, size_t> sortedBoxes;

  size_t numBoxes = ws->getNPoints();

  // --------- Count the overall signal density -----------------------------
  progress(0.10, "Counting Total Signal");
  double totalSignal = 0;
  for (size_t i = 0; i < numBoxes; i++)
    totalSignal += ws->getSignalAt(i);
  // Calculate the threshold density
  double thresholdDensity =
      (totalSignal * ws->getInverseVolume() / double(numBoxes)) *
      DensityThresholdFactor * m_densityScaleFactor;
  if (!std::isfinite(thresholdDensity)) {
    g_log.warning() << "Infinite or NaN overall density found. Your input data "
                       "may be invalid. Using a 0 threshold instead.\n";
    thresholdDensity = 0;
  }
  g_log.information() << "Threshold signal density: " << thresholdDensity
                      << '\n';

  // -------------- Sort and Filter by Density -----------------------------
  progress(0.20, "Sorting Boxes by Density");
  for (size_t i = 0; i < numBoxes; i++) {
    double density = ws->getSignalNormalizedAt(i) * m_densityScaleFactor;
    // Skip any boxes with too small a signal density.
    if (density > thresholdDensity)
      sortedBoxes.insert(dens_box(density, i));
  }

  // --------------- Find Peak Boxes -----------------------------
  // List of chosen possible peak boxes.
  std::vector<size_t> peakBoxes;

  prog = std::make_unique<Progress>(this, 0.30, 0.95, m_maxPeaks);

  int64_t numBoxesFound = 0;
  // Now we go (backwards) through the map
  // e.g. from highest density down to lowest density.
  std::multimap<double, size_t>::reverse_iterator it2;
  auto it2_end = sortedBoxes.rend();
  for (it2 = sortedBoxes.rbegin(); it2 != it2_end; ++it2) {
    signal_t density = it2->first;
    size_t index = it2->second;
    // Get the center of the box
    VMD boxCenter = ws->getCenter(index);

    // Compare to all boxes already picked.
    bool badBox = false;
    for (auto &peakBoxe : peakBoxes) {
      VMD otherCenter = ws->getCenter(peakBoxe);

      // Distance between this box and a box we already put in.
      coord_t distSquared = 0.0;
      for (size_t d = 0; d < nd; d++) {
        coord_t dist = otherCenter[d] - boxCenter[d];
        distSquared += (dist * dist);
      }

      // Reject this box if it is too close to another previously found box.
      if (distSquared < peakRadiusSquared) {
        badBox = true;
        break;
      }
    }

    // The box was not rejected for another reason.
    if (!badBox) {
      if (numBoxesFound++ >= m_maxPeaks) {
        g_log.notice() << "Number of peaks found exceeded the limit of "
                       << m_maxPeaks << ". Stopping peak finding.\n";
        break;
      }

      peakBoxes.emplace_back(index);
      g_log.debug() << "Found box at index " << index;
      g_log.debug() << "; Density = " << density << '\n';
      // Report progres for each box found.
      prog->report("Finding Peaks");
    }
  }

  if (ws->getNumExperimentInfo() == 0) {
    // --- Convert the "boxes" to peaks ----
    for (auto index : peakBoxes) {
      // The center of the box = Q in the lab frame
      VMD boxCenter = ws->getCenter(index);

      // Q of the centroid of the box
      V3D Q(boxCenter[0], boxCenter[1], boxCenter[2]);

      // The "bin count" used will be the box density.
      double binCount = ws->getSignalNormalizedAt(index) * m_densityScaleFactor;

      // Create the peak
      addLeanElasticPeak(Q, binCount);

      // Report progres for each box found.
      prog->report("Adding Peaks");

    } // for each box found
  } else {
    for (uint16_t iexp = 0; iexp < ws->getNumExperimentInfo(); iexp++) {
      ExperimentInfo_sptr ei = ws->getExperimentInfo(iexp);
      this->readExperimentInfo(ei);
      Geometry::InstrumentRayTracer tracer(inst);

      // Copy the instrument, sample, run to the peaks workspace.
      peakWS->copyExperimentInfoFrom(ei.get());

      // --- Convert the "boxes" to peaks ----
      for (auto index : peakBoxes) {
        // The center of the box = Q in the lab frame
        VMD boxCenter = ws->getCenter(index);

        // Q of the centroid of the box
        V3D Q(boxCenter[0], boxCenter[1], boxCenter[2]);

        // The "bin count" used will be the box density.
        double binCount =
            ws->getSignalNormalizedAt(index) * m_densityScaleFactor;

        // Create the peak
        if (m_leanElasticPeak)
          addLeanElasticPeak(Q, binCount, true);
        else
          addPeak(Q, binCount, tracer);

        // Report progres for each box found.
        prog->report("Adding Peaks");

      } // for each box found
    }
  }
  g_log.notice() << "Number of peaks found: " << peakWS->getNumberPeaks()
                 << '\n';
} // namespace MDAlgorithms

//----------------------------------------------------------------------------------------------
/** Execute the algorithm.
 */
void FindPeaksMD::exec() {

  // The MDEventWorkspace as input
  IMDWorkspace_sptr inWS = getProperty("InputWorkspace");
  checkWorkspaceDims(inWS);
  MDHistoWorkspace_sptr inMDHW =
      std::dynamic_pointer_cast<MDHistoWorkspace>(inWS);
  IMDEventWorkspace_sptr inMDEW =
      std::dynamic_pointer_cast<IMDEventWorkspace>(inWS);

  bool AppendPeaks = getProperty("AppendPeaks");

  uint16_t nexp = 0;
  if (inMDHW)
    nexp = inMDHW->getNumExperimentInfo();
  else if (inMDEW)
    nexp = inMDEW->getNumExperimentInfo();

  std::string peakType = getProperty("OutputPeakType");

  determineOutputPeakType(peakType, nexp);

  // Output peaks workspace, create if needed
  peakWS = getProperty("OutputWorkspace");

  if (!peakWS || !AppendPeaks) {
    if (m_leanElasticPeak)
      peakWS = LeanElasticPeaksWorkspace_sptr(new LeanElasticPeaksWorkspace());
    else
      peakWS = PeaksWorkspace_sptr(new PeaksWorkspace());
  }

  // Other parameters
  double PeakDistanceThreshold = getProperty("PeakDistanceThreshold");
  peakRadiusSquared =
      static_cast<coord_t>(PeakDistanceThreshold * PeakDistanceThreshold);

  DensityThresholdFactor = getProperty("DensityThresholdFactor");
  m_signalThresholdFactor = getProperty("SignalThresholdFactor");

  std::string strategy = getProperty("PeakFindingStrategy");
  m_useNumberOfEventsNormalization = strategy == numberOfEventsNormalization;

  m_maxPeaks = getProperty("MaxPeaks");
  m_edge = this->getProperty("EdgePixels");

  // Execute the proper algo based on the type of workspace
  if (inMDHW) {
    this->findPeaksHisto(inMDHW);
  } else if (inMDEW) {
    CALL_MDEVENT_FUNCTION3(this->findPeaks, inMDEW);
  } else {
    throw std::runtime_error("This algorithm can only find peaks on a "
                             "MDHistoWorkspace or a MDEventWorkspace; it does "
                             "not work on a regular MatrixWorkspace.");
  }

  // Do a sort by bank name and then descending bin count (intensity)
  std::vector<std::pair<std::string, bool>> criteria;
  criteria.emplace_back("RunNumber", true);
  criteria.emplace_back("BankName", true);
  criteria.emplace_back("bincount", false);
  peakWS->sort(criteria);

  for (auto i = 0; i != peakWS->getNumberPeaks(); ++i) {
    Mantid::Geometry::IPeak &p = peakWS->getPeak(i);
    p.setPeakNumber(i + 1);
  }
  // Save the output
  setProperty("OutputWorkspace", peakWS);
}

//----------------------------------------------------------------------------------------------
/** Validate the inputs.
 */
std::map<std::string, std::string> FindPeaksMD::validateInputs() {
  std::map<std::string, std::string> result;
  // Check for number of event normalzation
  std::string strategy = getProperty("PeakFindingStrategy");

  const bool useNumberOfEventsNormalization =
      strategy == numberOfEventsNormalization;
  IMDWorkspace_sptr inWS = getProperty("InputWorkspace");
  IMDEventWorkspace_sptr inMDEW =
      std::dynamic_pointer_cast<IMDEventWorkspace>(inWS);
  MDHistoWorkspace_sptr inMDHW =
      std::dynamic_pointer_cast<MDHistoWorkspace>(inWS);

  if (useNumberOfEventsNormalization && !inMDEW) {
    result["PeakFindingStrategy"] = "The NumberOfEventsNormalization selection "
                                    "can only be used with an MDEventWorkspace "
                                    "as the input.";
  }

  uint16_t nexp = 0;
  if (inMDHW)
    nexp = inMDHW->getNumExperimentInfo();
  else if (inMDEW)
    nexp = inMDEW->getNumExperimentInfo();

  std::string peakType = getProperty("OutputPeakType");

  if (peakType == "Peak" && nexp == 0)
    result["OutputPeakType"] =
        "The InputWorkspace doesn't contain any experiment information so the "
        "OutputPeakType cannot be Peak.";

  return result;
}

} // namespace MDAlgorithms
} // namespace Mantid

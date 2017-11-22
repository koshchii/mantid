#include "MantidDataHandling/LoadILLReflectometry.h"

#include "MantidAPI/Axis.h"
#include "MantidAPI/FileProperty.h"
#include "MantidAPI/FunctionFactory.h"
#include "MantidAPI/IPeakFunction.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidAPI/RegisterFileLoader.h"
#include "MantidAPI/SpectrumInfo.h"
#include "MantidAPI/WorkspaceFactory.h"
#include "MantidDataObjects/TableWorkspace.h"
#include "MantidHistogramData/LinearGenerator.h"
#include "MantidKernel/BoundedValidator.h"
#include "MantidKernel/EnabledWhenProperty.h"
#include "MantidKernel/ListValidator.h"
#include "MantidKernel/OptionalBool.h"
#include "MantidKernel/Quat.h"
#include "MantidKernel/UnitFactory.h"
#include "MantidDataObjects/Workspace2D.h"
#include "MantidDataObjects/WorkspaceCreation.h"

namespace {
/// Component coordinates for Figaro, in meter.
namespace Figaro {
// TODO: Consider moving these to the IPF.
constexpr double detectorRestY{0.509};
constexpr double DH1X{1.135}; //Motor DH1 horzontal position
constexpr double DH2X{2.077}; //Motor DH2 horzontal position
}

/// A struct for information needed for detector angle calibration.
struct PeakInfo {
  double detectorAngle;
  double detectorDistance;
  double peakCentre;
};

/** Convert degrees to radians.
 *  @param x an angle in degrees
 *  @return the angle in radians
 */
constexpr double inRad(const double x) { return x * M_PI / 180; }

/** Convert radians to degrees.
 *  @param x an angle in radians
 *  @return the angle in degrees
 */
constexpr double inDeg(const double x) { return x * 180 / M_PI; }

/** Convert millimeters to meters.
 *  @param x a distance in millimeters
 *  @return the distance in meters
 */
constexpr double inMeter(const double x) { return x * 1e-3; }

/** Create a table with data needed for detector angle calibration.
 * @param info data to be written to the table
 * @return a TableWorkspace containing the beam position info
 */
Mantid::API::ITableWorkspace_sptr
createPeakPositionTable(const PeakInfo &info) {
  auto table = Mantid::API::WorkspaceFactory::Instance().createTable();
  table->addColumn("double", "DetectorAngle");
  table->addColumn("double", "DetectorDistance");
  table->addColumn("double", "PeakCentre");
  table->appendRow();
  auto col = table->getColumn("DetectorAngle");
  col->cell<double>(0) = info.detectorAngle;
  col = table->getColumn("DetectorDistance");
  col->cell<double>(0) = info.detectorDistance;
  col = table->getColumn("PeakCentre");
  col->cell<double>(0) = info.peakCentre;
  return table;
}

/** Strip monitors from the beginning and end of a workspace.
 *  @param ws a workspace to work on
 *  @return begin and end ws indices for non-monitor histograms
 */
std::pair<size_t, size_t>
fitIntegrationWSIndexRange(const Mantid::API::MatrixWorkspace &ws) {
  const size_t nHisto = ws.getNumberHistograms();
  size_t begin = 0;
  const auto &spectrumInfo = ws.spectrumInfo();
  for (size_t i = 0; i < nHisto; ++i) {
    if (!spectrumInfo.isMonitor(i)) {
      break;
    }
    ++begin;
  }
  size_t end = nHisto - 1;
  for (ptrdiff_t i = static_cast<ptrdiff_t>(nHisto) - 1; i != 0; --i) {
    if (!spectrumInfo.isMonitor(i)) {
      break;
    }
    --end;
  }
  return std::pair<size_t, size_t>{begin, end};
}

/** Construct a DirectBeamMeasurement object from a beam position table.
 *  @param table a beam position TableWorkspace
 *  @return a DirectBeamMeasurement object corresonding to the table parameter.
 */
PeakInfo parseBeamPositionTable(const Mantid::API::ITableWorkspace &table) {
  if (table.rowCount() != 1) {
    throw std::runtime_error("BeamPosition table should have a single row.");
  }
  PeakInfo p;
  auto col = table.getColumn("DetectorAngle");
  p.detectorAngle = col->cell<double>(0);
  col = table.getColumn("DetectorDistance");
  p.detectorDistance = col->cell<double>(0);
  col = table.getColumn("PeakCentre");
  p.peakCentre = col->cell<double>(0);
  return p;
}

/** Calculate the offset angle between detector center and peak.
 *  @param peakCentre peak centre in pixels.
 *  @param detectorCentre detector centre in pixels.
 *  @param pixelWidth pixel width in meters.
 *  @param detectorDistance detector-sample distance in meters.
 *  @return the offset angle.
 */
double offsetAngle(const double peakCentre, const double detectorCentre,
                   const double pixelWidth, const double detectorDistance) {
  const double offsetWidth = (detectorCentre - peakCentre) * pixelWidth;
  return inDeg(std::atan2(offsetWidth, detectorDistance));
}

/** Fill the X values of the first histogram of ws with values 0, 1, 2,...
 *  @param ws a workspace to modify
 */
void rebinIntegralWorkspace(Mantid::API::MatrixWorkspace &ws) {
  auto &xs = ws.mutableX(0);
  std::iota(xs.begin(), xs.end(), 0.0);
}

/// Enumerations to define the rotation plane of the detector.
enum class RotationPlane { horizontal, vertical };

/** Calculate the detector position from given parameters.
 *  @param plane rotation plane of the detector
 *  @param distance sample to detector centre distance in meters
 *  @param angle an angle between the Z axis and the detector in degrees
 *  @return a vector pointing to the new detector centre
 */
Mantid::Kernel::V3D detectorPosition(const RotationPlane plane,
                                     const double distance,
                                     const double angle) {
  const double a = inRad(angle);
  double x, y, z;
  switch (plane) {
  case RotationPlane::horizontal:
    x = distance * std::sin(a);
    y = 0;
    z = distance * std::cos(a);
    break;
  case RotationPlane::vertical:
    x = 0;
    y = distance * std::sin(a);
    z = distance * std::cos(a);
    break;
  }
  return Mantid::Kernel::V3D(x, y, z);
}

/** Calculates the detector rotation such that it faces the origin.
 *  @param plane rotation plane of the detectorPosition
 *  @param angle an angle between the Z axis and the detector in degrees
 *  @return the calculated rotation transformation
 */
Mantid::Kernel::Quat detectorFaceRotation(const RotationPlane plane,
                                          const double angle) {
  const Mantid::Kernel::V3D axis = [plane]() {
    double x, y;
    switch (plane) {
    case RotationPlane::horizontal:
      x = 0;
      y = 1;
      break;
    case RotationPlane::vertical:
      x = -1;
      y = 0;
      break;
    }
    return Mantid::Kernel::V3D(x, y, 0);
  }();
  return Mantid::Kernel::Quat(angle, axis);
}
} // anonymous namespace

namespace Mantid {
namespace DataHandling {

using namespace Kernel;
using namespace API;
using namespace NeXus;

// Register the algorithm into the AlgorithmFactory
DECLARE_NEXUS_FILELOADER_ALGORITHM(LoadILLReflectometry)

/**
 * Return the confidence with this algorithm can load the file
 * @param descriptor A descriptor for the file
 * @returns An integer specifying the confidence level. 0 indicates it will not
 * be used
 */
int LoadILLReflectometry::confidence(
    Kernel::NexusDescriptor &descriptor) const {

  // fields existent only at the ILL
  if ((descriptor.pathExists("/entry0/wavelength") || // ILL D17
       descriptor.pathExists("/entry0/theta"))        // ILL Figaro
      &&
      descriptor.pathExists("/entry0/experiment_identifier") &&
      descriptor.pathExists("/entry0/mode") &&
      (descriptor.pathExists("/entry0/instrument/VirtualChopper") || // ILL D17
       descriptor.pathExists("/entry0/instrument/Theta")) // ILL Figaro
      )
    return 80;
  else
    return 0;
}

/// Initialize the algorithm's properties.
void LoadILLReflectometry::init() {
  declareProperty(Kernel::make_unique<FileProperty>("Filename", std::string(),
                                                    FileProperty::Load, ".nxs",
                                                    Direction::Input),
                  "Name of the Nexus file to load");

  declareProperty(Kernel::make_unique<WorkspaceProperty<>>(
                      "OutputWorkspace", std::string(), Direction::Output),
                  "Name of the output workspace");

  declareProperty(Kernel::make_unique<WorkspaceProperty<ITableWorkspace>>(
                      "OutputBeamPosition", std::string(), Direction::Output,
                      PropertyMode::Optional),
                  "Name of the fitted beam position output workspace");

  declareProperty(Kernel::make_unique<WorkspaceProperty<ITableWorkspace>>(
                      "BeamPosition", std::string(), Direction::Input,
                      PropertyMode::Optional),
                  "A workspace defining the beam position; used to calculate "
                  "the Bragg angle");

  declareProperty("BraggAngle", EMPTY_DBL(),
                  "User defined Bragg angle in degrees");
  const std::vector<std::string> availableUnits{"Wavelength", "TimeOfFlight"};
  declareProperty("XUnit", "Wavelength",
                  boost::make_shared<StringListValidator>(availableUnits),
                  "X unit of the OutputWorkspace");
}

/**
 * Validate inputs
 * @returns a string map containing the error messages
 */
std::map<std::string, std::string> LoadILLReflectometry::validateInputs() {
  std::map<std::string, std::string> result;
  return result;
}

/// Execute the algorithm.
void LoadILLReflectometry::exec() {
  // open the root node
  NeXus::NXRoot root(getPropertyValue("Filename"));
  NXEntry firstEntry{root.openFirstEntry()};
  // load Monitor details: n. monitors x monitor contents
  std::vector<std::vector<int>> monitorsData{loadMonitors(firstEntry)};
  // set instrument specific names of Nexus file entries
  initNames(firstEntry);
  // load Data details (number of tubes, channels, etc)
  loadDataDetails(firstEntry);
  // initialise workspace
  initWorkspace(monitorsData);
  // load the instrument from the IDF if it exists
  loadInstrument();
  // get properties
  loadNexusEntriesIntoProperties();
  // load data into the workspace
  loadData(firstEntry, monitorsData, getXValues());
  root.close();
  firstEntry.close();
  // Move components as if the sample was at the origin (it usually is).
  placeSource();
  placeDetector();
  // When other components are in-place
  placeSample();
  convertTofToWavelength();
  // Set the output workspace property
  setProperty("OutputWorkspace", m_localWorkspace);
} // exec

/// Run the Child Algorithm LoadInstrument.
void LoadILLReflectometry::loadInstrument() {
  // execute the Child Algorithm. Catch and log any error, but don't stop.
  g_log.debug("Loading instrument definition...");
  try {
    IAlgorithm_sptr loadInst = createChildAlgorithm("LoadInstrument");
    loadInst->setPropertyValue("InstrumentName", m_instrumentName);
    loadInst->setProperty("RewriteSpectraMap",
                          Mantid::Kernel::OptionalBool(true));
    loadInst->setProperty<MatrixWorkspace_sptr>("Workspace", m_localWorkspace);
    loadInst->executeAsChildAlg();
  } catch (std::runtime_error &e) {
    g_log.information()
        << "Unable to succesfully run LoadInstrument child algorithm: "
        << e.what() << '\n';
  }
}

/**
  * Init names of member variables based on instrument specific NeXus file
  * entries
  *
  * @param entry :: the NeXus file entry
  */
void LoadILLReflectometry::initNames(NeXus::NXEntry &entry) {
  std::string instrumentNamePath = m_loader.findInstrumentNexusPath(entry);
  m_instrumentName = entry.getString(instrumentNamePath.append("/name"));
  if (m_instrumentName.empty())
    throw std::runtime_error(
        "Cannot set the instrument name from the Nexus file!");
  // In NeXus files names are: D17 and figaro. The instrument
  // definition is independent and names start with a capital letter. This
  // loader follows its convention.
  boost::to_lower(m_instrumentName);
  m_instrumentName[0] = char((std::toupper(m_instrumentName[0])));
  g_log.debug() << "Instrument name: " << m_instrumentName << '\n';
  if (m_instrumentName == "D17") {
    m_detectorDistanceName = "det";
    m_detectorAngleName = "dan.value";
    m_sampleAngleName = "san.value";
    m_offsetFrom = "VirtualChopper";
    m_offsetName = "open_offset";
    m_pixelCentre = 127.5;
    m_chopper1Name = "Chopper1";
    m_chopper2Name = "Chopper2";
  } else if (m_instrumentName == "Figaro") {
    // For Figaro, the DTR field contains the sample-to-detector distance
    // when the detector is at the horizontal position (angle = 0).
    m_detectorDistanceName = "DTR";
    m_detectorAngleName = "VirtualAxis.DAN_actual_angle";
    m_sampleAngleName = "CollAngle.actual_coll_angle";
    m_offsetFrom = "CollAngle";
    m_offsetName = "openOffset";
    m_pixelCentre = 127.5;
    // Figaro: find out which of the four choppers are used
    NXFloat firstChopper =
        entry.openNXFloat("instrument/ChopperSetting/firstChopper");
    firstChopper.load();
    NXFloat secondChopper =
        entry.openNXFloat("instrument/ChopperSetting/secondChopper");
    secondChopper.load();
    m_chopper1Name = "CH" + std::to_string(int(firstChopper[0]));
    m_chopper2Name = "CH" + std::to_string(int(secondChopper[0]));
  }
  // get acquisition mode
  NXInt acqMode = entry.openNXInt("acquisition_mode");
  acqMode.load();
  m_acqMode = acqMode[0];
  m_acqMode ? g_log.debug("TOF mode") : g_log.debug("Monochromatic Mode");
}

/// Call child algorithm ConvertUnits for conversion from TOF to wavelength
void LoadILLReflectometry::convertTofToWavelength() {
  if (m_acqMode && (getPropertyValue("XUnit") == "Wavelength")) {
    auto convertToWavelength =
        createChildAlgorithm("ConvertUnits", -1, -1, true);
    convertToWavelength->initialize();
    convertToWavelength->setProperty<MatrixWorkspace_sptr>("InputWorkspace",
                                                           m_localWorkspace);
    convertToWavelength->setProperty<MatrixWorkspace_sptr>("OutputWorkspace",
                                                           m_localWorkspace);
    convertToWavelength->setPropertyValue("Target", "Wavelength");
    convertToWavelength->executeAsChildAlg();
  }
}

/**
 * Creates the workspace and initialises member variables with
 * the corresponding values
 *
 * @param monitorsData :: Monitors data already loaded
 */
void LoadILLReflectometry::initWorkspace(
    const std::vector<std::vector<int>> &monitorsData) {

  g_log.debug() << "Number of monitors: " << monitorsData.size() << '\n';
  for (size_t i = 0; i < monitorsData.size(); ++i) {
    if (monitorsData[i].size() != m_numberOfChannels)
      g_log.debug() << "Data size of monitor ID " << i << " is "
                    << monitorsData[i].size() << '\n';
  }
  // create the workspace
  try {
    m_localWorkspace = WorkspaceFactory::Instance().create(
        "Workspace2D", m_numberOfHistograms + monitorsData.size(),
        m_numberOfChannels + 1, m_numberOfChannels);
  } catch (std::out_of_range &) {
    throw std::runtime_error(
        "Workspace2D cannot be created, check number of histograms (" +
        std::to_string(m_numberOfHistograms) + "), monitors (" +
        std::to_string(monitorsData.size()) + "), and channels (" +
        std::to_string(m_numberOfChannels) + '\n');
  }
  if (m_acqMode)
    m_localWorkspace->getAxis(0)->unit() =
        UnitFactory::Instance().create("TOF");
  m_localWorkspace->setYUnitLabel("Counts");
  m_localWorkspace->mutableRun().addProperty("Facility", std::string("ILL"));
}

/**
 * Load Data details (number of tubes, channels, etc)
 *
 * @param entry First entry of nexus file
 */
void LoadILLReflectometry::loadDataDetails(NeXus::NXEntry &entry) {
  // PSD data D17 256 x 1 x 1000
  // PSD data Figaro 1 x 256 x 1000

  if (m_acqMode) {
    NXFloat timeOfFlight = entry.openNXFloat("instrument/PSD/time_of_flight");
    timeOfFlight.load();
    m_channelWidth = static_cast<double>(timeOfFlight[0]);
    m_numberOfChannels = size_t(timeOfFlight[1]);
    m_tofDelay = timeOfFlight[2];
    if (m_instrumentName == "Figaro") {
      NXFloat eDelay = entry.openNXFloat("instrument/Theta/edelay_delay");
      eDelay.load();
      m_tofDelay += static_cast<double>(eDelay[0]);
    }
  } else { // monochromatic mode
    m_numberOfChannels = 1;
  }

  NXInt nChannels = entry.openNXInt("instrument/PSD/detsize");
  nChannels.load();
  m_numberOfHistograms = nChannels[0];

  std::string widthName;
  if (m_instrumentName == "D17")
    widthName = "mppx";
  else if (m_instrumentName == "Figaro")
    widthName = "mppy";

  NXFloat pixelWidth = entry.openNXFloat("instrument/PSD/" + widthName);
  pixelWidth.load();
  m_pixelWidth = inMeter(static_cast<double>(pixelWidth[0]));

  g_log.debug()
      << "Please note that ILL reflectometry instruments have "
         "several tubes, after integration one "
         "tube remains in the Nexus file.\n Number of tubes (banks): 1\n";
  g_log.debug() << "Number of pixels per tube (number of detectors and number "
                   "of histograms): " << m_numberOfHistograms << '\n';
  g_log.debug() << "Number of time channels: " << m_numberOfChannels << '\n';
  g_log.debug() << "Channel width: " << m_channelWidth << " 1e-6 sec\n";
  g_log.debug() << "TOF delay: " << m_tofDelay << '\n';
  g_log.debug() << "Pixel width: " << m_pixelWidth << '\n';
}

double LoadILLReflectometry::doubleFromRun(const std::string &entryName) const {
  return m_localWorkspace->run().getPropertyValueAsType<double>(entryName);
}

/**
 * Load single monitor
 *
 * @param entry :: The Nexus entry
 * @param monitor_data :: A std::string containing the Nexus path to the monitor
 *data
 * @return monitor :: A std::vector containing monitor values
 */
std::vector<int>
LoadILLReflectometry::loadSingleMonitor(NeXus::NXEntry &entry,
                                        const std::string &monitor_data) {
  NXData dataGroup = entry.openNXData(monitor_data);
  NXInt data = dataGroup.openIntData();
  // load counts
  data.load();
  return std::vector<int>(data(), data() + data.size());
}

/**
 * Load monitor data
 *
 * @param entry :: The Nexus entry
 * @return :: A std::vector of vectors of monitors containing monitor values
 */
std::vector<std::vector<int>>
LoadILLReflectometry::loadMonitors(NeXus::NXEntry &entry) {
  g_log.debug("Read monitor data...");
  // vector of monitors with one entry
  const std::vector<std::vector<int>> monitors{
      loadSingleMonitor(entry, "monitor1/data"),
      loadSingleMonitor(entry, "monitor2/data")};
  return monitors;
}

/**
 * Determine x values (unit time-of-flight)
 *
 * @return :: vector holding the x values
 */
std::vector<double> LoadILLReflectometry::getXValues() {
  std::vector<double> xVals;             // no initialisation
  xVals.reserve(m_numberOfChannels + 1); // reserve memory
  try {
    if (m_acqMode) {
      std::string chopper{"Chopper"};
      double chop1Speed{0.0}, chop1Phase{0.0}, chop2Speed{0.0}, chop2Phase{0.0};
      if (m_instrumentName == "D17") {
        chop1Speed = doubleFromRun("VirtualChopper.chopper1_speed_average");
        chop1Phase = doubleFromRun("VirtualChopper.chopper1_phase_average");
        chop2Speed = doubleFromRun("VirtualChopper.chopper2_speed_average");
        chop2Phase = doubleFromRun("VirtualChopper.chopper2_phase_average");
      } else if (m_instrumentName == "Figaro") {
        chop1Phase = doubleFromRun(m_chopper1Name + ".phase");
        // Chopper 1 phase on Figaro is set to an arbitrary value (999.9)
        if (chop1Phase > 360.0)
          chop1Phase = 0.0;
      }
      const double POFF = doubleFromRun(m_offsetFrom + ".poff");
      const double openOffset =
          doubleFromRun(m_offsetFrom + "." + m_offsetName);
      if (m_instrumentName == "D17" && chop1Speed != 0.0 && chop2Speed != 0.0 &&
          chop2Phase != 0.0) {
        // virtual chopper entries are valid
        chopper = "Virtual chopper";
      } else {
        // use chopper values
        chop1Speed = doubleFromRun(m_chopper1Name + ".rotation_speed");
        chop2Speed = doubleFromRun(m_chopper2Name + ".rotation_speed");
        chop2Phase = doubleFromRun(m_chopper2Name + ".phase");
      }
      // logging
      g_log.debug() << "Poff: " << POFF << '\n';
      g_log.debug() << "Open offset: " << openOffset << '\n';
      g_log.debug() << "Chopper 1 phase: " << chop1Phase << '\n';
      g_log.debug() << chopper << " 1 speed: " << chop1Speed << '\n';
      g_log.debug() << chopper << " 2 phase: " << chop2Phase << '\n';
      g_log.debug() << chopper << " 2 speed: " << chop2Speed << '\n';

      if (chop1Speed <= 0.0) {
        g_log.error() << "First chopper velocity " << chop1Speed
                      << ". Check you NeXus file.\n";
      }
      const double chopWindow = 45.0;
      const double t_TOF2 = m_tofDelay -
                            1.e+6 * 60.0 * (POFF - chopWindow + chop2Phase -
                                            chop1Phase + openOffset) /
                                (2.0 * 360 * chop1Speed);
      g_log.debug() << "t_TOF2: " << t_TOF2 << '\n';
      // compute tof values
      xVals.emplace_back(t_TOF2 - 0.5 * m_channelWidth);
      for (int channelIndex = 0;
           channelIndex < static_cast<int>(m_numberOfChannels);
           ++channelIndex) {
        const double t_TOF1 = (channelIndex + 0.5) * m_channelWidth;
        xVals.push_back(t_TOF1 + t_TOF2);
      }
    } else {
      g_log.debug("Time channel index for axis description \n");
      for (size_t t = 0; t <= m_numberOfChannels; ++t)
        xVals.push_back(double(t));
    }
  } catch (std::runtime_error &e) {
    g_log.information() << "Unable to access NeXus file entry: " << e.what()
                        << '\n';
  }
  return xVals;
}

/**
 * Load data from nexus file
 *
 * @param entry :: The Nexus file entry
 * @param monitorsData :: Monitors data already loaded
 * @param xVals :: X values
 */
void LoadILLReflectometry::loadData(
    NeXus::NXEntry &entry, const std::vector<std::vector<int>> &monitorsData,
    const std::vector<double> &xVals) {
  g_log.debug("Loading data...");
  NXData dataGroup = entry.openNXData("data");
  NXInt data = dataGroup.openIntData();
  // load the counts from the file into memory
  data.load();
  const size_t nb_monitors = monitorsData.size();
  Progress progress(this, 0, 1, m_numberOfHistograms + nb_monitors);

  // write monitors
  if (!xVals.empty()) {
    HistogramData::BinEdges binEdges(xVals);
    // write data
    for (size_t j = 0; j < m_numberOfHistograms; ++j) {
      const int *data_p = &data(0, static_cast<int>(j), 0);
      const HistogramData::Counts counts(data_p, data_p + m_numberOfChannels);
      m_localWorkspace->setHistogram(j, binEdges, std::move(counts));
      progress.report();
      for (size_t im = 0; im < nb_monitors; ++im) {
        const int *monitor_p = monitorsData[im].data();
        const HistogramData::Counts counts(monitor_p,
                                           monitor_p + m_numberOfChannels);
        m_localWorkspace->setHistogram(im + m_numberOfHistograms, binEdges,
                                       std::move(counts));
        progress.report();
      }
    }
  } else
    g_log.debug("Vector of x values is empty");
} // LoadILLIndirect::loadData

/**
 * Use the LoadHelper utility to load most of the nexus entries into workspace
 * sample log properties
 */
void LoadILLReflectometry::loadNexusEntriesIntoProperties() {
  g_log.debug("Building properties...");
  // Open NeXus file
  const std::string filename{getPropertyValue("Filename")};
  NXhandle nxfileID;
  NXstatus stat = NXopen(filename.c_str(), NXACC_READ, &nxfileID);
  if (stat == NX_ERROR)
    throw Kernel::Exception::FileError("Unable to open File:", filename);
  m_loader.addNexusFieldsToWsRun(nxfileID, m_localWorkspace->mutableRun());
  stat = NXclose(&nxfileID);
}

/**
  * Gaussian fit to determine peak position.
  *
  * @return :: detector position of the peak: Gaussian fit and position
  * of the maximum (serves as start value for the optimization)
  */
double LoadILLReflectometry::fitReflectometryPeak() {
  size_t startIndex;
  size_t endIndex;
  std::tie(startIndex, endIndex) =
      fitIntegrationWSIndexRange(*m_localWorkspace);
  IAlgorithm_sptr integration = createChildAlgorithm("Integration");
  integration->initialize();
  integration->setProperty("InputWorkspace", m_localWorkspace);
  integration->setProperty("OutputWorkspace", "__unused_for_child");
  integration->setProperty("StartWorkspaceIndex", static_cast<int>(startIndex));
  integration->setProperty("EndWorkspaceIndex", static_cast<int>(endIndex));
  integration->execute();
  MatrixWorkspace_sptr integralWS = integration->getProperty("OutputWorkspace");
  IAlgorithm_sptr transpose = createChildAlgorithm("Transpose");
  transpose->initialize();
  transpose->setProperty("InputWorkspace", integralWS);
  transpose->setProperty("OutputWorkspace", "__unused_for_child");
  transpose->execute();
  integralWS = transpose->getProperty("OutputWorkspace");
  rebinIntegralWorkspace(*integralWS);
  // determine initial height: maximum value
  const auto maxValueIt =
      std::max_element(integralWS->y(0).cbegin(), integralWS->y(0).cend());
  const double height = *maxValueIt;
  // determine initial centre: index of the maximum value
  const size_t maxIndex = std::distance(integralWS->y(0).cbegin(), maxValueIt);
  const double centreByMax = static_cast<double>(maxIndex);
  g_log.debug() << "Peak maximum position: " << centreByMax << '\n';
  // determine sigma
  auto lessThanHalfMax = [height](const double x) { return x < 0.5 * height; };
  using IterType = HistogramData::HistogramY::const_iterator;
  std::reverse_iterator<IterType> revMaxValueIt{maxValueIt};
  auto revMinFwhmIt =
      std::find_if(revMaxValueIt, integralWS->y(0).crend(), lessThanHalfMax);
  auto maxFwhmIt =
      std::find_if(maxValueIt, integralWS->y(0).cend(), lessThanHalfMax);
  std::reverse_iterator<IterType> revMaxFwhmIt{maxFwhmIt};
  const double fwhm =
      static_cast<double>(std::distance(revMaxFwhmIt, revMinFwhmIt) + 1);
  g_log.debug() << "Initial fwhm (fixed window at half maximum): " << fwhm
                << '\n';
  // generate Gaussian
  auto func = API::FunctionFactory::Instance().createFunction("Gaussian");
  auto initialGaussian = boost::dynamic_pointer_cast<API::IPeakFunction>(func);
  initialGaussian->setHeight(height);
  initialGaussian->setCentre(centreByMax);
  initialGaussian->setFwhm(fwhm);
  // call Fit child algorithm
  API::IAlgorithm_sptr fitGaussian = createChildAlgorithm("Fit");
  fitGaussian->initialize();
  fitGaussian->setProperty(
      "Function", boost::dynamic_pointer_cast<API::IFunction>(initialGaussian));
  fitGaussian->setProperty("InputWorkspace", integralWS);
  bool success = fitGaussian->execute();
  if (!success)
    g_log.warning("Fit not successful, using position of max value.\n");
  else
    g_log.debug() << "Sigma: " << initialGaussian->fwhm() << '\n';
  const double centre = success ? initialGaussian->centre() : centreByMax;
  g_log.debug() << "Estimated peak position: " << centre << '\n';
  return centre;
}

/** Compute the detector and the Bragg angles.
 *  @return a pair where first is detector angle and second the Bragg angle.
 */
std::pair<double, double> LoadILLReflectometry::detectorAndBraggAngles() {
  ITableWorkspace_const_sptr posTable = getProperty("BeamPosition");
  const double peakCentre = fitReflectometryPeak();
  g_log.debug() << "Using detector angle (degrees): " << m_detectorAngle << '\n';
  if (!isDefault("OutputBeamPosition")) {
    PeakInfo p;
    p.detectorAngle = m_detectorAngle;
    p.detectorDistance = m_detectorDistance;
    p.peakCentre = peakCentre;
    setProperty("OutputBeamPosition", createPeakPositionTable(p));
  }
  const double userAngle = getProperty("BraggAngle");
  const double offset = offsetAngle(peakCentre, m_pixelCentre, m_pixelWidth,
                                    m_detectorDistance);
  m_log.debug() << "Beam offset angle: " << offset << '\n';
  if (userAngle != EMPTY_DBL()) {
    if (posTable) {
      g_log.notice() << "Ignoring BeamPosition, using BraggAngle instead.";
    }
    const double userDetectorAngle = 2 * userAngle - offset;
    return std::make_pair(userDetectorAngle, userAngle);
  }
  if (!posTable) {
    const double bragg = (m_detectorAngle + offset) / 2;
    return std::make_pair(m_detectorAngle, bragg);
  }
  const auto dbPeak = parseBeamPositionTable(*posTable);
  const double dbOffset = offsetAngle(dbPeak.peakCentre, m_pixelCentre,
                                      m_pixelWidth, dbPeak.detectorDistance);
  m_log.debug() << "Direct beam offset angle: " << dbOffset << '\n';
  const double detectorAngle =
      m_detectorAngle - dbPeak.detectorAngle - dbOffset;
  m_log.debug() << "Direct beam calibrated detector angle: " << detectorAngle
                << '\n';
  const double bragg = (detectorAngle + offset) / 2;
  return std::make_pair(detectorAngle, bragg);
}

/// Update detector position according to data file
void LoadILLReflectometry::placeDetector() {
  g_log.debug("Move the detector bank \n");
  m_detectorDistance = sampleDetectorDistance();
  m_detectorAngle = detectorAngle();
  g_log.debug() << "Sample-detector distance: " << m_detectorDistance
                << "m.\n";
  double detectorAngle;
  double braggAngle;
  std::tie(detectorAngle, braggAngle) = detectorAndBraggAngles();
  // incident angle for using the algorithm ConvertToReflectometryQ
  // TODO Doesn't seem to work with ConvertToReflectometryQ. Maybe they
  //      expect a time series?
  // TODO They are moving to 2theta in ISIS reflectometry algorithms.
  m_localWorkspace->mutableRun().addProperty("stheta", inRad(braggAngle));
  const std::string componentName = "detector";
  const RotationPlane rotPlane = [this]() {
    if (m_instrumentName == "D17")
      return RotationPlane::horizontal;
    else if (m_instrumentName == "Figaro")
      return RotationPlane::vertical;
    else
      return RotationPlane::horizontal;
  }();
  const auto newpos =
      detectorPosition(rotPlane, m_detectorDistance, detectorAngle);
  m_loader.moveComponent(m_localWorkspace, componentName, newpos);
  // apply a local rotation to stay perpendicular to the beam
  const auto rotation = detectorFaceRotation(rotPlane, detectorAngle);
  m_loader.rotateComponent(m_localWorkspace, componentName, rotation);
}

/// Update sample position.
void LoadILLReflectometry::placeSample() {
  if (m_instrumentName != "Figaro") {
    // Accept the sample position defined in the IDF.
    return;
  }
  m_sampleZOffset = inMeter(doubleFromRun("Theta.sampleHorizontalOffset"));
  const V3D newPos{0.0, 0.0, m_sampleZOffset};
  m_loader.moveComponent(m_localWorkspace, "sample_position", newPos);
}

/// Update source position.
void LoadILLReflectometry::placeSource() {
  const double dist = sourceSampleDistance();
  g_log.debug() << "Source-sample distance " << dist << "m.\n";
  const std::string source = "chopper1";
  const V3D newPos{0.0, 0.0, -dist};
  m_loader.moveComponent(m_localWorkspace, source, newPos);
}

double LoadILLReflectometry::detectorAngle() const
{
  if (m_instrumentName != "Figaro") {
    return doubleFromRun(m_detectorAngleName);
  }
  const double DH1Y = inMeter(doubleFromRun("DH1.value"));
  const double DH2Y = inMeter(doubleFromRun("DH2.value"));
  return inDeg(std::atan2(DH2Y - DH1Y, Figaro::DH2X - Figaro::DH1X));
}

/** Return the sample to detector distance for the current instrument.
 *  @return the sample to detector distance in meters
 */
double LoadILLReflectometry::sampleDetectorDistance() const {
  if (m_instrumentName != "Figaro") {
    return inMeter(doubleFromRun(m_detectorDistanceName + ".value"));
  }
  const double detectorRestX = inMeter(doubleFromRun(m_detectorDistanceName + ".value"));
  // Motor DH1 vertical coordinate.
  const double DH1Y = inMeter(doubleFromRun("DH1.value"));
  const double detAngle = detectorAngle();
  const double detectorX = std::cos(inRad(detAngle)) * (detectorRestX - Figaro::DH1X) + Figaro::DH1X;
  const double detectorY = std::sin(inRad(detAngle)) * (detectorRestX - Figaro::DH1X) + DH1Y - Figaro::detectorRestY;
  const double pixelOffset = Figaro::detectorRestY - 0.5 * m_pixelWidth;
  const double beamX = detectorX - pixelOffset * std::sin(inRad(detAngle));
  const double beamY = detectorY + pixelOffset * std::cos(inRad(detAngle));
  return std::hypot(beamX, beamY);
}

/** Return the source to sample distance for the current instrument.
 *  @return the source to sample distance in meters
 */
double LoadILLReflectometry::sourceSampleDistance() const {
  if (m_instrumentName == "D17") {
    const double pairCentre = doubleFromRun("VirtualChopper.dist_chop_samp");
    // Chopper pair separation is in cm in sample logs.
    const double pairSeparation = doubleFromRun("Distance.ChopperGap") / 100;
    return pairCentre - 0.5 * pairSeparation;
  } else if (m_instrumentName == "Figaro") {
    return inMeter(doubleFromRun("ChopperSetting.chopperpair_sample_distance"));
  }
  std::ostringstream out;
  out << "sourceSampleDistance: unknown instrument " << m_instrumentName;
  throw std::runtime_error(out.str());
}

} // namespace DataHandling
} // namespace Mantid

// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +

#include "MantidDataHandling/LoadILLSingleCrystal.h"
#include "MantidAPI/FileProperty.h"
#include "MantidAPI/Run.h"
#include "MantidAPI/Sample.h"
#include "MantidAPI/WorkspaceFactory.h"
#include "MantidDataHandling/H5Util.h"
#include "MantidDataObjects/MDHistoWorkspace.h"
#include "MantidGeometry/Crystal/CrystalStructure.h"
#include "MantidGeometry/Crystal/IsotropicAtomBraggScatterer.h"
#include "MantidGeometry/Crystal/PointGroupFactory.h"
#include "MantidGeometry/Crystal/SpaceGroupFactory.h"
#include "MantidGeometry/Instrument/DetectorInfo.h"
#include "MantidGeometry/MDGeometry/MDFrame.h"
#include "MantidGeometry/MDGeometry/MDFrameFactory.h"
#include "MantidGeometry/MDGeometry/UnknownFrame.h"

/*
 @class LoadILLSingleCrystal
 @author T. Weber, ILL
 @date 9/2021
 */
namespace Mantid {
namespace DataHandling {
using Mantid::API::WorkspaceProperty;
using Mantid::Kernel::Direction;

// Register the algorithm into the AlgorithmFactory
DECLARE_ALGORITHM(LoadILLSingleCrystal)

//----------------------------------------------------------------------------------------------

/// Algorithms name for identification. @see Algorithm::name
const std::string LoadILLSingleCrystal::name() const { return "LoadILLSingleCrystal"; }

/// Algorithm's version for identification. @see Algorithm::version
int LoadILLSingleCrystal::version() const { return 1; }

/// Algorithm's category for identification. @see Algorithm::category
const std::string LoadILLSingleCrystal::category() const { return "DataHandling\\Nexus;ILL\\Diffraction"; }

/// Algorithm's summary for use in the GUI and help. @see Algorithm::summary
const std::string LoadILLSingleCrystal::summary() const { return "Loads ILL single-crystal diffraction nexus files."; }

//----------------------------------------------------------------------------------------------
/** Helper functions
 */
template <class T> static std::vector<T> get_values(NeXus::File &file, const std::string &key) {
  std::vector<T> values;

  try {
    file.readData(key, values);
  } catch (const std::exception &e) {
  }

  return values;
}

template <> std::vector<std::string> get_values(NeXus::File &file, const std::string &key) {
  std::vector<std::string> values;

  file.openData(key);

  // TODO

  file.closeData();

  return values;
}

template <class T> static boost::optional<T> get_value(NeXus::File &file, const std::string &key) {
  std::vector<T> values = get_values<T>(file, key);

  if (!values.size())
    return boost::none;

  return values[0];
}

template <> boost::optional<std::string> get_value(NeXus::File &file, const std::string &key) {
  std::string str;

  try {
    file.readData(key, str);
  } catch (const std::exception &e) {
  }

  return str;
}

static std::vector<int64_t> get_dims(NeXus::File &file, const std::string &key) {
  file.openData(key);
  NeXus::Info infos = file.getInfo();
  file.closeData();

  return infos.dims;
}

//----------------------------------------------------------------------------------------------
/** Initialize the algorithm's properties.
 */
void LoadILLSingleCrystal::init() {
  // input file
  declareProperty(std::make_unique<API::FileProperty>("Filename", "", API::FileProperty::Load, ".nxs"),
                  "File path of the data file to load");

  // verbose log flag
  // declareProperty(std::make_unique<Kernel::PropertyWithValue<bool>>("Verbose Logging", false,
  //                Direction::Input), "Show verbose logs.");

  // output workspace
  declareProperty(std::make_unique<WorkspaceProperty<API::IMDHistoWorkspace>>("OutputWorkspace", "", Direction::Output),
                  "An output workspace.");
}

//----------------------------------------------------------------------------------------------
/** Execute the algorithm.
 */
void LoadILLSingleCrystal::exec() {
  // m_log.setEnabled(getProperty("Verbose Logging"));

  // get input file name
  m_filename = getPropertyValue("Filename");

  // open the corresponding nexus file
  m_file = std::make_unique<NeXus::File>(m_filename, NXACC_READ);
  if (!m_file)
    throw Kernel::Exception::FileError("Can not open file ", m_filename);

  // --------------------------------------------------------------------------
  // root group
  m_file->openGroup("/entry0", "NXentry");
  auto entries = m_file->getEntries();

  m_log.information() << "\nRoot entry group:" << std::endl;
  for (const auto &entry : entries)
    m_log.information() << entry.first << " = " << entry.second << std::endl;

  // title
  std::string title = *get_value<std::string>(*m_file, "title");

  // run number
  int numor = *get_value<int>(*m_file, "run_number");

  // wavelength
  float wavelength = *get_value<float>(*m_file, "wavelength");
  // --------------------------------------------------------------------------

  // --------------------------------------------------------------------------
  // instrument group
  m_file->openGroup("instrument", "NXinstrument");
  entries = m_file->getEntries();

  m_log.information() << "\nInstrument group:" << std::endl;
  for (const auto &entry : entries)
    m_log.information() << entry.first << " = " << entry.second << std::endl;

  // instrument name
  std::string instr_name = *get_value<std::string>(*m_file, "name");

  // initial crystal angles
  m_file->openGroup("chi", "NXpositioner");
  float chi = *get_value<float>(*m_file, "value");
  m_file->closeGroup();
  m_file->openGroup("omega", "NXpositioner");
  float omega = *get_value<float>(*m_file, "value");
  m_file->closeGroup();
  m_file->openGroup("phi", "NXpositioner");
  float phi = *get_value<float>(*m_file, "value");
  m_file->closeGroup();

  m_log.information() << "Initial sample angles: "
                      << "chi=" << chi << ", omega=" << omega << ", phi=" << phi << std::endl;

  // crystal
  m_file->openGroup("SingleCrystalSettings", "NXcrystal");
  auto _cell_a = get_value<float>(*m_file, "unit_cell_a");
  auto _cell_b = get_value<float>(*m_file, "unit_cell_b");
  auto _cell_c = get_value<float>(*m_file, "unit_cell_c");
  auto _cell_alpha = get_value<float>(*m_file, "unit_cell_alpha");
  auto _cell_beta = get_value<float>(*m_file, "unit_cell_beta");
  auto _cell_gamma = get_value<float>(*m_file, "unit_cell_gamma");

  // if a is not given, set it to 0
  float cell_a = _cell_a ? *_cell_a : 0;
  // if b is not given, set it to a
  float cell_b = _cell_b ? *_cell_b : cell_a;
  // if c is not given, set it to a
  float cell_c = _cell_c ? *_cell_c : cell_a;
  // if alpha is not given, set it to 90 deg
  float cell_alpha = _cell_alpha ? *_cell_alpha : 90;
  // if beta is not given, set it to alpha
  float cell_beta = _cell_beta ? *_cell_beta : cell_alpha;
  // if gamma is not given, set it to alpha
  float cell_gamma = _cell_gamma ? *_cell_gamma : cell_alpha;

  m_log.information() << "Unit cell: a=" << cell_a << ", b=" << cell_b << ", c=" << cell_c << ", alpha=" << cell_alpha
                      << ", beta=" << cell_beta << ", gamma=" << cell_gamma << std::endl;

  auto UB = get_values<float>(*m_file, "orientation_matrix");
  m_log.information() << "UB matrix:\n"
                      << "\t" << std::setw(15) << UB[0] << std::setw(15) << UB[1] << std::setw(15) << UB[2] << "\n"
                      << "\t" << std::setw(15) << UB[3] << std::setw(15) << UB[4] << std::setw(15) << UB[5] << "\n"
                      << "\t" << std::setw(15) << UB[6] << std::setw(15) << UB[7] << std::setw(15) << UB[8]
                      << std::endl;

  m_file->closeGroup(); // SingleCrystalSettings

  m_file->closeGroup(); // instrument
  // --------------------------------------------------------------------------

  // --------------------------------------------------------------------------
  // sample group
  m_file->openGroup("sample", "NXsample");
  entries = m_file->getEntries();

  m_log.information() << "\nSample group:" << std::endl;
  for (const auto &entry : entries)
    m_log.information() << entry.first << " = " << entry.second << std::endl;

  m_file->closeGroup(); // sample
  // --------------------------------------------------------------------------

  // --------------------------------------------------------------------------
  // user group
  m_file->openGroup("user", "NXuser");
  entries = m_file->getEntries();

  m_log.information() << "\nUser group:" << std::endl;
  for (const auto &entry : entries)
    m_log.information() << entry.first << " = " << entry.second << std::endl;

  m_file->closeGroup(); // user
  // --------------------------------------------------------------------------

  // --------------------------------------------------------------------------
  // monitor group
  m_file->openGroup("monitor", "NXmonitor");
  entries = m_file->getEntries();

  m_log.information() << "\nMonitor group:" << std::endl;
  for (const auto &entry : entries)
    m_log.information() << entry.first << " = " << entry.second << std::endl;

  // total monitor counts
  float monitor_sum = *get_value<float>(*m_file, "monsum");

  m_file->closeGroup(); // monitor
  // --------------------------------------------------------------------------

  // --------------------------------------------------------------------------
  // data_scan group
  m_file->openGroup("data_scan", "NXdata");
  entries = m_file->getEntries();

  m_log.information() << "\nData group:" << std::endl;
  for (const auto &entry : entries)
    m_log.information() << entry.first << " = " << entry.second << std::endl;

  // steps
  int total_steps = *get_value<int>(*m_file, "total_steps");
  int actual_step = *get_value<int>(*m_file, "actual_step");

  m_log.information() << "Total steps=" << total_steps << ", actual step=" << actual_step << std::endl;

  // detector data
  m_file->openGroup("detector_data", "ILL_detectors_data_scan");

  auto detector_data_dims = get_dims(*m_file, "data");
  auto detector_data = get_values<uint32_t>(*m_file, "data");

  if (detector_data_dims.size() < 3)
    throw std::runtime_error("Detector data dimension < 3.");

  m_log.information() << "Detector data dimensions: ";
  for (std::int16_t i = 0; i < detector_data_dims.size(); ++i) {
    m_log.information() << detector_data_dims[i];
    if (i < detector_data_dims.size() - 1)
      m_log.information() << " x ";
  }
  m_log.information() << std::endl;

  m_file->closeGroup(); // detector_data

  // scanned variables
  m_file->openGroup("scanned_variables", "ILL_data_scan_vars");
  auto scanned_vars_entries = m_file->getEntries();

  m_log.information() << "\nData group -> Scanned variables:" << std::endl;
  for (const auto &entry : scanned_vars_entries)
    m_log.information() << entry.first << " = " << entry.second << std::endl;

  auto scanned_variables_values = get_values<double>(*m_file, "data");

  m_file->openGroup("variables_names", "ILL_data_scan_vars");
  auto variables_names_entries = m_file->getEntries();

  m_log.information() << "\nData group -> Scanned variables -> Variables names:" << std::endl;
  for (const auto &entry : variables_names_entries)
    m_log.information() << entry.first << " = " << entry.second << std::endl;

  auto scanned_variables_names = get_values<std::string>(*m_file, "name");
  auto scanned_variables = get_values<std::uint8_t>(*m_file, "scanned");

  // get index of scanned variable
  std::uint8_t scanned_idx;
  for (std::size_t idx = 0; idx < scanned_variables.size(); ++idx) {
    if (scanned_variables[idx]) {
      scanned_idx = idx;
      break;
    }
  }
  m_log.information() << "Scanned variable index: " << (int)scanned_idx << std::endl;

  m_file->closeGroup(); // variable names

  m_file->closeGroup(); // scanned_variables

  m_file->closeGroup(); // data_scan
  // --------------------------------------------------------------------------

  m_file->closeGroup();

  // create workspace
  Geometry::GeneralFrame frame_x("x", "");
  Geometry::GeneralFrame frame_y("y", "");
  Geometry::GeneralFrame frame_scanned("scanned", "");

  std::vector<Mantid::Geometry::MDHistoDimension_sptr> dimensions{{
      std::make_shared<Geometry::MDHistoDimension>("y", "y", frame_y, 0, detector_data_dims[2] - 1,
                                                   detector_data_dims[2]),
      std::make_shared<Geometry::MDHistoDimension>("x", "x", frame_x, 0, detector_data_dims[1] - 1,
                                                   detector_data_dims[1]),
      std::make_shared<Geometry::MDHistoDimension>("Scanned Variable", "scanned", frame_scanned, 0,
                                                   detector_data_dims[0] - 1, detector_data_dims[0]),
  }};

  m_workspace = std::make_shared<DataObjects::MDHistoWorkspace>(dimensions, API::NoNormalization);

  // set the metadata
  auto info = std::make_shared<API::ExperimentInfo>();
  info->mutableRun().addProperty("Filename", m_filename);
  info->mutableRun().addProperty("Instrument", instr_name);
  info->mutableRun().addProperty("Wavelength", wavelength);
  info->mutableRun().addProperty("Numor", numor);
  info->mutableRun().addProperty("Monitor", monitor_sum);
  info->mutableRun().addProperty("Sample_a", cell_a);
  info->mutableRun().addProperty("Sample_b", cell_b);
  info->mutableRun().addProperty("Sample_c", cell_c);
  info->mutableRun().addProperty("Sample_alpha", cell_alpha);
  info->mutableRun().addProperty("Sample_beta", cell_beta);
  info->mutableRun().addProperty("Sample_gamma", cell_gamma);
  for (std::size_t i = 0; i < UB.size(); ++i)
    info->mutableRun().addProperty("Sample_UB_" + std::to_string(i), UB[i]);

  Geometry::CrystalStructure crys(
      // unit cell
      Geometry::UnitCell(cell_a, cell_b, cell_c, cell_alpha, cell_beta, cell_gamma, Geometry::angDegrees),
      // dummy space group
      Geometry::SpaceGroupFactory::Instance().createSpaceGroup("P 1"),
      // dummy scatterer
      Geometry::CompositeBraggScatterer::create(Geometry::IsotropicAtomBraggScattererParser("")()));
  info->mutableSample().setCrystalStructure(crys);

  m_workspace->addExperimentInfo(info);
  m_workspace->setTitle(title);

  // size of one detector image in pixels
  std::int64_t frame_size = detector_data_dims[1] * detector_data_dims[2];

  // copy the detector data to the workspace
  std::int64_t last_scan_step = -1;
  for (std::size_t idx = 0; idx < detector_data.size(); ++idx) {
    // completed all pixels of a frame?
    std::int64_t cur_scan_step = idx / frame_size;
    if (last_scan_step != cur_scan_step) {
      double scanned_variable_value = scanned_variables_values[detector_data_dims[0] * scanned_idx + cur_scan_step];

      // std::cout << scanned_variable_value << std::endl;
      // TODO: associate scanned_variable_value with this frame (i.e. cur_scan_step)

      last_scan_step = cur_scan_step;
    }

    m_workspace->setSignalAt(idx, detector_data[idx]);
  }

  // set ouput workspace
  setProperty("OutputWorkspace", std::dynamic_pointer_cast<API::IMDHistoWorkspace>(m_workspace));
}

} // namespace DataHandling
} // namespace Mantid

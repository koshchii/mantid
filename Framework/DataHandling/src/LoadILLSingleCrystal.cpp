// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +

#include "MantidDataHandling/LoadILLSingleCrystal.h"
#include "MantidAPI/FileProperty.h"
#include "MantidAPI/WorkspaceFactory.h"
#include "MantidDataHandling/H5Util.h"
#include "MantidDataObjects/MDHistoWorkspace.h"
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
template <class T> static boost::optional<T> get_value(NeXus::File &file, const std::string &key) {
  std::vector<T> values;
  file.readData(key, values);

  if (!values.size())
    return boost::none;

  return values[0];
}

template <> boost::optional<std::string> get_value(NeXus::File &file, const std::string &key) {
  std::string str;
  file.readData(key, str);
  return str;
}

template <class T> static std::vector<T> get_values(NeXus::File &file, const std::string &key) {
  std::vector<T> values;
  file.readData(key, values);

  return values;
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

  // output workspace
  declareProperty(std::make_unique<WorkspaceProperty<API::IMDHistoWorkspace>>("OutputWorkspace", "", Direction::Output),
                  "An output workspace.");
}

//----------------------------------------------------------------------------------------------
/** Execute the algorithm.
 */
void LoadILLSingleCrystal::exec() {
  // get input file name
  m_filename = getPropertyValue("Filename");

  // open the corresponding nexus file
  m_file = std::make_unique<NeXus::File>(m_filename, NXACC_READ);

  // --------------------------------------------------------------------------
  // root group
  m_file->openGroup("/entry0", "NXentry");
  auto entries = m_file->getEntries();

  std::cout << "\nRoot entry group:" << std::endl;
  for (const auto &entry : entries)
    std::cout << entry.first << " = " << entry.second << std::endl;

  // run number
  int numor = *get_value<int>(*m_file, "run_number");

  // wavelength
  float wavelength = *get_value<float>(*m_file, "wavelength");
  // --------------------------------------------------------------------------

  // --------------------------------------------------------------------------
  // instrument group
  m_file->openGroup("instrument", "NXinstrument");
  entries = m_file->getEntries();

  std::cout << "\nInstrument group:" << std::endl;
  for (const auto &entry : entries)
    std::cout << entry.first << " = " << entry.second << std::endl;

  // instrument name
  std::string instr_name = *get_value<std::string>(*m_file, "name");

  m_file->closeGroup(); // instrument
  // --------------------------------------------------------------------------

  // --------------------------------------------------------------------------
  // sample group
  m_file->openGroup("sample", "NXsample");
  entries = m_file->getEntries();

  std::cout << "\nSample group:" << std::endl;
  for (const auto &entry : entries)
    std::cout << entry.first << " = " << entry.second << std::endl;

  m_file->closeGroup(); // sample
  // --------------------------------------------------------------------------

  // --------------------------------------------------------------------------
  // user group
  m_file->openGroup("user", "NXuser");
  entries = m_file->getEntries();

  std::cout << "\nUser group:" << std::endl;
  for (const auto &entry : entries)
    std::cout << entry.first << " = " << entry.second << std::endl;

  m_file->closeGroup(); // user
  // --------------------------------------------------------------------------

  // --------------------------------------------------------------------------
  // monitor group
  m_file->openGroup("monitor", "NXmonitor");
  entries = m_file->getEntries();

  std::cout << "\nMonitor group:" << std::endl;
  for (const auto &entry : entries)
    std::cout << entry.first << " = " << entry.second << std::endl;

  // total monitor counts
  float monitor_sum = *get_value<float>(*m_file, "monsum");

  m_file->closeGroup(); // monitor
  // --------------------------------------------------------------------------

  // --------------------------------------------------------------------------
  // data_scan group
  m_file->openGroup("data_scan", "NXdata");
  entries = m_file->getEntries();

  std::cout << "\nData group:" << std::endl;
  for (const auto &entry : entries)
    std::cout << entry.first << " = " << entry.second << std::endl;

  // steps
  int total_steps = *get_value<int>(*m_file, "total_steps");
  int actual_step = *get_value<int>(*m_file, "actual_step");

  // detector data
  m_file->openGroup("detector_data", "ILL_detectors_data_scan");

  auto detector_data_dims = get_dims(*m_file, "data");
  auto detector_data = get_values<uint32_t>(*m_file, "data");

  m_file->closeGroup(); // detector_data

  // scanned variables
  m_file->openGroup("scanned_variables", "ILL_data_scan_vars");
  auto scanned_vars_entries = m_file->getEntries();

  std::cout << "\nData group -> Scanned variables:" << std::endl;
  for (const auto &entry : scanned_vars_entries)
    std::cout << entry.first << " = " << entry.second << std::endl;

  auto scanned_variables_values = get_values<double>(*m_file, "data");

  m_file->openGroup("variables_names", "ILL_data_scan_vars");
  auto variables_names_entries = m_file->getEntries();
  std::cout << "\nData group -> Scanned variables -> Variables names:" << std::endl;
  for (const auto &entry : variables_names_entries)
    std::cout << entry.first << " = " << entry.second << std::endl;

  // auto scanned_variables_names = get_values<std::string>(*m_file, "names");
  auto scanned_variables = get_values<unsigned char>(*m_file, "scanned");

  m_file->closeGroup(); // variable names

  m_file->closeGroup(); // scanned_variables

  m_file->closeGroup(); // data_scan
  // --------------------------------------------------------------------------

  m_file->closeGroup();

  // create workspace
  Geometry::GeneralFrame frame_x("x", "");
  Geometry::GeneralFrame frame_y("y", "");
  Geometry::GeneralFrame frame_z("z", "");

  std::vector<Mantid::Geometry::MDHistoDimension_sptr> dimensions{{
      std::make_shared<Geometry::MDHistoDimension>("x", "x", frame_x, 0, detector_data_dims[2] - 1,
                                                   detector_data_dims[2]),
      std::make_shared<Geometry::MDHistoDimension>("y", "y", frame_y, 0, detector_data_dims[1] - 1,
                                                   detector_data_dims[1]),
      std::make_shared<Geometry::MDHistoDimension>("z", "z", frame_z, 0, detector_data_dims[0] - 1,
                                                   detector_data_dims[0]),
  }};
  m_workspace = std::make_shared<DataObjects::MDHistoWorkspace>(dimensions, API::NoNormalization);

  // copy the detector data to the workspace
  for (std::size_t idx = 0; idx < detector_data.size(); ++idx)
    m_workspace->setSignalAt(idx, detector_data[idx]);

  // set ouput workspace
  setProperty("OutputWorkspace", std::dynamic_pointer_cast<API::IMDHistoWorkspace>(m_workspace));
}

} // namespace DataHandling
} // namespace Mantid

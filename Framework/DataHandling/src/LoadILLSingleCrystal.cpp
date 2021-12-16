// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +

#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>
#include <hdf5.h>

#include "MantidDataHandling/H5Util.h"
#include "MantidDataHandling/LoadILLSingleCrystal.h"

#include "MantidAPI/FileProperty.h"
#include "MantidAPI/RegisterFileLoader.h"
#include "MantidAPI/Run.h"
#include "MantidAPI/Sample.h"
#include "MantidAPI/WorkspaceFactory.h"

#include "MantidKernel/OptionalBool.h"

#include "MantidDataObjects/MDEvent.h"
#include "MantidDataObjects/MDEventWorkspace.h"
#include "MantidDataObjects/MDHistoWorkspace.h"

#include "MantidGeometry/Crystal/CrystalStructure.h"
#include "MantidGeometry/Crystal/IsotropicAtomBraggScatterer.h"
#include "MantidGeometry/Crystal/OrientedLattice.h"
#include "MantidGeometry/Crystal/PointGroupFactory.h"
#include "MantidGeometry/Crystal/SpaceGroupFactory.h"
#include "MantidGeometry/Instrument.h"
#include "MantidGeometry/Instrument/Component.h"
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

// sample angle names
#define SAMPLE_CHI_NAME "chi"
#define SAMPLE_OMEGA_NAME "omega"
#define SAMPLE_PHI_NAME "phi"

// Register the algorithm into the AlgorithmFactory
// DECLARE_ALGORITHM(LoadILLSingleCrystal);
DECLARE_NEXUS_FILELOADER_ALGORITHM(LoadILLSingleCrystal)

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
/** Helper functions:  file handling
 */

// get numerical values from a nexus file
template <class T> static std::vector<T> get_values(NeXus::File &file, const std::string &key) {
  std::vector<T> values;

  try {
    file.readData(key, values);
  } catch (const std::exception &e) {
  }

  return values;
}

// get strings from a nexus file (apparently not supported by the nexus library!)
/*template <> std::vector<std::string> get_values(NeXus::File &file, const std::string &key) {
  std::vector<std::string> values;
  file.openData(key);

  std::string valuestr = file.getStrData();
  //NeXus::Info info = file.getInfo();

  // TODO: add this function (and remove the raw hdf5 code) as soon as the nexus library supports it

  file.closeData();
  return values;
}*/

// get a numerical value from a nexus file
template <class T> static boost::optional<T> get_value(NeXus::File &file, const std::string &key) {
  std::vector<T> values = get_values<T>(file, key);

  if (!values.size())
    return boost::none;

  return values[0];
}

// get a string from a nexus file
template <> boost::optional<std::string> get_value(NeXus::File &file, const std::string &key) {
  std::string str;

  try {
    file.readData(key, str);
  } catch (const std::exception &e) {
  }

  return str;
}

// get the dimensions of a dataset from a nexus file
static std::vector<std::int64_t> get_dims(NeXus::File &file, const std::string &key) {
  file.openData(key);
  NeXus::Info infos = file.getInfo();
  file.closeData();

  return infos.dims;
}

// get the dimensions of a dataset from a hdf5 dataset
static std::vector<hsize_t> get_dims(hid_t data) {
  std::vector<hsize_t> dims;

  hid_t space = H5Dget_space(data);
  if (space < 0)
    return dims;

  int rank = H5Sget_simple_extent_ndims(space);
  dims.resize(rank);
  H5Sget_simple_extent_dims(space, dims.data(), 0);

  H5Sclose(space);
  return dims;
}

// get strings from a raw hdf5 file
static std::vector<std::string> get_values(hid_t file, const std::string &key) {
  std::vector<std::string> values;

  hid_t data = H5Dopen(file, key.c_str(), H5P_DEFAULT);
  if (data < 0)
    return values;

  auto dims = get_dims(data);
  if (dims.size() == 1) {
    hid_t type_id = H5Tcopy(H5T_C_S1);
    H5Tset_size(type_id, H5T_VARIABLE);

    hsize_t num_strings = dims[0];
    if (num_strings) {
      const char **str = new const char *[num_strings];
      if (H5Dread(data, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, str) >= 0) {
        for (hsize_t i = 0; i < num_strings; ++i)
          values.push_back(str[i]);
      }
      delete[] str;
    }
  }

  H5Dclose(data);
  return values;
}
//----------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------
/** Helper function for rotations
 */
template <class t_real> static inline void rotate_around_axis(int axis, t_real &x, t_real &y, t_real &z, t_real angle) {
  // save the old angles
  const t_real _x = x;
  const t_real _y = y;
  const t_real _z = z;

  const t_real s = std::sin(angle);
  const t_real c = std::cos(angle);

  // rotate around x axis
  if (axis == 0) {
    y = _y * c - _z * s;
    z = _y * s + _z * c;
  }
  // rotate around y axis
  else if (axis == 1) {
    x = _x * c + _z * s;
    z = -_x * s + _z * c;
  }
  // rotate around z axis
  else if (axis == 2) {
    x = _x * c - _y * s;
    y = _x * s + _y * c;
  }
}

static int get_axis_index(std::int32_t dir_x, std::int32_t dir_y, std::int32_t dir_z) {
  int axis = 0;

  if (dir_x)
    axis = 0;
  else if (dir_y)
    axis = 1;
  else if (dir_z)
    axis = 2;

  return axis;
}
//----------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------
/** B matrix (https://en.wikipedia.org/wiki/Fractional_coordinates)
 * code from tlibs 2, doi: 10.5281/zenodo.5717779
 */
template <class t_real> static void B_matrix(const t_real *cell, const t_real *angles, t_real *B) {
  const t_real ca = std::cos(angles[0] / t_real(180. * M_PI));
  const t_real cb = std::cos(angles[1] / t_real(180. * M_PI));
  const t_real cc = std::cos(angles[2] / t_real(180. * M_PI));
  const t_real sc = std::sin(angles[2] / t_real(180. * M_PI));
  const t_real rr = std::sqrt(t_real(1) + t_real(2) * ca * cb * cc - (ca * ca + cb * cb + cc * cc));

  B[0] = t_real(1) / cell[0];
  B[1] = 0;
  B[2] = 0;
  B[3] = t_real(-1) / cell[0] * cc / sc;
  B[4] = 1.f / cell[1] * 1.f / sc;
  B[5] = 0;
  B[6] = (cc * ca - cb) / (cell[0] * sc * rr);
  B[7] = (cb * cc - ca) / (cell[1] * sc * rr);
  B[8] = sc / (cell[2] * rr);
}

/** position on spherical detector
 */
template <class t_real> static void get_spherical_pos(t_real theta, t_real phi, t_real *vec) {
  theta = -theta;
  theta += t_real(M_PI) * t_real(0.5);
  phi += t_real(M_PI) * t_real(0.5);

  vec[0] = std::cos(phi) * std::sin(theta);
  vec[1] = std::sin(phi) * std::sin(theta);
  vec[2] = std::cos(theta);
}

/** position on cylindrical detector
 */
template <class t_real> static void get_cylindrical_pos(t_real z, t_real phi, t_real rad, t_real *vec) {
  phi += t_real(M_PI) * t_real(0.5);

  vec[0] = rad * std::cos(phi);
  vec[1] = rad * std::sin(phi);
  vec[2] = z;

  t_real len = std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
  vec[0] /= len;
  vec[1] /= len;
  vec[2] /= len;
}
//----------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------
/** Initialize the algorithm's properties.
 */
void LoadILLSingleCrystal::init() {
  // input file
  declareProperty(std::make_unique<API::FileProperty>("Filename", "", API::FileProperty::Load, ".nxs"),
                  "File path of the data file to load");

  // create Q events workspace
  declareProperty(
      std::make_unique<Kernel::PropertyWithValue<bool>>("CreateEventWorkspace", false, Kernel::Direction::Input),
      "Create the workspace with Q events.");

  declareProperty(std::make_unique<Kernel::PropertyWithValue<std::size_t>>("QxBins", 64, Kernel::Direction::Input),
                  "Number of bins along the Qx direction.");
  declareProperty(std::make_unique<Kernel::PropertyWithValue<std::size_t>>("QyBins", 64, Kernel::Direction::Input),
                  "Number of bins along the Qy direction.");
  declareProperty(std::make_unique<Kernel::PropertyWithValue<std::size_t>>("QzBins", 64, Kernel::Direction::Input),
                  "Number of bins along the Qz direction.");

  // output workspaces
  // use a histogram workspace for pixel data
  declareProperty(std::make_unique<API::WorkspaceProperty<API::IMDHistoWorkspace>>(
                      "OutputHistoWorkspace", "out_histo", Kernel::Direction::Output, API::PropertyMode::Mandatory),
                  "An output histogram workspace for the detector pixel values.");
  // use an event workspace for Q coordinates
  declareProperty(std::make_unique<API::WorkspaceProperty<API::IMDEventWorkspace>>(
                      "OutputEventWorkspace", "out_event", Kernel::Direction::Output, API::PropertyMode::Optional),
                  "An output event workspace in momentum coordinates.");
}

//----------------------------------------------------------------------------------------------
/**
 * confidence that we can load the given file
 */
int LoadILLSingleCrystal::confidence(Kernel::NexusDescriptor &descr) const {
  bool det_data_exists = descr.pathExists("/entry0/data_scan/detector_data/data");

  return det_data_exists ? 75 : 0;
}

bool LoadILLSingleCrystal::LoadScannedVariables() {
  // --------------------------------------------------------------------------
  // hack: open the raw hdf5 file because not all data can be accessed
  // via the Nexus interface, e.g. string arrays.
  hid_t file_raw = H5Fopen(m_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_raw < 0)
    throw Kernel::Exception::FileError("Cannot open hdf5 file ", m_filename);

  m_scanned_variables_names = get_values(file_raw, "/entry0/data_scan/scanned_variables/variables_names/name");

  // the raw file has to be closed before reopening it using nexus
  H5Fclose(file_raw);
  return true;
  // --------------------------------------------------------------------------
}

bool LoadILLSingleCrystal::LoadInstrumentGroup() {
  namespace fs = boost::filesystem;

  m_file->openGroup("instrument", "NXinstrument");
  auto entries = m_file->getEntries();

  m_log.information() << "\nInstrument group:" << std::endl;
  for (const auto &entry : entries)
    m_log.information() << entry.first << " = " << entry.second << std::endl;

  // instrument name
  m_instr_name = *get_value<decltype(m_instr_name)>(*m_file, "name");

  // find the instrument definition file
  std::string instr_dir = Kernel::ConfigService::Instance().getInstrumentDirectory();
  fs::path instr_file = fs::path(instr_dir) / fs::path(m_instr_name + "_Definition.xml");
  m_instr_file = instr_file.string();
  if (!fs::exists(instr_file)) {
    m_log.error() << "No definition file was found for instrument \"" << m_instr_name << "\"" << std::endl;
  }

  // detector
  m_file->openGroup("Det1", "NXdetector");

  // detector dimensions
  m_det_num_rows = *get_value<decltype(m_det_num_rows)>(*m_file, "nrows");
  m_det_num_cols = *get_value<decltype(m_det_num_cols)>(*m_file, "ncols");

  // in mm
  t_real det_pixel_height = *get_value<t_real>(*m_file, "height");

  // there's a bug that just puts 0 as heigth value
  if (std::abs(det_pixel_height) < std::numeric_limits<t_real>::epsilon())
    det_pixel_height = 0.4f;

  // in mm
  m_det_height = det_pixel_height * t_real(m_det_num_rows);

  // in mm
  m_dist_sample_det = *get_value<decltype(m_dist_sample_det)>(*m_file, "sample_distance");

  // in deg
  m_det_angular_width = *get_value<decltype(m_det_angular_width)>(*m_file, "angular_width");
  m_det_angular_height = std::abs(std::atan(m_det_height / m_dist_sample_det));
  m_file->closeGroup();

  m_log.information() << "Detector: "
                      << "dimensions=" << m_det_num_rows << "x" << m_det_num_cols << ", "
                      << "angular width=" << m_det_angular_width << ", "
                      << "sample-detector distance=" << m_dist_sample_det << std::endl;

  m_det_angular_width = m_det_angular_width / 180.f * t_real(M_PI);

  // detector zero position
  m_file->openGroup("gamma", "NXpositioner");
  m_det_gamma_deg = *get_value<decltype(m_det_gamma_deg)>(*m_file, "value");
  std::int32_t gamma_cw = *get_value<std::int32_t>(*m_file, "clockwise");
  if (gamma_cw)
    m_det_gamma_deg = -m_det_gamma_deg;
  std::int32_t gamma_dir_x = *get_value<std::int32_t>(*m_file, "dir_x");
  std::int32_t gamma_dir_y = *get_value<std::int32_t>(*m_file, "dir_y");
  std::int32_t gamma_dir_z = *get_value<std::int32_t>(*m_file, "dir_z");
  m_file->closeGroup();

  // get detector rotation axis
  m_det_gamma_axis = get_axis_index(gamma_dir_x, gamma_dir_y, gamma_dir_z);

  m_log.information() << "Detector zero position: "
                      << "gamma=" << m_det_gamma_deg << ", axis: " << m_det_gamma_axis << std::endl;

  // initial crystal angles
  m_file->openGroup(SAMPLE_CHI_NAME, "NXpositioner");
  m_chi_deg = *get_value<decltype(m_chi_deg)>(*m_file, "value");
  std::int32_t chi_cw = *get_value<std::int32_t>(*m_file, "clockwise");
  if (chi_cw)
    m_chi_deg = -m_chi_deg;
  std::int32_t chi_dir_x = *get_value<std::int32_t>(*m_file, "dir_x");
  std::int32_t chi_dir_y = *get_value<std::int32_t>(*m_file, "dir_y");
  std::int32_t chi_dir_z = *get_value<std::int32_t>(*m_file, "dir_z");
  m_file->closeGroup();

  m_file->openGroup(SAMPLE_OMEGA_NAME, "NXpositioner");
  m_omega_deg = *get_value<decltype(m_omega_deg)>(*m_file, "value");
  std::int32_t omega_cw = *get_value<std::int32_t>(*m_file, "clockwise");
  if (omega_cw)
    m_omega_deg = -m_omega_deg;
  std::int32_t omega_dir_x = *get_value<std::int32_t>(*m_file, "dir_x");
  std::int32_t omega_dir_y = *get_value<std::int32_t>(*m_file, "dir_y");
  std::int32_t omega_dir_z = *get_value<std::int32_t>(*m_file, "dir_z");
  m_file->closeGroup();

  m_file->openGroup(SAMPLE_PHI_NAME, "NXpositioner");
  m_phi_deg = *get_value<decltype(m_phi_deg)>(*m_file, "value");
  std::int32_t phi_cw = *get_value<std::int32_t>(*m_file, "clockwise");
  if (phi_cw)
    m_phi_deg = -m_phi_deg;
  std::int32_t phi_dir_x = *get_value<std::int32_t>(*m_file, "dir_x");
  std::int32_t phi_dir_y = *get_value<std::int32_t>(*m_file, "dir_y");
  std::int32_t phi_dir_z = *get_value<std::int32_t>(*m_file, "dir_z");
  m_file->closeGroup();

  // get axis indices
  m_chi_axis = get_axis_index(chi_dir_x, chi_dir_y, chi_dir_z);
  m_omega_axis = get_axis_index(omega_dir_x, omega_dir_y, omega_dir_z);
  m_phi_axis = get_axis_index(phi_dir_x, phi_dir_y, phi_dir_z);

  m_log.information() << "Initial sample angles: "
                      << "chi=" << m_chi_deg << ", omega=" << m_omega_deg << ", phi=" << m_phi_deg << std::endl;
  m_log.information() << "Sample axes: "
                      << "chi=" << m_chi_axis << ", omega=" << m_omega_axis << ", phi=" << m_phi_axis << std::endl;

  // crystal
  m_file->openGroup("SingleCrystalSettings", "NXcrystal");
  auto _cell_a = get_value<t_real>(*m_file, "unit_cell_a");
  auto _cell_b = get_value<t_real>(*m_file, "unit_cell_b");
  auto _cell_c = get_value<t_real>(*m_file, "unit_cell_c");
  auto _cell_alpha = get_value<t_real>(*m_file, "unit_cell_alpha");
  auto _cell_beta = get_value<t_real>(*m_file, "unit_cell_beta");
  auto _cell_gamma = get_value<t_real>(*m_file, "unit_cell_gamma");

  // if a is not given, set it to 0
  m_cell[0] = _cell_a ? *_cell_a : 0;
  // if b is not given, set it to a
  m_cell[1] = _cell_b ? *_cell_b : m_cell[0];
  // if c is not given, set it to a
  m_cell[2] = _cell_c ? *_cell_c : m_cell[0];
  // if alpha is not given, set it to 90 deg
  m_cell_angles[0] = _cell_alpha ? *_cell_alpha : 90;
  // if beta is not given, set it to alpha
  m_cell_angles[1] = _cell_beta ? *_cell_beta : m_cell_angles[0];
  // if gamma is not given, set it to alpha
  m_cell_angles[2] = _cell_gamma ? *_cell_gamma : m_cell_angles[0];

  m_log.information() << "Unit cell: a=" << m_cell[0] << ", b=" << m_cell[1] << ", c=" << m_cell[2]
                      << ", alpha=" << m_cell_angles[0] << ", beta=" << m_cell_angles[1]
                      << ", gamma=" << m_cell_angles[2] << std::endl;

  m_UB = get_values<typename decltype(m_UB)::value_type>(*m_file, "orientation_matrix");

  // calculate B matrix from given lattice
  if (m_UB.size() != 9)
    m_UB.resize(9);
  // own calculation for verification
  // B_matrix(m_cell, m_cell_angles, m_UB.data());

  m_log.information() << "UB matrix:\n"
                      << "\t" << std::setw(15) << m_UB[0] << std::setw(15) << m_UB[1] << std::setw(15) << m_UB[2]
                      << "\n"
                      << "\t" << std::setw(15) << m_UB[3] << std::setw(15) << m_UB[4] << std::setw(15) << m_UB[5]
                      << "\n"
                      << "\t" << std::setw(15) << m_UB[6] << std::setw(15) << m_UB[7] << std::setw(15) << m_UB[8]
                      << std::endl;

  m_file->closeGroup(); // SingleCrystalSettings
  m_file->closeGroup(); // instrument

  return true;
}

bool LoadILLSingleCrystal::LoadSampleGroup() {
  m_file->openGroup("sample", "NXsample");
  auto entries = m_file->getEntries();

  m_log.information() << "\nSample group:" << std::endl;
  for (const auto &entry : entries)
    m_log.information() << entry.first << " = " << entry.second << std::endl;

  m_file->closeGroup(); // sample

  return true;
}

bool LoadILLSingleCrystal::LoadUserGroup() {
  m_file->openGroup("user", "NXuser");
  auto entries = m_file->getEntries();

  m_log.information() << "\nUser group:" << std::endl;
  for (const auto &entry : entries)
    m_log.information() << entry.first << " = " << entry.second << std::endl;

  m_file->closeGroup(); // user

  return true;
}

bool LoadILLSingleCrystal::LoadMonitorGroup() {
  m_file->openGroup("monitor", "NXmonitor");
  auto entries = m_file->getEntries();

  m_log.information() << "\nMonitor group:" << std::endl;
  for (const auto &entry : entries)
    m_log.information() << entry.first << " = " << entry.second << std::endl;

  // total monitor counts
  m_monitor_sum = *get_value<decltype(m_monitor_sum)>(*m_file, "monsum");

  m_file->closeGroup(); // monitor

  return true;
}

bool LoadILLSingleCrystal::LoadDataScanGroup() {
  m_file->openGroup("data_scan", "NXdata");
  auto entries = m_file->getEntries();

  m_log.information() << "\nData group:" << std::endl;
  for (const auto &entry : entries)
    m_log.information() << entry.first << " = " << entry.second << std::endl;

  // steps
  std::int64_t total_steps = *get_value<std::int64_t>(*m_file, "total_steps");
  std::int64_t actual_step = *get_value<std::int64_t>(*m_file, "actual_step");

  m_log.information() << "Total steps=" << total_steps << ", actual step=" << actual_step << std::endl;

  // detector data
  m_file->openGroup("detector_data", "ILL_detectors_data_scan");

  m_detector_data_dims = get_dims(*m_file, "data");
  m_detector_data = get_values<typename decltype(m_detector_data)::value_type>(*m_file, "data");

  if (m_detector_data_dims.size() < 3)
    throw std::runtime_error("Detector data dimension < 3.");

  m_log.information() << "Detector data dimensions: ";
  for (std::size_t i = 0; i < m_detector_data_dims.size(); ++i) {
    m_log.information() << m_detector_data_dims[i];
    if (i < m_detector_data_dims.size() - 1)
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

  m_scanned_variables_values = get_values<typename decltype(m_scanned_variables_values)::value_type>(*m_file, "data");

  m_file->openGroup("variables_names", "ILL_data_scan_vars");
  auto variables_names_entries = m_file->getEntries();

  m_log.information() << "\nData group -> Scanned variables -> Variables names:" << std::endl;
  for (const auto &entry : variables_names_entries)
    m_log.information() << entry.first << " = " << entry.second << std::endl;

  // already got m_scanned_variables_names via raw hdf5 interface
  // auto_m_scanned_variables_names = get_values<std::string>(*m_file, "name");
  auto scanned_variables = get_values<std::uint8_t>(*m_file, "scanned");

  // get indices of scanned variables
  m_first_scanned_idx = 0;
  m_chi_idx = 0;
  m_omega_idx = 0;
  m_phi_idx = 0;

  m_first_scanned_idx_set = false;
  m_chi_idx_set = false;
  m_omega_idx_set = false;
  m_phi_idx_set = false;

  for (std::size_t idx = 0; idx < scanned_variables.size(); ++idx) {
    if (scanned_variables[idx]) {
      if (!m_first_scanned_idx_set) {
        m_first_scanned_idx = idx;
        m_first_scanned_idx_set = true;
      }

      if (!m_chi_idx_set && m_scanned_variables_names[idx] == SAMPLE_CHI_NAME) {
        m_chi_idx = idx;
        m_chi_idx_set = true;
      }

      if (!m_omega_idx_set && m_scanned_variables_names[idx] == SAMPLE_OMEGA_NAME) {
        m_omega_idx = idx;
        m_omega_idx_set = true;
      }

      if (!m_phi_idx_set && m_scanned_variables_names[idx] == SAMPLE_PHI_NAME) {
        m_phi_idx = idx;
        m_phi_idx = true;
      }
    }
  }

  m_log.information() << "First scanned variable: " << m_scanned_variables_names[m_first_scanned_idx]
                      << ", index: " << m_first_scanned_idx << std::endl;

  m_file->closeGroup(); // variable names
  m_file->closeGroup(); // scanned_variables
  m_file->closeGroup(); // data_scan

  return true;
}

//----------------------------------------------------------------------------------------------
/** Execute the algorithm.
 */
void LoadILLSingleCrystal::exec() {
  // data type for event workspace
  using t_event = DataObjects::MDEvent<3>;
  using t_eventworkspace = DataObjects::MDEventWorkspace<t_event, 3>;

  bool create_event_workspace = getProperty("CreateEventWorkspace");
  if (getProperty("OutputEventWorkspace").operator std::string() == "")
    create_event_workspace = false;

  // get input file name
  m_filename = getPropertyValue("Filename");

  LoadScannedVariables();

  // open the corresponding nexus file
  m_file = std::make_unique<NeXus::File>(m_filename, NXACC_READ);
  if (!m_file)
    throw Kernel::Exception::FileError("Cannot open Nexus file ", m_filename);

  // --------------------------------------------------------------------------
  // load root group
  m_file->openGroup("/entry0", "NXentry");
  auto entries = m_file->getEntries();

  m_log.information() << "\nRoot entry group:" << std::endl;
  for (const auto &entry : entries)
    m_log.information() << entry.first << " = " << entry.second << std::endl;

  // title
  std::string title = *get_value<std::string>(*m_file, "title");

  // run number
  std::int64_t numor = *get_value<std::int64_t>(*m_file, "run_number");

  // wavelength
  t_real wavelength = *get_value<t_real>(*m_file, "wavelength");
  t_real wavenumber = t_real(2. * M_PI) / wavelength;
  m_log.information() << "Wavelength: " << wavelength << ", wavenumber: " << wavenumber << std::endl;

  LoadInstrumentGroup();
  LoadSampleGroup();
  LoadUserGroup();
  LoadMonitorGroup();
  LoadDataScanGroup();

  m_file->closeGroup(); // entry0
  // --------------------------------------------------------------------------

  // create histogram workspace
  m_log.information() << "\nCreating histogram workspace..." << std::endl;

  Geometry::GeneralFrame frame_x("x", "");
  Geometry::GeneralFrame frame_y("y", "");
  Geometry::GeneralFrame frame_scanned("scanned", "");

  std::vector<Mantid::Geometry::MDHistoDimension_sptr> dimensions{{
      std::make_shared<Geometry::MDHistoDimension>("y", "y", frame_y, 0, m_detector_data_dims[2] - 1,
                                                   m_detector_data_dims[2]),
      std::make_shared<Geometry::MDHistoDimension>("x", "x", frame_x, 0, m_detector_data_dims[1] - 1,
                                                   m_detector_data_dims[1]),
      std::make_shared<Geometry::MDHistoDimension>(m_scanned_variables_names[m_first_scanned_idx],
                                                   m_scanned_variables_names[m_first_scanned_idx], frame_scanned, 0,
                                                   m_detector_data_dims[0] - 1, m_detector_data_dims[0]),
  }};

  m_workspace_histo = std::make_shared<DataObjects::MDHistoWorkspace>(dimensions, API::NoNormalization);

  // load instrument definition file into a dummy workspace
  auto workspace_instr = API::WorkspaceFactory::Instance().create("Workspace2D", 1, 1, 1);
  if (m_instr_file != "") {
    API::Algorithm_sptr load_instr = createChildAlgorithm("LoadInstrument");

    if (load_instr) {
      load_instr->setPropertyValue("Filename", m_instr_file);
      load_instr->setProperty("Workspace", workspace_instr);
      load_instr->setProperty<Kernel::OptionalBool>("RewriteSpectraMap", true);

      if (load_instr->execute())
        m_log.information() << "Successfully loaded instrument parameters from \"" << m_instr_file << "\"" << std::endl;
      else
        m_log.error() << "Failed to loaded instrument parameters from \"" << m_instr_file << "\"" << std::endl;

      // TODO: load detector dimensions, etc.
    }
  }

  // set the metadata
  auto info = std::make_shared<API::ExperimentInfo>();
  info->mutableRun().addProperty("Filename", m_filename);
  info->mutableRun().addProperty("Instrument", m_instr_name);
  info->mutableRun().addProperty("Wavelength", wavelength);
  info->mutableRun().addProperty("Wavenumber", wavenumber);
  info->mutableRun().addProperty("Numor", numor);
  info->mutableRun().addProperty("Monitor", m_monitor_sum);
  info->mutableRun().addProperty("Sample_a", m_cell[0]);
  info->mutableRun().addProperty("Sample_b", m_cell[1]);
  info->mutableRun().addProperty("Sample_c", m_cell[2]);
  info->mutableRun().addProperty("Sample_alpha", m_cell_angles[0]);
  info->mutableRun().addProperty("Sample_beta", m_cell_angles[1]);
  info->mutableRun().addProperty("Sample_gamma", m_cell_angles[2]);

  // set the metadata units
  info->mutableRun().getProperty("Wavelength")->setUnits("Angstrom");
  info->mutableRun().getProperty("Wavenumber")->setUnits("1/Angstrom");
  info->mutableRun().getProperty("Sample_a")->setUnits("Angstrom");
  info->mutableRun().getProperty("Sample_b")->setUnits("Angstrom");
  info->mutableRun().getProperty("Sample_c")->setUnits("Angstrom");
  info->mutableRun().getProperty("Sample_alpha")->setUnits("degree");
  info->mutableRun().getProperty("Sample_beta")->setUnits("degree");
  info->mutableRun().getProperty("Sample_gamma")->setUnits("degree");

  // set the loaded instrument infos, needed by PredictPeaks.cpp
  auto instr = workspace_instr->getInstrument();
  info->setInstrument(instr);

  for (std::size_t i = 0; i < m_UB.size(); ++i)
    info->mutableRun().addProperty("Sample_UB_" + std::to_string(i), m_UB[i]);

  Geometry::CrystalStructure crys(
      // unit cell
      Geometry::UnitCell(m_cell[0], m_cell[1], m_cell[2], m_cell_angles[0], m_cell_angles[1], m_cell_angles[2],
                         Geometry::angDegrees),
      // dummy space group
      Geometry::SpaceGroupFactory::Instance().createSpaceGroup("P 1"),
      // dummy scatterer
      Geometry::CompositeBraggScatterer::create(Geometry::IsotropicAtomBraggScattererParser("")()));
  info->mutableSample().setCrystalStructure(crys);

  auto lattice = std::make_unique<Geometry::OrientedLattice>(m_cell[0], m_cell[1], m_cell[2], m_cell_angles[0],
                                                             m_cell_angles[1], m_cell_angles[2]);
  info->mutableSample().setOrientedLattice(std::move(lattice));

  m_workspace_histo->addExperimentInfo(info);
  m_workspace_histo->setTitle(title);

  // minimum and maximum Qs
  t_real Qxminmax[2] = {std::numeric_limits<t_real>::max(), std::numeric_limits<t_real>::lowest()};
  t_real Qyminmax[2] = {std::numeric_limits<t_real>::max(), std::numeric_limits<t_real>::lowest()};
  t_real Qzminmax[2] = {std::numeric_limits<t_real>::max(), std::numeric_limits<t_real>::lowest()};

  // size of one detector image in pixels
  std::int64_t frame_size = m_detector_data_dims[1] * m_detector_data_dims[2];

  std::vector<t_event> events;
  events.reserve(m_detector_data.size());

  // copy the detector data to the workspace
  std::int64_t last_scan_step = -1;
  API::Progress progress(this, 0., create_event_workspace ? 0.5 : 1., m_detector_data.size() / frame_size);
  bool cancelled = false;

  for (std::size_t idx = 0; idx < m_detector_data.size(); ++idx) {
    // completed all pixels of a frame?
    std::int64_t cur_scan_step = idx / frame_size;
    std::int64_t cur_frame_step = idx % frame_size;
    std::int64_t cur_y = cur_frame_step / m_detector_data_dims[1];
    std::int64_t cur_x = cur_frame_step % m_detector_data_dims[1];

    // save pixels in histogram workspace
    Mantid::signal_t intensity = m_detector_data[idx];
    m_workspace_histo->setSignalAt(idx, intensity);

    // only update this once per frame
    if (last_scan_step != cur_scan_step) {
      // get sample angles
      if (m_chi_idx_set)
        m_chi_deg = t_real(m_scanned_variables_values[m_detector_data_dims[0] * m_chi_idx + cur_scan_step]);
      if (m_omega_idx_set)
        m_omega_deg = t_real(m_scanned_variables_values[m_detector_data_dims[0] * m_omega_idx + cur_scan_step]);
      if (m_phi_idx_set)
        m_phi_deg = t_real(m_scanned_variables_values[m_detector_data_dims[0] * m_phi_idx + cur_scan_step]);

      last_scan_step = cur_scan_step;

      progress.report(cur_scan_step);
      if (progress.hasCancellationBeenRequested()) {
        cancelled = true;
        break;
      }
    }

    t_real chi = t_real(m_chi_deg * M_PI / 180.);
    t_real omega = t_real(m_omega_deg * M_PI / 180.);
    t_real phi = t_real(m_phi_deg * M_PI / 180.);

    if (create_event_workspace) {
      t_real pix_coord[2] = {t_real(cur_x), t_real(cur_y)};

      // in-plane scattering angle
      t_real twotheta = pix_coord[0] / t_real(m_detector_data_dims[0]) * m_det_angular_width;
      // rotate to detector zero position
      // twotheta -= m_det_angular_width * 0.5f;
      twotheta += m_det_gamma_deg / t_real(180. * M_PI);

      // out-of-plane scattering angle
      t_real angle_oop = pix_coord[1] / t_real(m_detector_data_dims[1]) * m_det_angular_height;
      angle_oop -= m_det_angular_height * 0.5f;

      t_real pos_oop = pix_coord[1] / t_real(m_detector_data_dims[1]) * m_det_height;
      pos_oop -= m_det_height * 0.5f;

      // convert pixels to Q coordinates and insert events
      t_real ki[3] = {0, 1, 0};
      t_real kf[3] = {0, 1, 0};
      // get_spherical_pos(angle_oop, twotheta, kf);
      get_cylindrical_pos(pos_oop, twotheta, m_dist_sample_det, kf);
      t_real Q_coord[3] = {(ki[0] - kf[0]) * wavenumber, (ki[1] - kf[1]) * wavenumber, (ki[2] - kf[2]) * wavenumber};

      // rotate by sample euler angles
      // TODO: check order
      rotate_around_axis(m_omega_axis, Q_coord[0], Q_coord[1], Q_coord[2], omega);
      rotate_around_axis(m_chi_axis, Q_coord[0], Q_coord[1], Q_coord[2], chi);
      rotate_around_axis(m_phi_axis, Q_coord[0], Q_coord[1], Q_coord[2], phi);

      // calculate Q ranges
      Qxminmax[0] = std::min(Q_coord[0], Qxminmax[0]);
      Qxminmax[1] = std::max(Q_coord[0], Qxminmax[1]);
      Qyminmax[0] = std::min(Q_coord[1], Qyminmax[0]);
      Qyminmax[1] = std::max(Q_coord[1], Qyminmax[1]);
      Qzminmax[0] = std::min(Q_coord[2], Qzminmax[0]);
      Qzminmax[1] = std::max(Q_coord[2], Qzminmax[1]);

      // store all events
      events.emplace_back(t_event(intensity, std::sqrt(intensity), Q_coord));
    }
  }

  // create event workspace
  if (create_event_workspace && !cancelled) {
    m_log.information() << "\nCreating event workspace..." << std::endl;

    m_log.information() << "Q_x range: [" << Qxminmax[0] << ", " << Qxminmax[1] << "]." << std::endl;
    m_log.information() << "Q_y range: [" << Qyminmax[0] << ", " << Qyminmax[1] << "]." << std::endl;
    m_log.information() << "Q_z range: [" << Qzminmax[0] << ", " << Qzminmax[1] << "]." << std::endl;

    std::size_t progress_div = events.size() / 100;
    API::Progress progress_evt(this, 0.5, 1., events.size() / progress_div);

    std::size_t num_Qx_bins = getProperty("QxBins");
    std::size_t num_Qy_bins = getProperty("QyBins");
    std::size_t num_Qz_bins = getProperty("QzBins");

    Geometry::QLab frame_Qx, frame_Qy, frame_Qz;

    std::vector<Mantid::Geometry::MDHistoDimension_sptr> dimensions_Q{{
        std::make_shared<Geometry::MDHistoDimension>("Qx", "Qx", frame_Qx, Qxminmax[0], Qxminmax[1], num_Qx_bins),
        std::make_shared<Geometry::MDHistoDimension>("Qy", "Qy", frame_Qy, Qyminmax[0], Qyminmax[1], num_Qy_bins),
        std::make_shared<Geometry::MDHistoDimension>("Qz", "Qz", frame_Qz, Qzminmax[0], Qzminmax[1], num_Qz_bins),
    }};

    m_workspace_event = std::make_shared<t_eventworkspace>(API::NoNormalization, API::NoNormalization);

    m_workspace_event->addDimension(dimensions_Q[0]);
    m_workspace_event->addDimension(dimensions_Q[1]);
    m_workspace_event->addDimension(dimensions_Q[2]);
    m_workspace_event->setCoordinateSystem(Kernel::QLab);
    m_workspace_event->initialize();

    m_workspace_event->addExperimentInfo(info);
    m_workspace_event->setTitle(title);

    // add all events
    for (std::size_t evt_idx = 0; evt_idx < events.size(); ++evt_idx) {
      if (evt_idx % progress_div == 0) {
        progress_evt.report(evt_idx / progress_div);
        if (progress_evt.hasCancellationBeenRequested()) {
          cancelled = true;
          break;
        }
      }

      const t_event &event = events[evt_idx];
      std::dynamic_pointer_cast<t_eventworkspace>(m_workspace_event)->addEvent(event);
    }

    m_workspace_event->refreshCache();
  }

  if (!cancelled) {
    // set ouput workspaces
    if (m_workspace_histo) {
      setProperty("OutputHistoWorkspace", std::dynamic_pointer_cast<API::IMDHistoWorkspace>(m_workspace_histo));
    }

    if (m_workspace_event) {
      setProperty("OutputEventWorkspace", std::dynamic_pointer_cast<API::IMDEventWorkspace>(m_workspace_event));
    }
  }
}

} // namespace DataHandling
} // namespace Mantid

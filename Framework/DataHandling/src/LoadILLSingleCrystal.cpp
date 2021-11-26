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
#include <hdf5.h>

// sample angle names
#define SAMPLE_CHI_NAME "chi"
#define SAMPLE_OMEGA_NAME "omega"
#define SAMPLE_PHI_NAME "phi"

/*
 @class LoadILLSingleCrystal
 @author T. Weber, ILL
 @date 9/2021
 */
namespace Mantid {
namespace DataHandling {

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
  declareProperty(std::make_unique<API::WorkspaceProperty<API::IMDHistoWorkspace>>(
                      "OutputHistogramWorkspace", "out_histo", Kernel::Direction::Output, API::PropertyMode::Mandatory),
                  "An output histogram workspace for the detector pixel values.");
  // we do not really use events, but some of the following algorithms expect an event workspace
  declareProperty(std::make_unique<API::WorkspaceProperty<API::IMDEventWorkspace>>(
                      "OutputEventWorkspace", "out_event", Kernel::Direction::Output, API::PropertyMode::Optional),
                  "An output event workspace in momentum coordinates.");
}

//----------------------------------------------------------------------------------------------
/** Execute the algorithm.
 */
void LoadILLSingleCrystal::exec() {
  // number data types
  using t_real = coord_t;
  using t_data = double;

  // data type for event workspace
  using t_event = DataObjects::MDEvent<3>;
  using t_eventworkspace = DataObjects::MDEventWorkspace<t_event, 3>;

  bool create_event_workspace = getProperty("CreateEventWorkspace");
  if (getProperty("OutputEventWorkspace").operator std::string() == "")
    create_event_workspace = false;

  // get input file name
  m_filename = getPropertyValue("Filename");

  // --------------------------------------------------------------------------
  // hack: open the raw hdf5 file because not all data can be accessed
  // via the Nexus interface, e.g. string arrays.
  hid_t file_raw = H5Fopen(m_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_raw < 0)
    throw Kernel::Exception::FileError("Cannot open hdf5 file ", m_filename);

  auto scanned_variables_names = get_values(file_raw, "/entry0/data_scan/scanned_variables/variables_names/name");

  // the raw file has to be closed before reopening it using nexus
  H5Fclose(file_raw);
  // --------------------------------------------------------------------------

  // open the corresponding nexus file
  m_file = std::make_unique<NeXus::File>(m_filename, NXACC_READ);
  if (!m_file)
    throw Kernel::Exception::FileError("Cannot open Nexus file ", m_filename);

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
  std::int64_t numor = *get_value<std::int64_t>(*m_file, "run_number");

  // wavelength
  t_real wavelength = *get_value<t_real>(*m_file, "wavelength");
  t_real wavenumber = t_real(2. * M_PI) / wavelength;
  m_log.information() << "Wavelength: " << wavelength << std::endl;
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

  // detector
  m_file->openGroup("Det1", "NXdetector");
  // in mm
  t_real det_height = 441; // TODO: get this from the file
  // in mm
  t_real dist_sample_det = *get_value<t_real>(*m_file, "sample_distance");
  // in deg
  t_real det_angular_width = *get_value<t_real>(*m_file, "angular_width");
  t_real det_angular_height = std::abs(std::atan(det_height / dist_sample_det));
  m_file->closeGroup();

  m_log.information() << "Detector: "
                      << "angular width=" << det_angular_width << ", sample-detector distance=" << dist_sample_det
                      << std::endl;

  det_angular_width = det_angular_width / 180.f * t_real(M_PI);

  // initial crystal angles
  m_file->openGroup(SAMPLE_CHI_NAME, "NXpositioner");
  t_real chi_deg = *get_value<t_real>(*m_file, "value");
  std::int32_t chi_dir_x = *get_value<std::int32_t>(*m_file, "dir_x");
  std::int32_t chi_dir_y = *get_value<std::int32_t>(*m_file, "dir_y");
  std::int32_t chi_dir_z = *get_value<std::int32_t>(*m_file, "dir_z");
  m_file->closeGroup();
  m_file->openGroup(SAMPLE_OMEGA_NAME, "NXpositioner");
  t_real omega_deg = *get_value<t_real>(*m_file, "value");
  std::int32_t omega_dir_x = *get_value<std::int32_t>(*m_file, "dir_x");
  std::int32_t omega_dir_y = *get_value<std::int32_t>(*m_file, "dir_y");
  std::int32_t omega_dir_z = *get_value<std::int32_t>(*m_file, "dir_z");
  m_file->closeGroup();
  m_file->openGroup(SAMPLE_PHI_NAME, "NXpositioner");
  t_real phi_deg = *get_value<t_real>(*m_file, "value");
  std::int32_t phi_dir_x = *get_value<std::int32_t>(*m_file, "dir_x");
  std::int32_t phi_dir_y = *get_value<std::int32_t>(*m_file, "dir_y");
  std::int32_t phi_dir_z = *get_value<std::int32_t>(*m_file, "dir_z");
  m_file->closeGroup();

  // get axis indices
  int chi_axis = get_axis_index(chi_dir_x, chi_dir_y, chi_dir_z);
  int omega_axis = get_axis_index(omega_dir_x, omega_dir_y, omega_dir_z);
  int phi_axis = get_axis_index(phi_dir_x, phi_dir_y, phi_dir_z);

  m_log.information() << "Initial sample angles: "
                      << "chi=" << chi_deg << ", omega=" << omega_deg << ", phi=" << phi_deg << std::endl;
  m_log.information() << "Sample axes: "
                      << "chi=" << chi_axis << ", omega=" << omega_axis << ", phi=" << phi_axis << std::endl;

  // crystal
  m_file->openGroup("SingleCrystalSettings", "NXcrystal");
  auto _cell_a = get_value<t_real>(*m_file, "unit_cell_a");
  auto _cell_b = get_value<t_real>(*m_file, "unit_cell_b");
  auto _cell_c = get_value<t_real>(*m_file, "unit_cell_c");
  auto _cell_alpha = get_value<t_real>(*m_file, "unit_cell_alpha");
  auto _cell_beta = get_value<t_real>(*m_file, "unit_cell_beta");
  auto _cell_gamma = get_value<t_real>(*m_file, "unit_cell_gamma");

  // if a is not given, set it to 0
  t_real cell_a = _cell_a ? *_cell_a : 0;
  // if b is not given, set it to a
  t_real cell_b = _cell_b ? *_cell_b : cell_a;
  // if c is not given, set it to a
  t_real cell_c = _cell_c ? *_cell_c : cell_a;
  // if alpha is not given, set it to 90 deg
  t_real cell_alpha = _cell_alpha ? *_cell_alpha : 90;
  // if beta is not given, set it to alpha
  t_real cell_beta = _cell_beta ? *_cell_beta : cell_alpha;
  // if gamma is not given, set it to alpha
  t_real cell_gamma = _cell_gamma ? *_cell_gamma : cell_alpha;

  m_log.information() << "Unit cell: a=" << cell_a << ", b=" << cell_b << ", c=" << cell_c << ", alpha=" << cell_alpha
                      << ", beta=" << cell_beta << ", gamma=" << cell_gamma << std::endl;

  auto UB = get_values<t_real>(*m_file, "orientation_matrix");
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
  t_real monitor_sum = *get_value<t_real>(*m_file, "monsum");

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
  std::int64_t total_steps = *get_value<std::int64_t>(*m_file, "total_steps");
  std::int64_t actual_step = *get_value<std::int64_t>(*m_file, "actual_step");

  m_log.information() << "Total steps=" << total_steps << ", actual step=" << actual_step << std::endl;

  // detector data
  m_file->openGroup("detector_data", "ILL_detectors_data_scan");

  auto detector_data_dims = get_dims(*m_file, "data");
  auto detector_data = get_values<std::uint32_t>(*m_file, "data");

  if (detector_data_dims.size() < 3)
    throw std::runtime_error("Detector data dimension < 3.");

  m_log.information() << "Detector data dimensions: ";
  for (std::size_t i = 0; i < detector_data_dims.size(); ++i) {
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

  auto scanned_variables_values = get_values<t_data>(*m_file, "data");

  m_file->openGroup("variables_names", "ILL_data_scan_vars");
  auto variables_names_entries = m_file->getEntries();

  m_log.information() << "\nData group -> Scanned variables -> Variables names:" << std::endl;
  for (const auto &entry : variables_names_entries)
    m_log.information() << entry.first << " = " << entry.second << std::endl;

  // already got scanned_variables_names via raw hdf5 interface
  // auto scanned_variables_names = get_values<std::string>(*m_file, "name");
  auto scanned_variables = get_values<std::uint8_t>(*m_file, "scanned");

  // get indices of scanned variables
  std::size_t first_scanned_idx = 0;
  std::size_t chi_idx = 0;
  std::size_t omega_idx = 0;
  std::size_t phi_idx = 0;

  bool first_scanned_idx_set = false;
  bool chi_idx_set = false;
  bool omega_idx_set = false;
  bool phi_idx_set = false;

  for (std::size_t idx = 0; idx < scanned_variables.size(); ++idx) {
    if (scanned_variables[idx]) {
      if (!first_scanned_idx_set) {
        first_scanned_idx = idx;
        first_scanned_idx_set = true;
      }

      if (!chi_idx_set && scanned_variables_names[idx] == SAMPLE_CHI_NAME) {
        chi_idx = idx;
        chi_idx_set = true;
      }

      if (!omega_idx_set && scanned_variables_names[idx] == SAMPLE_OMEGA_NAME) {
        omega_idx = idx;
        omega_idx_set = true;
      }

      if (!phi_idx_set && scanned_variables_names[idx] == SAMPLE_PHI_NAME) {
        phi_idx = idx;
        phi_idx = true;
      }
    }
  }

  m_log.information() << "First scanned variable: " << scanned_variables_names[first_scanned_idx]
                      << ", index: " << first_scanned_idx << std::endl;

  m_file->closeGroup(); // variable names
  m_file->closeGroup(); // scanned_variables
  m_file->closeGroup(); // data_scan
  // --------------------------------------------------------------------------

  m_file->closeGroup();

  // create histogram workspace
  Geometry::GeneralFrame frame_x("x", "");
  Geometry::GeneralFrame frame_y("y", "");
  Geometry::GeneralFrame frame_scanned("scanned", "");

  std::vector<Mantid::Geometry::MDHistoDimension_sptr> dimensions{{
      std::make_shared<Geometry::MDHistoDimension>("y", "y", frame_y, 0, detector_data_dims[2] - 1,
                                                   detector_data_dims[2]),
      std::make_shared<Geometry::MDHistoDimension>("x", "x", frame_x, 0, detector_data_dims[1] - 1,
                                                   detector_data_dims[1]),
      std::make_shared<Geometry::MDHistoDimension>(scanned_variables_names[first_scanned_idx],
                                                   scanned_variables_names[first_scanned_idx], frame_scanned, 0,
                                                   detector_data_dims[0] - 1, detector_data_dims[0]),
  }};

  m_workspace_histo = std::make_shared<DataObjects::MDHistoWorkspace>(dimensions, API::NoNormalization);

  // set the metadata
  auto info = std::make_shared<API::ExperimentInfo>();
  info->mutableRun().addProperty("Filename", m_filename);
  info->mutableRun().addProperty("Instrument", instr_name);
  info->mutableRun().addProperty("Wavelength", wavelength);
  info->mutableRun().addProperty("Wavenumber", wavenumber);
  info->mutableRun().addProperty("Numor", numor);
  info->mutableRun().addProperty("Monitor", monitor_sum);
  info->mutableRun().addProperty("Sample_a", cell_a);
  info->mutableRun().addProperty("Sample_b", cell_b);
  info->mutableRun().addProperty("Sample_c", cell_c);
  info->mutableRun().addProperty("Sample_alpha", cell_alpha);
  info->mutableRun().addProperty("Sample_beta", cell_beta);
  info->mutableRun().addProperty("Sample_gamma", cell_gamma);

  // set the metadata units
  info->mutableRun().getProperty("Wavelength")->setUnits("Angstrom");
  info->mutableRun().getProperty("Wavenumber")->setUnits("1/Angstrom");
  info->mutableRun().getProperty("Sample_a")->setUnits("Angstrom");
  info->mutableRun().getProperty("Sample_b")->setUnits("Angstrom");
  info->mutableRun().getProperty("Sample_c")->setUnits("Angstrom");
  info->mutableRun().getProperty("Sample_alpha")->setUnits("degree");
  info->mutableRun().getProperty("Sample_beta")->setUnits("degree");
  info->mutableRun().getProperty("Sample_gamma")->setUnits("degree");

  // set the instrument with a sample and source, needed by PredictPeaks.cpp
  auto instr = std::make_shared<Geometry::Instrument>();

  // add a sample
  auto *sample = new Geometry::Component("sample", instr.get());
  sample->setPos(0, 0, 0); // TODO
  instr->add(sample);
  instr->markAsSamplePos(sample);

  // add a source
  auto *source = new Geometry::Component("source", instr.get());
  source->setPos(0, 0, -10); // TODO
  instr->add(source);
  instr->markAsSource(source);

  info->setInstrument(instr);

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

  auto lattice = std::make_unique<Geometry::OrientedLattice>(cell_a, cell_b, cell_c, cell_alpha, cell_beta, cell_gamma);
  info->mutableSample().setOrientedLattice(std::move(lattice));

  m_workspace_histo->addExperimentInfo(info);
  m_workspace_histo->setTitle(title);

  // minimum and maximum Qs
  t_real Qxminmax[2] = {std::numeric_limits<t_real>::max(), std::numeric_limits<t_real>::lowest()};
  t_real Qyminmax[2] = {std::numeric_limits<t_real>::max(), std::numeric_limits<t_real>::lowest()};
  t_real Qzminmax[2] = {std::numeric_limits<t_real>::max(), std::numeric_limits<t_real>::lowest()};

  // size of one detector image in pixels
  std::int64_t frame_size = detector_data_dims[1] * detector_data_dims[2];

  std::vector<t_event> events;
  events.reserve(detector_data.size());

  // copy the detector data to the workspace
  std::int64_t last_scan_step = -1;
  API::Progress progress(this, 0., create_event_workspace ? 0.5 : 1., detector_data.size() / frame_size);
  bool cancelled = false;

  for (std::size_t idx = 0; idx < detector_data.size(); ++idx) {
    // completed all pixels of a frame?
    std::int64_t cur_scan_step = idx / frame_size;
    std::int64_t cur_frame_step = idx % frame_size;
    std::int64_t cur_y = cur_frame_step / detector_data_dims[1];
    std::int64_t cur_x = cur_frame_step % detector_data_dims[1];

    // save pixels in histogram workspace
    Mantid::signal_t intensity = detector_data[idx];
    m_workspace_histo->setSignalAt(idx, intensity);

    // only update this once per frame
    if (last_scan_step != cur_scan_step) {
      // get sample angles
      if (chi_idx_set)
        chi_deg = t_real(scanned_variables_values[detector_data_dims[0] * chi_idx + cur_scan_step]);
      if (omega_idx_set)
        omega_deg = t_real(scanned_variables_values[detector_data_dims[0] * omega_idx + cur_scan_step]);
      if (phi_idx_set)
        phi_deg = t_real(scanned_variables_values[detector_data_dims[0] * phi_idx + cur_scan_step]);

      last_scan_step = cur_scan_step;

      progress.report(cur_scan_step);
      if (progress.hasCancellationBeenRequested()) {
        cancelled = true;
        break;
      }
    }

    t_real chi = t_real(chi_deg * M_PI / 180.);
    t_real omega = t_real(omega_deg * M_PI / 180.);
    t_real phi = t_real(phi_deg * M_PI / 180.);

    if (create_event_workspace) {
      t_real pix_coord[2] = {t_real(cur_x), t_real(cur_y)};

      // TODO: find zero positions on detector
      // in-plane scattering angle
      t_real twotheta = pix_coord[0] / t_real(detector_data_dims[0]) * det_angular_width;
      twotheta -= det_angular_width * 0.5f;

      // out-of-plane scattering angle
      t_real angle_oop = pix_coord[1] / t_real(detector_data_dims[1]) * det_angular_height;
      angle_oop -= det_angular_height * 0.5f;

      // convert pixels to Q coordinates and insert events
      t_real Q_coord[3] = {
          wavenumber * std::sin(angle_oop) * std::cos(twotheta),
          wavenumber * std::sin(angle_oop) * std::sin(twotheta),
          wavenumber * std::cos(angle_oop),
      };

      // rotate by sample euler angles
      // TODO: check order
      rotate_around_axis(omega_axis, Q_coord[0], Q_coord[1], Q_coord[2], omega);
      rotate_around_axis(chi_axis, Q_coord[0], Q_coord[1], Q_coord[2], chi);
      rotate_around_axis(phi_axis, Q_coord[0], Q_coord[1], Q_coord[2], phi);

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

  m_log.information() << "Q_x range: [" << Qxminmax[0] << ", " << Qxminmax[1] << "]." << std::endl;
  m_log.information() << "Q_y range: [" << Qyminmax[0] << ", " << Qyminmax[1] << "]." << std::endl;
  m_log.information() << "Q_z range: [" << Qzminmax[0] << ", " << Qzminmax[1] << "]." << std::endl;

  // create event workspace
  if (create_event_workspace && !cancelled) {
    std::size_t progress_div = events.size() / 100;
    API::Progress progress_evt(this, 0.5, 1., events.size() / progress_div);

    std::size_t num_Qx_bins = getProperty("QxBins");
    std::size_t num_Qy_bins = getProperty("QyBins");
    std::size_t num_Qz_bins = getProperty("QzBins");

    Geometry::QSample frame_Qx, frame_Qy, frame_Qz;

    std::vector<Mantid::Geometry::MDHistoDimension_sptr> dimensions_Q{{
        std::make_shared<Geometry::MDHistoDimension>("Qx", "Qx", frame_Qx, Qxminmax[0], Qxminmax[1], num_Qx_bins),
        std::make_shared<Geometry::MDHistoDimension>("Qy", "Qy", frame_Qy, Qyminmax[0], Qyminmax[1], num_Qy_bins),
        std::make_shared<Geometry::MDHistoDimension>("Qz", "Qz", frame_Qz, Qzminmax[0], Qzminmax[1], num_Qz_bins),
    }};

    m_workspace_event = std::make_shared<t_eventworkspace>(API::NoNormalization, API::NoNormalization);

    m_workspace_event->addDimension(dimensions_Q[0]);
    m_workspace_event->addDimension(dimensions_Q[1]);
    m_workspace_event->addDimension(dimensions_Q[2]);
    m_workspace_event->setCoordinateSystem(Kernel::QSample);
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
      setProperty("OutputHistogramWorkspace", std::dynamic_pointer_cast<API::IMDHistoWorkspace>(m_workspace_histo));
    }

    if (m_workspace_event) {
      setProperty("OutputEventWorkspace", std::dynamic_pointer_cast<API::IMDEventWorkspace>(m_workspace_event));
    }
  }
}

} // namespace DataHandling
} // namespace Mantid

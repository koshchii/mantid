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

  // create Q events workspace
  declareProperty(
      std::make_unique<Kernel::PropertyWithValue<bool>>("Create Event Workspace", false, Kernel::Direction::Input),
      "Create the workspace with Q events.");

  // output workspaces
  declareProperty(std::make_unique<API::WorkspaceProperty<API::IMDHistoWorkspace>>(
                      "Output Histogram Workspace", "", Kernel::Direction::Output, API::PropertyMode::Mandatory),
                  "An output histogram workspace for the detector pixel values.");
  // we do not really use events, but some of the following algorithms expect an event workspace
  declareProperty(std::make_unique<API::WorkspaceProperty<API::IMDEventWorkspace>>(
                      "Output Event Workspace", "", Kernel::Direction::Output, API::PropertyMode::Optional),
                  "An output event workspace in momentum coordinates.");
}

//----------------------------------------------------------------------------------------------
/** Execute the algorithm.
 */
void LoadILLSingleCrystal::exec() {
  bool create_event_workspace = getProperty("Create Event Workspace");
  if (getProperty("Output Event Workspace").operator std::string() == "")
    create_event_workspace = false;

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
  float wavenumber = float(2. * M_PI) / wavelength;
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
  float dist_sample_det = *get_value<float>(*m_file, "sample_distance");
  float det_angular_width = *get_value<float>(*m_file, "angular_width");
  m_file->closeGroup();

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

  auto scanned_variables_values = get_values<double>(*m_file, "data");

  m_file->openGroup("variables_names", "ILL_data_scan_vars");
  auto variables_names_entries = m_file->getEntries();

  m_log.information() << "\nData group -> Scanned variables -> Variables names:" << std::endl;
  for (const auto &entry : variables_names_entries)
    m_log.information() << entry.first << " = " << entry.second << std::endl;

  auto scanned_variables_names = get_values<std::string>(*m_file, "name");
  auto scanned_variables = get_values<std::uint8_t>(*m_file, "scanned");

  // get index of scanned variable
  std::size_t scanned_idx = 0;
  for (std::size_t idx = 0; idx < scanned_variables.size(); ++idx) {
    if (scanned_variables[idx]) {
      scanned_idx = idx;
      break;
    }
  }
  m_log.information() << "Scanned variable index: " << scanned_idx << std::endl;

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
      std::make_shared<Geometry::MDHistoDimension>("Scanned Variable", "scanned", frame_scanned, 0,
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

  m_workspace_histo->addExperimentInfo(info);
  m_workspace_histo->setTitle(title);

  // minimum and maximum Qs
  coord_t Qxminmax[2] = {std::numeric_limits<coord_t>::max(), std::numeric_limits<coord_t>::lowest()};
  coord_t Qyminmax[2] = {std::numeric_limits<coord_t>::max(), std::numeric_limits<coord_t>::lowest()};
  coord_t Qzminmax[2] = {std::numeric_limits<coord_t>::max(), std::numeric_limits<coord_t>::lowest()};

  // size of one detector image in pixels
  std::int64_t frame_size = detector_data_dims[1] * detector_data_dims[2];

  // data type for event workspace
  using t_event = DataObjects::MDEvent<3>;
  using t_eventworkspace = DataObjects::MDEventWorkspace<t_event, 3>;
  std::vector<t_event> events;
  events.reserve(detector_data.size());

  // copy the detector data to the workspace
  std::int64_t last_scan_step = -1;
  for (std::size_t idx = 0; idx < detector_data.size(); ++idx) {
    // completed all pixels of a frame?
    std::int64_t cur_scan_step = idx / frame_size;
    std::int64_t cur_frame_step = idx % frame_size;
    std::int64_t cur_y = cur_frame_step / detector_data_dims[1];
    std::int64_t cur_x = cur_frame_step % detector_data_dims[1];

    // save pixels in histogram workspace
    Mantid::signal_t intensity = detector_data[idx];
    m_workspace_histo->setSignalAt(idx, intensity);

    double scanned_variable_value = -1;
    if (last_scan_step != cur_scan_step) {
      // only update this once per frame
      scanned_variable_value = scanned_variables_values[detector_data_dims[0] * scanned_idx + cur_scan_step];
      last_scan_step = cur_scan_step;
    }

    // TODO: check if scanned_variable_value really corresponds to omega
    coord_t omega = coord_t(scanned_variable_value * M_PI / 180.);

    if (create_event_workspace) {
      coord_t pix_coord[2] = {coord_t(cur_x), coord_t(cur_y)};

      // in-plane scattering angle
      coord_t phi = pix_coord[0] / coord_t(detector_data_dims[0]) * det_angular_width;
      phi = phi / 180.f * coord_t(M_PI);

      // out-of-plane scattering angle
      float det_angular_height = 40.f; // TODO
      coord_t theta = pix_coord[1] / coord_t(detector_data_dims[1]) * det_angular_height;
      theta -= det_angular_height * 0.5;
      theta = theta / 180.f * coord_t(M_PI);

      // convert pixels to Q coordinates and insert events
      coord_t Q_coord[3] = {
          wavenumber * std::sin(theta) * std::cos(phi),
          wavenumber * std::sin(theta) * std::sin(phi),
          wavenumber * std::cos(theta),
      };

      // rotate by sample omega angle
      coord_t Q_coord_rot[3] = {
          Q_coord[0] * std::cos(omega) - Q_coord[1] * std::sin(omega),
          Q_coord[0] * std::sin(omega) + Q_coord[1] * std::cos(omega),
          Q_coord[2],
      };

      // calculate Q ranges
      Qxminmax[0] = std::min(Q_coord_rot[0], Qxminmax[0]);
      Qxminmax[1] = std::max(Q_coord_rot[0], Qxminmax[1]);
      Qyminmax[0] = std::min(Q_coord_rot[1], Qyminmax[0]);
      Qyminmax[1] = std::max(Q_coord_rot[1], Qyminmax[1]);
      Qzminmax[0] = std::min(Q_coord_rot[2], Qzminmax[0]);
      Qzminmax[1] = std::max(Q_coord_rot[2], Qzminmax[1]);

      // store all events
      events.emplace_back(t_event(intensity, std::sqrt(intensity), Q_coord_rot));
    }
  }

  m_log.information() << "Q_x range: [" << Qxminmax[0] << ", " << Qxminmax[1] << "]." << std::endl;
  m_log.information() << "Q_y range: [" << Qyminmax[0] << ", " << Qyminmax[1] << "]." << std::endl;
  m_log.information() << "Q_z range: [" << Qzminmax[0] << ", " << Qzminmax[1] << "]." << std::endl;

  // create event workspace
  if (create_event_workspace) {
    Geometry::QSample frame_Qx, frame_Qy, frame_Qz;

    std::size_t num_Q_bins = 32;
    std::vector<Mantid::Geometry::MDHistoDimension_sptr> dimensions_Q{{
        std::make_shared<Geometry::MDHistoDimension>("Qx", "Qx", frame_Qx, Qxminmax[0], Qxminmax[1], num_Q_bins),
        std::make_shared<Geometry::MDHistoDimension>("Qy", "Qy", frame_Qy, Qyminmax[0], Qyminmax[1], num_Q_bins),
        std::make_shared<Geometry::MDHistoDimension>("Qz", "Qz", frame_Qz, Qzminmax[0], Qzminmax[1], num_Q_bins),
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
    for (const t_event &event : events)
      std::dynamic_pointer_cast<t_eventworkspace>(m_workspace_event)->addEvent(event);
  }

  // set ouput workspaces
  if (m_workspace_histo) {
    setProperty("Output Histogram Workspace", std::dynamic_pointer_cast<API::IMDHistoWorkspace>(m_workspace_histo));
  }

  if (m_workspace_event) {
    setProperty("Output Event Workspace", std::dynamic_pointer_cast<API::IMDEventWorkspace>(m_workspace_event));
  }
}

} // namespace DataHandling
} // namespace Mantid

// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidAPI/Algorithm.h"
#include "MantidAPI/IFileLoader.h"
#include "MantidAPI/IMDEventWorkspace.h"
#include "MantidAPI/IMDHistoWorkspace.h"
#include "MantidDataHandling/DllConfig.h"
#include "MantidKernel/NexusDescriptor.h"
#include "MantidKernel/System.h"

#include <nexus/NeXusFile.hpp>

namespace Mantid {
namespace DataHandling {

/** LoadILLSingleCrystal
 */
class MANTID_DATAHANDLING_DLL LoadILLSingleCrystal
    : public /*API::Algorithm*/ API::IFileLoader<Kernel::NexusDescriptor> {
public:
  // number data types
  using t_real = coord_t;
  using t_data = double;

public:
  virtual const std::string name() const override;
  virtual int version() const override;
  virtual const std::string category() const override;
  virtual const std::string summary() const override;
  virtual int confidence(Kernel::NexusDescriptor &descr) const override;

protected:
  bool LoadScannedVariables();
  bool LoadInstrumentGroup();
  bool LoadSampleGroup();
  bool LoadUserGroup();
  bool LoadMonitorGroup();
  bool LoadDataScanGroup();

private:
  virtual void init() override;
  virtual void exec() override;

  std::string m_filename;
  std::unique_ptr<NeXus::File> m_file;

  API::IMDHistoWorkspace_sptr m_workspace_histo;
  API::IMDEventWorkspace_sptr m_workspace_event;

  // --------------------------------------------------------------------------
  // instrument group
  // --------------------------------------------------------------------------
  std::string m_instr_name;

  t_real m_cell[3];
  t_real m_cell_angles[3];
  std::vector<t_real> m_UB;

  t_real m_chi_deg, m_omega_deg, m_phi_deg;
  int m_chi_axis, m_omega_axis, m_phi_axis;

  t_real m_dist_sample_det;
  t_real m_det_angular_width, m_det_angular_height;
  // --------------------------------------------------------------------------

  // --------------------------------------------------------------------------
  // monitor group
  // --------------------------------------------------------------------------
  t_real m_monitor_sum;
  // --------------------------------------------------------------------------

  // --------------------------------------------------------------------------
  // data scan group
  // --------------------------------------------------------------------------
  std::vector<std::int64_t> m_detector_data_dims;
  std::vector<std::uint32_t> m_detector_data;

  std::vector<std::string> m_scanned_variables_names;
  std::vector<t_data> m_scanned_variables_values;

  std::size_t m_first_scanned_idx;
  std::size_t m_chi_idx, m_omega_idx, m_phi_idx;
  bool m_first_scanned_idx_set;
  bool m_chi_idx_set, m_omega_idx_set, m_phi_idx_set;
  // --------------------------------------------------------------------------
};

} // namespace DataHandling
} // namespace Mantid

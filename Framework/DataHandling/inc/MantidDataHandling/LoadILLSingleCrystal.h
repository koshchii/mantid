// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidAPI/Algorithm.h"
#include "MantidAPI/IMDEventWorkspace.h"
#include "MantidAPI/IMDHistoWorkspace.h"
#include "MantidDataHandling/DllConfig.h"
#include "MantidKernel/System.h"
#include <nexus/NeXusFile.hpp>

namespace Mantid {
namespace DataHandling {

/** LoadILLSingleCrystal : TODO: DESCRIPTION
 */
class MANTID_DATAHANDLING_DLL LoadILLSingleCrystal : public API::Algorithm {
public:
  const std::string name() const override;
  int version() const override;
  const std::string category() const override;
  const std::string summary() const override;

private:
  void init() override;
  void exec() override;

  std::string m_filename;
  std::unique_ptr<NeXus::File> m_file;
  API::IMDHistoWorkspace_sptr m_workspace_histo;
  API::IMDEventWorkspace_sptr m_workspace_event;
};

} // namespace DataHandling
} // namespace Mantid

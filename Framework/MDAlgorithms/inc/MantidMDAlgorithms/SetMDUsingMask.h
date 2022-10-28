// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2011 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidAPI/Algorithm.h"
#include "MantidMDAlgorithms/DllConfig.h"

namespace Mantid {
namespace MDAlgorithms {

/** Algorithm to set a MDHistoWorkspace in points determined by a mask boolean
  MDHistoWorkspace.

  @date 2011-11-10
*/
class MANTID_MDALGORITHMS_DLL SetMDUsingMask final : public API::Algorithm {
public:
  const std::string name() const override;
  /// Summary of algorithms purpose
  const std::string summary() const override {
    return "Algorithm to set a MDHistoWorkspace in points determined by a mask "
           "boolean MDHistoWorkspace.";
  }

  int version() const override;
  const std::string category() const override;

private:
  void init() override;
  void exec() override;
};

} // namespace MDAlgorithms
} // namespace Mantid

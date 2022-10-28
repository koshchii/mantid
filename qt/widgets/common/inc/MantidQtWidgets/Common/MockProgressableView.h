// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidKernel/WarningSuppressions.h"
#include "MantidQtWidgets/Common/ProgressableView.h"
#include <gmock/gmock.h>

using namespace MantidQt::MantidWidgets;

class MockProgressableView : public ProgressableView {
public:
  GNU_DIAG_OFF_SUGGEST_OVERRIDE
  MOCK_METHOD1(setProgress, void(int));
  MOCK_METHOD0(clearProgress, void());
  MOCK_METHOD2(setProgressRange, void(int, int));
  MOCK_CONST_METHOD0(isPercentageIndicator, bool());
  MOCK_METHOD0(setAsPercentageIndicator, void());
  MOCK_METHOD0(setAsEndlessIndicator, void());
  GNU_DIAG_ON_SUGGEST_OVERRIDE
};

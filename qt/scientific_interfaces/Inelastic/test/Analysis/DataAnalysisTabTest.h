// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2019 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "Analysis/DataAnalysisTab.h"
#include "MantidAPI/FunctionFactory.h"
#include <cxxtest/TestSuite.h>
#include <gmock/gmock.h>

using namespace MantidQt::CustomInterfaces::IDA;

class DataAnalysisTabTest : public CxxTest::TestSuite {
public:
  static DataAnalysisTabTest *createSuite() { return new DataAnalysisTabTest(); }

  static void destroySuite(DataAnalysisTabTest *suite) { delete suite; }

  void test_that_single_function_correctly_identified() {
    std::string functionName = "ExpDecay";
    auto fitFunction = Mantid::API::FunctionFactory::Instance().createFunction(functionName);
    auto occurances = DataAnalysisTab::getNumberOfSpecificFunctionContained(functionName, fitFunction.get());
    TS_ASSERT_EQUALS(occurances, 1);
  }

  void test_that_single_layer_composite_function_handled_correctly() {
    std::string functionName = "name=ExpDecay;name=StretchExp";
    auto fitFunction = Mantid::API::FunctionFactory::Instance().createInitialized(functionName);
    auto occurances = DataAnalysisTab::getNumberOfSpecificFunctionContained("ExpDecay", fitFunction.get());
    auto stretchOccurances = DataAnalysisTab::getNumberOfSpecificFunctionContained("StretchExp", fitFunction.get());
    TS_ASSERT_EQUALS(occurances, 1);
    TS_ASSERT_EQUALS(stretchOccurances, 1);
  }

  void test_that_no_matched_name_is_correct() {
    std::string functionName = "name=ExpDecay;name=StretchExp";
    auto fitFunction = Mantid::API::FunctionFactory::Instance().createInitialized(functionName);
    auto occurances = DataAnalysisTab::getNumberOfSpecificFunctionContained("NotHere", fitFunction.get());
    TS_ASSERT_EQUALS(occurances, 0);
  }

  void test_that_multi_layer_composite_function_handled_correctly() {
    std::string functionName = "name=ExpDecay;name=ExpDecay;(composite="
                               "ProductFunction,NumDeriv=false;name=ExpDecay;"
                               "name=ExpDecay)";
    auto fitFunction = Mantid::API::FunctionFactory::Instance().createInitialized(functionName);
    auto occurances = DataAnalysisTab::getNumberOfSpecificFunctionContained("ExpDecay", fitFunction.get());
    TS_ASSERT_EQUALS(occurances, 4);
  }
};

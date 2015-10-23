#ifndef MANTID_MDALGORITHMS_ACCUMULATEMDTEST_H_
#define MANTID_MDALGORITHMS_ACCUMULATEMDTEST_H_

#include <cxxtest/TestSuite.h>
#include "MantidMDAlgorithms/AccumulateMD.h"
#include "MantidDataObjects/MDHistoWorkspace.h"
#include "MantidTestHelpers/MDEventsTestHelper.h"
#include "MantidTestHelpers/WorkspaceCreationHelper.h"
#include "MantidAPI/IMDHistoWorkspace.h"
#include "MantidKernel/ConfigService.h"
#include <Poco/Path.h>
#include <Poco/File.h>

using Mantid::MDAlgorithms::AccumulateMD;
using namespace Mantid::API;
using namespace Mantid::DataObjects;

class AccumulateMDTest : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static AccumulateMDTest *createSuite() { return new AccumulateMDTest(); }

  static void destroySuite(AccumulateMDTest *suite) { delete suite; }

  void test_Init() {
    AccumulateMD alg;
    TS_ASSERT_THROWS_NOTHING(alg.initialize());
  }

  void test_filter_to_existing_sources_nothing_exist() {
    // Create vector of data_sources to filter
    std::vector<std::string> data_sources;

    // Add absolute path to a file which doesn't exist
    Poco::Path filepath =
        Poco::Path(Mantid::Kernel::ConfigService::Instance().getTempDir(),
                   "ACCUMULATEMDTEST_NONEXISTENTFILE");
    data_sources.push_back(filepath.toString());

    Mantid::Kernel::Logger logger("AccumulateMDTest");
    data_sources =
        Mantid::MDAlgorithms::filterToExistingSources(data_sources, logger);

    TS_ASSERT(data_sources.empty());
  }

  void test_filter_to_existing_sources_workspace_exist() {
    // Create vector of data_sources to filter
    std::vector<std::string> data_sources;

    // Create a cheap workspace
    std::string ws_name = "ACCUMULATEMDTEST_EXISTENTWORKSPACE";
    auto bkgWS = WorkspaceCreationHelper::Create1DWorkspaceRand(1);
    // add to ADS
    AnalysisDataService::Instance().add(ws_name, bkgWS);

    data_sources.push_back(ws_name);

    Mantid::Kernel::Logger logger("AccumulateMDTest");
    data_sources =
        Mantid::MDAlgorithms::filterToExistingSources(data_sources, logger);

    TS_ASSERT(!data_sources.empty());

    // Remove workspace from the data service.
    AnalysisDataService::Instance().remove(ws_name);
  }

  void test_filter_to_existing_sources_file_exist() {
    // Create vector of data_sources to filter
    std::vector<std::string> data_sources;

    // Create a temporary file to find
    Poco::Path filepath =
        Poco::Path(Mantid::Kernel::ConfigService::Instance().getTempDir(),
                   "ACCUMULATEMDTEST_EXISTENTFILE");
    Poco::File existent_file(filepath);
    existent_file.createFile();
    data_sources.push_back(filepath.toString());

    Mantid::Kernel::Logger logger("AccumulateMDTest");
    data_sources =
        Mantid::MDAlgorithms::filterToExistingSources(data_sources, logger);

    TS_ASSERT(!data_sources.empty());

    // Remove the temporary file
    existent_file.remove();
  }

  void test_filter_to_new_none_new() {
    std::vector<std::string> input_data;
    input_data.push_back("test1");
    input_data.push_back("test2");
    input_data.push_back("test3");
    std::vector<std::string> current_data = input_data;

    std::vector<std::string> result =
        Mantid::MDAlgorithms::filterToNew(input_data, current_data);

    // Two input vectors were identical, so we should get an empty vector back
    TS_ASSERT(result.empty());
  }

  void test_filter_to_new() {
    std::vector<std::string> input_data;
    input_data.push_back("test1");
    input_data.push_back("test2");
    input_data.push_back("test3");
    input_data.push_back("test4");
    input_data.push_back("test5");

    std::vector<std::string> current_data;
    current_data.push_back("test1");
    current_data.push_back("test3");
    current_data.push_back("test4");

    std::vector<std::string> result =
        Mantid::MDAlgorithms::filterToNew(input_data, current_data);

    // test2 and test5 is new data (it is in input_data but not current_data)
    // and so should be returned in the vector
    TS_ASSERT_EQUALS(result[0], "test2");
    TS_ASSERT_EQUALS(result[1], "test5");
  }
};

#endif /* MANTID_MDALGORITHMS_ACCUMULATEMDTEST_H_ */

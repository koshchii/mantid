// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2022 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include <cmath>
#include <cxxtest/TestSuite.h>

#include "MantidAPI/AnalysisDataService.h"
#include "MantidDataHandling/LoadILLSALSA.h"
#include "MantidGeometry/IComponent.h"
#include "MantidGeometry/Instrument.h"
#include "MantidKernel/V3D.h"

class LoadILLSALSATest : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static LoadILLSALSATest *createSuite() { return new LoadILLSALSATest(); }
  static void destroySuite(LoadILLSALSATest *suite) { delete suite; }

  LoadILLSALSATest() {
    Mantid::Kernel::ConfigService::Instance().appendDataSearchSubDir("ILL/SALSA/");
    Mantid::Kernel::ConfigService::Instance().setFacility("ILL");
  }

  void tearDown() override { Mantid::API::AnalysisDataService::Instance().clear(); }

  void test_name() {
    Mantid::DataHandling::LoadILLSALSA alg;
    TS_ASSERT_EQUALS(alg.name(), "LoadILLSALSA");
  }

  void test_version() {
    Mantid::DataHandling::LoadILLSALSA alg;
    TS_ASSERT_EQUALS(alg.version(), 1);
  }

  void test_init() {
    Mantid::DataHandling::LoadILLSALSA alg;
    TS_ASSERT_THROWS_NOTHING(alg.initialize())
    TS_ASSERT(alg.isInitialized())
  }

  void test_loadV1() {
    Mantid::DataHandling::LoadILLSALSA alg;
    alg.setChild(true);
    alg.initialize();
    TS_ASSERT_THROWS_NOTHING(alg.setPropertyValue("Filename", "046508.nxs"))
    TS_ASSERT_THROWS_NOTHING(alg.setPropertyValue("DetectorDistance", "3.5"))
    TS_ASSERT_THROWS_NOTHING(alg.setPropertyValue("OutputWorkspace", "__unused_for_child"))
    TS_ASSERT_THROWS_NOTHING(alg.execute())
    TS_ASSERT(alg.isExecuted())
    Mantid::API::MatrixWorkspace_const_sptr outputWS = alg.getProperty("OutputWorkspace");
    TS_ASSERT(outputWS)
    TS_ASSERT_EQUALS(outputWS->getNumberHistograms(), 256 * 256 + 1)
    TS_ASSERT_EQUALS(outputWS->blocksize(), 1)
    auto component = outputWS->getInstrument()->getComponentByName("detector");
    auto pos = component->getPos();
    TS_ASSERT_DELTA(sqrt(pow(pos.Z(), 2) + pow(pos.Y(), 2) + pow(pos.X(), 2)), 3.5, 1e-5);
    TS_ASSERT_DELTA(acos(pos.Z() / 3.5) * 180.0 / M_PI, 100.502, 1e-5)
  }

  void test_loadV2() {
    Mantid::DataHandling::LoadILLSALSA alg;
    alg.setChild(true);
    alg.initialize();
    TS_ASSERT_THROWS_NOTHING(alg.setPropertyValue("Filename", "046430.nxs"))
    TS_ASSERT_THROWS_NOTHING(alg.setPropertyValue("OutputWorkspace", "__unused_for_child"))
    TS_ASSERT_THROWS_NOTHING(alg.execute())
    TS_ASSERT(alg.isExecuted())
    Mantid::API::MatrixWorkspace_const_sptr outputWS = alg.getProperty("OutputWorkspace");
    TS_ASSERT(outputWS)
    TS_ASSERT_EQUALS(outputWS->getNumberHistograms(), 256 * 256 + 1)
    TS_ASSERT_EQUALS(outputWS->blocksize(), 10)
    auto component = outputWS->getInstrument()->getComponentByName("detector");
    auto pos = component->getPos();
    TS_ASSERT_DELTA(sqrt(pow(pos.Z(), 2) + pow(pos.Y(), 2) + pow(pos.X(), 2)), 1.0, 1e-5);
    TS_ASSERT_DELTA(acos(pos.Z() / 1.0) * 180.0 / M_PI, 85.603, 1e-5)
  }
};

class LoadILLSALSATestPerformance : public CxxTest::TestSuite {
public:
  void setUp() override {
    m_alg.initialize();
    m_alg.setChild(true);
    m_alg.setPropertyValue("Filename", "ILL/SALSA/046430.nxs");
    m_alg.setPropertyValue("OutputWorkspace", "__unused_for_child");
  }
  void testLoadILLSALSAPerformance() {
    for (int i = 0; i < 10; ++i) {
      TS_ASSERT_THROWS_NOTHING(m_alg.execute());
    }
  }
  void tearDown() override { Mantid::API::AnalysisDataService::Instance().clear(); }

private:
  Mantid::DataHandling::LoadILLSALSA m_alg;
};

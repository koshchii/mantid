#ifndef IMD_DIMENSIONFACTORYTEST_H_
#define IMD_DIMENSIONFACTORYTEST_H_

#include <cxxtest/TestSuite.h>
#include "MantidGeometry/MDGeometry/IMDDimensionFactory.h"
#include <Poco/DOM/DOMParser.h>
#include <Poco/DOM/Document.h>
#include <Poco/AutoPtr.h>

using namespace Mantid::Geometry;

class IMDDimensionFactoryTest: public CxxTest::TestSuite
{
private:

  std::string constructDimensionWithUnitsXMLString()
  {
    std::string xmlToParse = std::string("<Dimension ID=\"qz\">") + "<Name>Qz</Name>"
        + "<Units>Cubits</Units>"
        + "<UpperBounds>3</UpperBounds>" + "<LowerBounds>-3</LowerBounds>"
        + "<NumberOfBins>8</NumberOfBins>"
        + "<ReciprocalDimensionMapping>q3</ReciprocalDimensionMapping>" + "</Dimension>";

    return xmlToParse;
  }

  std::string constructDimensionWithoutUnitsXMLString()
  {
    std::string xmlToParse = std::string("<Dimension ID=\"qz\">") + "<Name>Qz</Name>"
        + "<UpperBounds>3</UpperBounds>" + "<LowerBounds>-3</LowerBounds>"
        + "<NumberOfBins>8</NumberOfBins>"
        + "<ReciprocalDimensionMapping>q3</ReciprocalDimensionMapping>" + "</Dimension>";

    return xmlToParse;
  }

  std::string constructNonReciprocalDimensionXMLString()
  {
    return std::string("<Dimension ID=\"en\">") + "<Name>Energy</Name>"
        + "<UpperBounds>150</UpperBounds>" + "<LowerBounds>0</LowerBounds>"
        + "<NumberOfBins>4</NumberOfBins>" + "</Dimension>";
  }

  Poco::AutoPtr<Poco::XML::Document> constructNonReciprocalDimensionXML()
  {
    std::string xmlToParse = constructNonReciprocalDimensionXMLString();
    Poco::XML::DOMParser pParser;
    return pParser.parseString(xmlToParse);
  }

public:

  void testCorrectGeneration()
  {
    IMDDimension_const_sptr dimension = IMDDimensionFactory::create(constructDimensionWithUnitsXMLString());
    TS_ASSERT_EQUALS("Cubits", dimension->getUnits());
    TS_ASSERT_EQUALS("Qz", dimension->getName());
    TS_ASSERT_EQUALS("qz", dimension->getDimensionId());
    TS_ASSERT_EQUALS(-3, dimension->getMinimum());
    TS_ASSERT_EQUALS(3, dimension->getMaximum());
    TS_ASSERT_EQUALS(8, dimension->getNBins());
  }

  void testCorrectGenerationWithoutUnits()
  {
    IMDDimension_const_sptr dimension = IMDDimensionFactory::create(constructDimensionWithoutUnits());
    TS_ASSERT_EQUALS("None", dimension->getUnits());
    TS_ASSERT_EQUALS("Qz", dimension->getName());
    TS_ASSERT_EQUALS("qz", dimension->getDimensionId());
    TS_ASSERT_EQUALS(-3, dimension->getMinimum());
    TS_ASSERT_EQUALS(3, dimension->getMaximum());
    TS_ASSERT_EQUALS(8, dimension->getNBins());
  }

  void testStaticCreation()
  {
    std::string xmlToParse = constructNonReciprocalDimensionXMLString();
    IMDDimension_const_sptr viaString = IMDDimensionFactory::create(xmlToParse);
    auto document = constructNonReciprocalDimensionXML();
    IMDDimension_const_sptr viaXML = IMDDimensionFactory::create(document->documentElement());

    //Constructed either way, the products should be equivalent
    TSM_ASSERT_EQUALS("Created through either route, the products should be equal", viaString->getDimensionId(), viaXML->getDimensionId());

  }
};
#endif

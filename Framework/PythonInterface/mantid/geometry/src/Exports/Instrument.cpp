// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidGeometry/Instrument.h"
#include "MantidGeometry/Instrument/ReferenceFrame.h"
#include "MantidKernel/WarningSuppressions.h"
#include "MantidPythonInterface/core/GetPointer.h"
#include "MantidPythonInterface/core/Policies/RemoveConst.h"

#include <boost/python/class.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/register_ptr_to_python.hpp>

using namespace Mantid::Geometry;
using Mantid::detid_t;
using namespace boost::python;
using Mantid::PythonInterface::Policies::RemoveConstSharedPtr;

GET_POINTER_SPECIALIZATION(Instrument)

namespace {

// Ignore -Wconversion warnings coming from boost::python
// Seen with GCC 7.1.1 and Boost 1.63.0
GNU_DIAG_OFF("conversion")

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Instrument_getNumberDetectors, Instrument::getNumberDetectors, 0, 1)

GNU_DIAG_ON("conversion")

} // namespace

void export_Instrument() {
  register_ptr_to_python<std::shared_ptr<Instrument>>();

  boost::python::enum_<Instrument::ContainsState>("ContainsState")
      .value("Full", Instrument::ContainsState::Full)
      .value("Partial", Instrument::ContainsState::Partial)
      .value("None", Instrument::ContainsState::None);

  class_<Instrument, bases<CompAssembly>, boost::noncopyable>("Instrument", no_init)
      .def("getSample", &Instrument::getSample, arg("self"), return_value_policy<RemoveConstSharedPtr>(),
           "Return the :class:`~mantid.geometry.Component` object that "
           "represents the sample")

      .def("getSource", &Instrument::getSource, arg("self"), return_value_policy<RemoveConstSharedPtr>(),
           "Return the :class:`~mantid.geometry.Component` object that "
           "represents the source")

      .def("getComponentByName",
           (std::shared_ptr<const IComponent>(Instrument::*)(const std::string &, int) const) &
               Instrument::getComponentByName,
           (arg("self"), arg("cname"), arg("nlevels") = 0), "Returns the named :class:`~mantid.geometry.Component`")

      .def("getDetector",
           (std::shared_ptr<const IDetector>(Instrument::*)(const detid_t &) const) & Instrument::getDetector,
           (arg("self"), arg("detector_id")), "Returns the :class:`~mantid.geometry.Detector` with the given ID")

      .def("getNumberDetectors", &Instrument::getNumberDetectors,
           Instrument_getNumberDetectors((arg("self"), arg("skipMonitors") = false)))

      .def("getReferenceFrame",
           (std::shared_ptr<const ReferenceFrame>(Instrument::*)()) & Instrument::getReferenceFrame, arg("self"),
           return_value_policy<RemoveConstSharedPtr>(),
           "Returns the :class:`~mantid.geometry.ReferenceFrame` attached that "
           "defines the instrument "
           "axes")

      .def("getValidFromDate", &Instrument::getValidFromDate, arg("self"),
           "Return the valid from :class:`~mantid.kernel.DateAndTime` of the "
           "instrument")

      .def("getValidToDate", &Instrument::getValidToDate, arg("self"),
           "Return the valid to :class:`~mantid.kernel.DateAndTime` of the "
           "instrument")

      .def("getBaseInstrument", &Instrument::baseInstrument, arg("self"), return_value_policy<RemoveConstSharedPtr>(),
           "Return reference to the base instrument")

      .def("containsRectDetectors", &Instrument::containsRectDetectors, arg("self"),
           "Return Full if all detectors are rect., Partial if some, None if none");
}

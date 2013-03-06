#ifndef MANTID_SLICEVIEWER_PEAKTRANSFORMQSAMPLE_H_
#define MANTID_SLICEVIEWER_PEAKTRANSFORMQSAMPLE_H_

#include "MantidQtSliceViewer/PeakTransform.h"
#include "MantidQtSliceViewer/ConcretePeakTransformFactory.h"

namespace MantidQt
{
  namespace SliceViewer
  {
    /**
    @class PeakTransformQSample
    Used to remap coordinates into a form consistent with an axis reordering.
    */
    class DLLExport PeakTransformQSample : public PeakTransform
    {
    public:
      /// Transform name.
      static std::string name(){ return "QLab"; }
      /// Constructor
      PeakTransformQSample();
      /// Constructor
      PeakTransformQSample(const std::string& xPlotLabel, const std::string& yPlotLabel);
      /// Destructor
      virtual ~PeakTransformQSample();
      /// Copy constructor
      PeakTransformQSample(const PeakTransformQSample& other);
      /// Assigment
      PeakTransformQSample & operator=(const PeakTransformQSample & other);
      /// Virtual constructor
      PeakTransform_sptr clone() const;
      /// Transform peak.
      Mantid::Kernel::V3D transformPeak(const Mantid::API::IPeak& peak) const;
      /// Getter for the transform name.
      virtual std::string getFriendlyName() const
      {
        return name();
      }
      /// Getter for the special coordinate representation of this transform type.
      Mantid::API::SpecialCoordinateSystem getCoordinateSystem() const;
    };

    /// Typedef a factory for type of PeaksTransform.
    typedef ConcretePeakTransformFactory<PeakTransformQSample> PeakTransformQSampleFactory;

  }
}

#endif /* MANTID_SLICEVIEWER_PEAKTRANSFORMQSAMPLE_H_ */

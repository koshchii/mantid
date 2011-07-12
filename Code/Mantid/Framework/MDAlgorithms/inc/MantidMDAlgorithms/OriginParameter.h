#ifndef ORIGINPARAMETER_H_
#define ORIGINPARAMETER_H_

//----------------------------------------------------------------------
// Includes
//----------------------------------------------------------------------
#include <vector>
#include "MantidKernel/System.h"
#include "MantidAPI/ImplicitFunctionParameter.h"

namespace Mantid
{
    namespace MDAlgorithms
    {
        /**

        OriginParameter. Wraps a vector expressing origin location.

        @author Owen Arnold, Tessella plc
        @date 01/10/2010

        Copyright &copy; 2010 ISIS Rutherford Appleton Laboratory & NScD Oak Ridge National Laboratory

        This file is part of Mantid.

        Mantid is free software; you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation; either version 3 of the License, or
        (at your option) any later version.

        Mantid is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more details.

        You should have received a copy of the GNU General Public License
        along with this program.  If not, see <http://www.gnu.org/licenses/>.

        File change history is stored at: <https://svn.mantidproject.org/mantid/trunk/Code/Mantid>
        Code Documentation is available at: <http://doxygen.mantidproject.org>
        */

        class DLLExport OriginParameter :public Mantid::API::ImplicitFunctionParameter
        {
        private:

           std::vector<double> m_origin;

        public:

            OriginParameter(double o1, double o2, double o3);
            
            OriginParameter();
            
            OriginParameter(const OriginParameter & other);
            
            OriginParameter& operator=(const OriginParameter& other);

            bool operator==(const OriginParameter &other) const;

            bool operator!=(const OriginParameter &other) const;

            bool isValid() const;

            std::string getName() const;

            void asVector(std::vector<double>& origin) const;

            OriginParameter* clone() const;

            double getX() const;

            double getY() const;

            double getZ() const;

            std::string toXMLString() const;

            ~OriginParameter();
            
            /*
            Treat the origin parameter as an array of scalar values.
            @param index : the index to fetch (index is one of x=0, y, z)
            @return the origin value at the specified index.
            */
            double& operator[] (int index)
            {
              return m_origin[index];
            }

            /*
            Getter for the parameter name associated with this type.
            @return the parameter name.
            */
            static std::string parameterName()
            {
                return "OriginParameter";
            }
        };
    }
}

#endif

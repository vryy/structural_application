/*
==============================================================================
KratosR1StructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: nagel $
//   Date:                $Date: 2007-10-16 12:48:37 $
//   Revision:            $Revision: 1.1 $
//
//
#if !defined(KRATOS_TIP_CONDITION_H_INCLUDED )
#define  KRATOS_TIP_CONDITION_H_INCLUDED


// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "utilities/math_utils.h"

namespace Kratos
{
/**
 * Tip condition to tie the pile to the ground
 */
class TipCondition : public Condition
{
    public:
        // Counted pointer of TipCondition
        KRATOS_CLASS_POINTER_DEFINITION(TipCondition);

        typedef Element::GeometryType GeometryType;
        typedef GeometryType::PointType NodeType;
        typedef NodeType::PointType PointType;

        /**
         * Default constructor.
         */
        TipCondition();
        TipCondition( IndexType NewId, GeometryType::Pointer pGeometry);

        TipCondition( IndexType NewId, GeometryType::Pointer pGeometry,
                       PropertiesType::Pointer pProperties
                     );

        TipCondition( IndexType NewId, GeometryType::Pointer pGeometry,
                      PropertiesType::Pointer pProperties,
                      Element::Pointer tip_soilElement,
                      Element::Pointer TipElement,
                      PointType& rTipSoilLocalPoint, PointType& rTipLocalPoint);

        /**
         * Destructor.
         */
        virtual ~TipCondition();

        /**
         * Operations.
         */

        Condition::Pointer Create( IndexType NewId,
                                   NodesArrayType const& ThisNodes,
                                   PropertiesType::Pointer pProperties) const;

        /**
         * Calculates the local system contributions for this contact element
         */
        void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                   VectorType& rRightHandSideVector,
                                   const ProcessInfo& rCurrentProcessInfo);

        void CalculateRightHandSide( VectorType& rRightHandSideVector,
                                     const ProcessInfo& rCurrentProcessInfo);
//          void DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo);

        void EquationIdVector( EquationIdVectorType& rResult,
                               const ProcessInfo& rCurrentProcessInfo) const;

        void GetDofList( DofsVectorType& ConditionalDofList,
                         const ProcessInfo& CurrentProcessInfo) const;

        void Initialize(const ProcessInfo& CurrentProcessInfo);
        /**
         * Turn back information as a string.
         * (DEACTIVATED)
         */
        //std::string Info();

        /**
         * Print information about this object.
         * (DEACTIVATED)
         */
        //virtual void PrintInfo(std::ostream& rOStream) const;

        /**
         * Print object's data.
         * (DEACTIVATED)
         */
        //virtual void PrintData(std::ostream& rOStream) const;

    protected:


    private:

        friend class Serializer;

        virtual void save ( Serializer& rSerializer ) const
        {
        rSerializer.save ( "name", "TipCondition" );
        KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, Condition )
        }

        virtual void load ( Serializer& rSerializer )
        {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, Condition )
        }

        void CalculateAll( MatrixType& rLeftHandSideMatrix,
                           VectorType& rRightHandSideVector,
                           const ProcessInfo& rCurrentProcessInfo,
                           bool CalculateStiffnessMatrixFlag,
                           bool CalculateResidualVectorFlag);


        PointType mTipLocalPoint;
        PointType mTipSoilLocalPoint;
        Element::Pointer mpTipElement;
        Element::Pointer mpTipSoilElement;
        PointType mTipGlobalPoint;
        PointType mTipSoilGlobalPoint;
//            Vector mTPileGlobalVector;
///     int iSoilIntegrationPointIndex;

    }; // Class TipCondition
}  // namespace Kratos.

#endif // KRATOS_TIP_CONDITION_H_INCLUDED defined


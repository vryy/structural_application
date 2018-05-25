/*
==============================================================================
KratosStructuralApplication
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
/* *********************************************************
*
*   Last Modified by:    $Author: Jelena$
*   Date:                $Date: 2014 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

#if !defined(KRATOS_POINT_POINT_LAGRANGE_CONDITION_H_INCLUDED )
#define  KRATOS_POINT_POINT_LAGRANGE_CONDITION_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/serializer.h"
//#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "custom_conditions/master_contact_point_2d.h"
#include "custom_conditions/slave_contact_point_2d.h"

namespace Kratos
{
/**
 * Joint condition to link 2 nodes by means of a 6x6 arbitrary stiffness matrix
 */
class PointPointLagrangeCondition : public Condition
{

public:
    typedef Condition BaseType;
    typedef BaseType::EquationIdVectorType EquationIdVectorType;
    typedef BaseType::MatrixType LHS_ContributionType;

    typedef Condition::GeometryType::Pointer PointerGeometryType;

    // Counted pointer of PointPointContactLink
    KRATOS_CLASS_POINTER_DEFINITION(PointPointLagrangeCondition);

    /**
     * Default constructor.
     */
    PointPointLagrangeCondition( IndexType NewId, GeometryType::Pointer pGeometry);

    PointPointLagrangeCondition( IndexType NewId, GeometryType::Pointer pGeometry,
                            PropertiesType::Pointer pProperties);

//     PointPointLagrangeCondition( IndexType NewId, NodesArrayType const& ThisNodes);

    /**
     * Custom constructor.
     */
    PointPointLagrangeCondition( IndexType NewId, Node<3>::Pointer const& node1, Node<3>::Pointer const& node2, PropertiesType::Pointer pProperties );
    
    PointPointLagrangeCondition( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties );


    /**
     * Destructor.
     */
    virtual ~PointPointLagrangeCondition();

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
                               ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSide( VectorType& rRightHandSideVector,
                                 ProcessInfo& rCurrentProcessInfo);

    void EquationIdVector( EquationIdVectorType& rResult,
                           ProcessInfo& rCurrentProcessInfo);

    void GetDofList( DofsVectorType& ConditionalDofList,
                     ProcessInfo& CurrentProcessInfo);

//     void SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo );

    
private:
    void CalculateAll( MatrixType& rLeftHandSideMatrix,
                       VectorType& rRightHandSideVector,
                       ProcessInfo& rCurrentProcessInfo,
                       bool CalculateStiffnessMatrixFlag,
                       bool CalculateResidualVectorFlag);

//     MatrixType mStiffnessMatrix;

    friend class Serializer;

    // A private default constructor necessary for serialization
    PointPointLagrangeCondition() {};

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
    }

}; // Class PointPointLagrangeCondition
}  // namespace Kratos.

#endif // KRATOS_POINT_POINT_LAGRANGE_CONDITION_H_INCLUDED  defined 

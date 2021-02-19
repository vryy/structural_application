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
*   Last Modified by:    $Author: hbui$
*   Date:                $Date: 6 Dec 2017 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

#if !defined(KRATOS_FACE_FACE_JOINT_CONDITION_H_INCLUDED )
#define  KRATOS_FACE_FACE_JOINT_CONDITION_H_INCLUDED



// System includes


// External includes


// Project includes
#include "includes/condition.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

namespace Kratos
{
/**
 * Joint condition between two surfaces/lines using Penalty method.
 * This condition requires that the two surfaces are matching geometrically and parametrically. If the surface does not match either geometrically or parametrically then the joint model using mortar shall be used.
 * The integration of the virtual works resulted by joint will be performed on the first surface.
 * This condition cannot be used in mdpa because it needs two surfaces from two conditions.
 * In post-processing, this condition will be represented as a line between two first nodes of the surfaces.
 */
class FaceFaceJointCondition : public Condition
{

public:
    typedef Condition BaseType;
    typedef BaseType::EquationIdVectorType EquationIdVectorType;
    typedef BaseType::GeometryType GeometryType;

    // Counted pointer of PointPointContactLink
    KRATOS_CLASS_POINTER_DEFINITION(FaceFaceJointCondition);

    /**
     * Custom constructor.
     */
    FaceFaceJointCondition( IndexType NewId, GeometryType::Pointer geom1, GeometryType::Pointer geom2, PropertiesType::Pointer pProperties );

    FaceFaceJointCondition( IndexType NewId, BaseType::Pointer cond1, BaseType::Pointer cond2, PropertiesType::Pointer pProperties );


    /**
     * Destructor.
     */
    virtual ~FaceFaceJointCondition();

    /**
     * Operations.
     */
    Condition::Pointer Create( IndexType NewId,
                               GeometryType::Pointer geom1, GeometryType::Pointer geom2,
                               PropertiesType::Pointer pProperties) const;

    Condition::Pointer Create( IndexType NewId,
                               BaseType::Pointer cond1, BaseType::Pointer cond2,
                               PropertiesType::Pointer pProperties) const;

    void Initialize(const ProcessInfo& rCurrentProcessInfo);

    void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                               VectorType& rRightHandSideVector,
                               const ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSide( VectorType& rRightHandSideVector,
                                 const ProcessInfo& rCurrentProcessInfo);

    void EquationIdVector( EquationIdVectorType& rResult,
                           const ProcessInfo& rCurrentProcessInfo) const;

    void GetDofList( DofsVectorType& ConditionalDofList,
                     const ProcessInfo& CurrentProcessInfo) const;


private:

    IntegrationMethod mThisIntegrationMethod;
    GeometryType::Pointer mpGeom1;
    GeometryType::Pointer mpGeom2;

    void CalculateAll( MatrixType& rLeftHandSideMatrix,
                       VectorType& rRightHandSideVector,
                       const ProcessInfo& rCurrentProcessInfo,
                       bool CalculateStiffnessMatrixFlag,
                       bool CalculateResidualVectorFlag);

    friend class Serializer;

    // A private default constructor necessary for serialization
    FaceFaceJointCondition() {};

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
    }

}; // Class FaceFaceJointCondition
}  // namespace Kratos.

#endif // KRATOS_FACE_FACE_JOINT_CONDITION_H_INCLUDED  defined

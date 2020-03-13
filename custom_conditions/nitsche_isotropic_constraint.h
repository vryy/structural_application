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
*   Date:                $Date: 9 Feb 2018 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

#if !defined(KRATOS_NITSCHE_ISOTROPIC_CONSTRAINT_CONDITION_H_INCLUDED )
#define  KRATOS_NITSCHE_ISOTROPIC_CONSTRAINT_CONDITION_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/serializer.h"
//#include "includes/ublas_interface.h"
#include "includes/variables.h"

namespace Kratos
{
/**
 * Weakly enforced boundary condition using Nitsche method for linear elasticity.
 * REF: TN7
 */
class NitscheIsotropicConstraint : public Condition
{

public:
    typedef Condition BaseType;
    typedef BaseType::EquationIdVectorType EquationIdVectorType;
    typedef BaseType::MatrixType LHS_ContributionType;

    typedef Condition::GeometryType::Pointer PointerGeometryType;
    typedef Condition::GeometryType::CoordinatesArrayType CoordinatesArrayType;

    // Counted pointer of PointPointContactLink
    KRATOS_CLASS_POINTER_DEFINITION(NitscheIsotropicConstraint);

    /**
     * Default constructor.
     */
    NitscheIsotropicConstraint( IndexType NewId, GeometryType::Pointer pGeometry);

    NitscheIsotropicConstraint( IndexType NewId, GeometryType::Pointer pGeometry,
                           PropertiesType::Pointer pProperties);

    /**
     * Destructor.
     */
    virtual ~NitscheIsotropicConstraint();

    /**
     * Operations.
     */
    Condition::Pointer Create( IndexType NewId,
                               NodesArrayType const& ThisNodes,
                               PropertiesType::Pointer pProperties) const;

    Condition::Pointer Create( IndexType NewId,
                               GeometryType::Pointer pGeom,
                               PropertiesType::Pointer pProperties) const;

    virtual void Initialize(const ProcessInfo& rCurrentProcessInfo);

    virtual void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                               VectorType& rRightHandSideVector,
                               ProcessInfo& rCurrentProcessInfo);

    virtual void CalculateRightHandSide( VectorType& rRightHandSideVector,
                                 ProcessInfo& rCurrentProcessInfo);

    virtual void EquationIdVector( EquationIdVectorType& rResult,
                           ProcessInfo& rCurrentProcessInfo);

    virtual void GetDofList( DofsVectorType& ConditionalDofList,
                     ProcessInfo& CurrentProcessInfo);

    virtual int Check( const ProcessInfo& rCurrentProcessInfo );

private:

    IntegrationMethod mThisIntegrationMethod;

    void CalculateAll( MatrixType& rLeftHandSideMatrix,
                       VectorType& rRightHandSideVector,
                       ProcessInfo& rCurrentProcessInfo,
                       bool CalculateStiffnessMatrixFlag,
                       bool CalculateResidualVectorFlag);

    void CalculateInterpolationOperator( const unsigned int& Dim, Matrix& N_Operator, Vector& Ncontainer );

    void CalculateNoperator( const unsigned int& Dim, Matrix& N_Operator, const Vector& n );

    void CalculateBoperator( Matrix& B_Operator, const Matrix& DN_DX );

    void CalculateElasticMatrix2DPlaneStrain( Matrix& C, const double& E, const double& NU );

    void CalculateElasticMatrix3D( Matrix& C, const double& E, const double& NU );

    friend class Serializer;

    // A private default constructor necessary for serialization
    NitscheIsotropicConstraint() {};

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
    }

}; // Class NitscheIsotropicConstraint
}  // namespace Kratos.

#endif // KRATOS_NITSCHE_ISOTROPIC_CONSTRAINT_CONDITION_H_INCLUDED  defined


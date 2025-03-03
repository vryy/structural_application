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
*   Date:                $Date: 17 Oct 2017 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

#if !defined(KRATOS_COLLINEAR_CONSTRAINT_CONDITION_H_INCLUDED )
#define  KRATOS_COLLINEAR_CONSTRAINT_CONDITION_H_INCLUDED



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
 * Collinear constraint to align 3 points in 2D.
 * This constraint ensures 3 collinear points will remain collinear.
 * Because the constraint is nonlinear, the Lagrange multiplier method is used. The Lagrange multiplier is stored at middle node.
 * The line description:
 *  (A)-----(C)-----(B)
 * The constraint:
 *  (yC-yA)/(xC-xA) = (yB-yA)/(xB-xA)  <=>  (xA*yB + xB*yC + xC*yA) - (xB*yA + xC*yB + xA*yC) = 0
 * The weak form:
 *  r = lambda*[(xA*yB + xB*yC + xC*yA) - (xB*yA + xC*yB + xA*yC)]
 * node C will carry the additional Lagrange multiplier dof
 */
class CollinearConstraint2D : public Condition
{

public:
    typedef Condition BaseType;
    typedef BaseType::EquationIdVectorType EquationIdVectorType;

    KRATOS_CLASS_POINTER_DEFINITION(CollinearConstraint2D);

    /**
     * Default constructor.
     */
    CollinearConstraint2D( IndexType NewId, GeometryType::Pointer pGeometry);

    CollinearConstraint2D( IndexType NewId, GeometryType::Pointer pGeometry,
                           PropertiesType::Pointer pProperties);

    /**
     * Custom constructor.
     */
    CollinearConstraint2D( IndexType NewId, Node<3>::Pointer const& nodeA, Node<3>::Pointer const& nodeC, Node<3>::Pointer const& nodeB, PropertiesType::Pointer pProperties );

    CollinearConstraint2D( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties );

    /**
     * Destructor.
     */
    virtual ~CollinearConstraint2D();

    /**
     * Operations.
     */
    Condition::Pointer Create( IndexType NewId,
                               NodesArrayType const& ThisNodes,
                               PropertiesType::Pointer pProperties ) const final;

    Condition::Pointer Create( IndexType NewId,
                               GeometryType::Pointer pGeom,
                               PropertiesType::Pointer pProperties ) const final;

    /**
     * Calculates the local system contributions for this contact element
     */
    void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                               VectorType& rRightHandSideVector,
                               const ProcessInfo& rCurrentProcessInfo ) final;

    void CalculateRightHandSide( VectorType& rRightHandSideVector,
                                 const ProcessInfo& rCurrentProcessInfo ) final;

    void EquationIdVector( EquationIdVectorType& rResult,
                           const ProcessInfo& rCurrentProcessInfo ) const final;

    void GetDofList( DofsVectorType& ConditionalDofList,
                     const ProcessInfo& CurrentProcessInfo ) const final;

    std::string Info() const final
    {
        return "CollinearConstraint";
    }

private:

    void CalculateAll( MatrixType& rLeftHandSideMatrix,
                       VectorType& rRightHandSideVector,
                       const ProcessInfo& rCurrentProcessInfo,
                       bool CalculateStiffnessMatrixFlag,
                       bool CalculateResidualVectorFlag );

    friend class Serializer;

    // A private default constructor necessary for serialization
    CollinearConstraint2D() {};

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
    }

}; // Class CollinearConstraint2D

}  // namespace Kratos.

#endif // KRATOS_COLLINEAR_CONSTRAINT_CONDITION_H_INCLUDED  defined

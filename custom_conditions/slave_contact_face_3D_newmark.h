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
*   Last Modified by:    $Author: janosch $
*   Date:                $Date: 2007-03-13 15:01:51 $
*   Revision:            $Revision: 1.5 $
*
* ***********************************************************/

#if !defined(KRATOS_SLAVE_CONTACT_FACE_3D_NEWMARK_CONDITION_H_INCLUDED )
#define  KRATOS_SLAVE_CONTACT_FACE_3D_NEWMARK_CONDITION_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "containers/pointer_vector_set.h"

#include "custom_conditions/master_contact_face_3D.h"


namespace Kratos
{

/**
 * Contact surface element for 3D contact problems
 * Defines a facet of a 3D-Element as a contact surface for
 * slave contact surfaces
 * adapted from face2D originally written by riccardo.
 */
class SlaveContactFace3DNewmark : public Condition
{
public:
    // Counted pointer of SlaveContactFace3DNewmark
    KRATOS_CLASS_POINTER_DEFINITION(SlaveContactFace3DNewmark);

    typedef Condition BaseType;
    typedef BaseType::GeometryType GeometryType;
    typedef PointerVectorSet<MasterContactFace3D, IndexedObject> ContactMasterContainerType;
    typedef GeometryType::PointsArrayType PointsArrayType;
    typedef BaseType::EquationIdVectorType EquationIdVectorType;
    typedef std::size_t IndexType;
    #ifdef SD_APP_FORWARD_COMPATIBILITY
    typedef Point PointType;
    #else
    typedef Point<3> PointType;
    #endif

    /**
     * REMOVED
     */
    //     typedef BaseType::SecondaryCondition SecondaryConditionType;
    //     typedef PointerVectorSet< SecondaryConditionType, IndexedObject> SecondaryConditionContainerType;

    /**
     * Default constructor.
     */
    SlaveContactFace3DNewmark( IndexType NewId, GeometryType::Pointer pGeometry);
    SlaveContactFace3DNewmark( IndexType NewId, GeometryType::Pointer pGeometry,
                               PropertiesType::Pointer pProperties);
    /**
     * Destructor.
     */
    virtual ~SlaveContactFace3DNewmark();

    /**
     * Operations.
     */

//             Matrix TangentialVectors( GeometryType::CoordinatesArrayType& rPoint );

//             Vector NormalVector( GeometryType::CoordinatesArrayType& rPoint );

    Condition::Pointer Create( IndexType NewId,
                               NodesArrayType const& ThisNodes,
                               PropertiesType::Pointer pProperties) const;

    Condition::Pointer Create( IndexType NewId,
                               GeometryType::Pointer pGeom,
                               PropertiesType::Pointer pProperties) const;

    /**
     * searches the contact partner for a given integration point of the current
     * slave surface
     * @param AllMasterElements a Set of possible master Elements for the current
     * slave surface
     * @param IntegrationPointIndex the index of the current integration point
     * @return true if current slave surface has found a partner in each
     * integration point
     */
    bool SearchPartner( ContactMasterContainerType& AllMasterElements,
                        IndexType IntegrationPointIndex,
                        PointType MasterContactLocalPoint,
                        PointType SlaveContactLocalPoint,
                        Condition::Pointer CurrentMaster,
                        double CurrentLambda
                      );

    /**
     * updates current lagrange multiplier
     * @return a Matrix containing each new Lambda and DeltaLambda for each
     * integration point
     */
    Matrix UpdateLambda();


    /**
     * Calculates the local system contributions for this contact element
     */
    void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                               VectorType& rRightHandSideVector,
                               const ProcessInfo& rCurrentProcessInfo);

    /**
     * REMOVED
     */
//             void CalculateCrossElementarySystemContributions( SecondaryConditionContainerType& SecondaryConditions, EquationIdVectorType& PrimaryEquationId, const ProcessInfo& rCurrentProcessInfo );

    /**
     * calculates the stiffness matrix contributions resulting from linearization
     * of normal vector
     * (REMOVED)
     */
//             void CalculateNormalLinearizationElementarySystemContributions( SecondaryConditionContainerType& SecondaryConditions, EquationIdVectorType& PrimaryEquationId, const ProcessInfo& rCurrentProcessInfo );


    void CalculateRightHandSide( VectorType& rRightHandSideVector,
                                 const ProcessInfo& rCurrentProcessInfo);

    void EquationIdVector( EquationIdVectorType& rResult,
                           const ProcessInfo& rCurrentProcessInfo) const;

    void GetDofList( DofsVectorType& ConditionalDofList,
                     const ProcessInfo& CurrentProcessInfo) const;

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& Output, const ProcessInfo& rCurrentProcessInfo);

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
    void CalculateAll( MatrixType& rLeftHandSideMatrix,
                       VectorType& rRightHandSideVector,
                       const ProcessInfo& rCurrentProcessInfo,
                       bool CalculateStiffnessMatrixFlag,
                       bool CalculateResidualVectorFlag);

    void CalculateAndAddKc( Matrix& K,
                            const Vector& N,
                            double weight,
                            double dA,
                            Vector v
                          );

    void CalculateAndAdd_PressureForce( Vector& residualvector,
                                        const Vector& N,
                                        Vector& v3,
                                        double pressure,
                                        double weight,
                                        double DetJ );

    Matrix TangentialVectors_inOrigin( const GeometryType::CoordinatesArrayType& rPoint );

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization
    SlaveContactFace3DNewmark() {};

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
    }

    /**
     * Assignment operator.
     * (DEACTIVATED)
     */
    //SlaveContactFace3DNewmark& operator=(const SlaveContactFace3DNewmark& rOther);

    /**
     * Copy constructor.
     * (DEACTIVATED)
     */
    //SlaveContactFace3DNewmark(const SlaveContactFace3DNewmark& rOther);

    ContactMasterContainerType::Pointer mpMasterElements;
    PointsArrayType mpContactPartnersGlobal;
    PointsArrayType mpContactPartnersLocal;
//             VectorType mLambdas;
//             VectorType mGaps;
    double penalty;
}; // Class SlaveContactFace3DNewmark
}  // namespace Kratos.

#endif // KRATOS_SLAVE_CONTACT_FACE_3D_CONDITION_H_INCLUDED  defined

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 9 Aug 2018 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_ELASTIC_FACE_SPRINGS_CONDITION_H_INCLUDED )
#define  KRATOS_ELASTIC_FACE_SPRINGS_CONDITION_H_INCLUDED



// System includes


// External includes


// Project includes
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"


namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/**
 * Surface springs normal to the surface
*/
class KRATOS_API(STRUCTURAL_APPLICATION) ElasticFaceSprings : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of ElasticFaceSprings
    KRATOS_CLASS_POINTER_DEFINITION(ElasticFaceSprings);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ElasticFaceSprings(IndexType NewId, GeometryType::Pointer pGeometry);
    ElasticFaceSprings(IndexType NewId, GeometryType::Pointer pGeometry,
                 PropertiesType::Pointer pProperties);

    /// Special Constructor for elastic point constraint
    ElasticFaceSprings( IndexType NewId, Node<3>::Pointer const& pNode, PropertiesType::Pointer pProperties );

    /// Destructor.
    ~ElasticFaceSprings() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    IntegrationMethod GetIntegrationMethod() const override;

    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo ) override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(DofsVectorType& rConditionalDofList, const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ElasticFaceSprings";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        BaseType::PrintInfo(rOStream);
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {}

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    friend class Serializer;

    // A private default constructor necessary for serialization
    ElasticFaceSprings() {};

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
    }

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                       const ProcessInfo& rCurrentProcessInfo, bool CalculateStiffnessMatrixFlag,
                       bool CalculateResidualVectorFlag );

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //ElasticFaceSprings& operator=(const ElasticFaceSprings& rOther);

    /// Copy constructor.
    //ElasticFaceSprings(const ElasticFaceSprings& rOther);


    ///@}

}; // Class ElasticFaceSprings

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

}  // namespace Kratos.

#endif // KRATOS_ELASTIC_FACE_SPRINGS_CONDITION_H_INCLUDED  defined

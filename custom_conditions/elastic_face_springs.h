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
#include "includes/node.h"
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
template<typename TNodeType>
class KRATOS_API(STRUCTURAL_APPLICATION) BaseElasticFaceSprings : public BaseCondition<TNodeType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of BaseElasticFaceSprings
    KRATOS_CLASS_POINTER_DEFINITION(BaseElasticFaceSprings);

    typedef BaseCondition<TNodeType> BaseType;

    typedef typename BaseType::ConditionType ConditionType;

    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef typename BaseType::NodesArrayType NodesArrayType;

    typedef typename BaseType::EquationIdVectorType EquationIdVectorType;

    typedef typename BaseType::DofsVectorType DofsVectorType;

    typedef typename BaseType::IndexType IndexType;

    typedef typename BaseType::DataType DataType;

    typedef typename BaseType::ValueType ValueType;

    typedef typename BaseType::VectorType VectorType;

    typedef typename BaseType::ZeroVectorType ZeroVectorType;

    typedef typename BaseType::MatrixType MatrixType;

    typedef typename BaseType::ZeroMatrixType ZeroMatrixType;

    typedef typename BaseType::GeometryType GeometryType;

    typedef typename BaseType::PropertiesType PropertiesType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BaseElasticFaceSprings(IndexType NewId, typename GeometryType::Pointer pGeometry);
    BaseElasticFaceSprings(IndexType NewId, typename GeometryType::Pointer pGeometry,
                           typename PropertiesType::Pointer pProperties);

    /// Special Constructor for elastic point constraint
    BaseElasticFaceSprings( IndexType NewId, typename TNodeType::Pointer const& pNode, typename PropertiesType::Pointer pProperties );

    /// Destructor.
    ~BaseElasticFaceSprings() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    IntegrationMethod GetIntegrationMethod() const override;

    typename ConditionType::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, typename PropertiesType::Pointer pProperties) const override;

    typename ConditionType::Pointer Create(IndexType NewId, typename GeometryType::Pointer pGeom, typename PropertiesType::Pointer pProperties) const override;

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
    BaseElasticFaceSprings() {};

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType );
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType );
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
    //BaseElasticFaceSprings& operator=(const BaseElasticFaceSprings& rOther);

    /// Copy constructor.
    //BaseElasticFaceSprings(const BaseElasticFaceSprings& rOther);


    ///@}

}; // Class BaseElasticFaceSprings

///@}

///@name Type Definitions
///@{

typedef BaseElasticFaceSprings<RealNode> ElasticFaceSprings;
typedef BaseElasticFaceSprings<ComplexNode> ComplexElasticFaceSprings;
typedef BaseElasticFaceSprings<GComplexNode> GComplexElasticFaceSprings;

///@}
///@name Input and output
///@{


///@}

}  // namespace Kratos.

#endif // KRATOS_ELASTIC_FACE_SPRINGS_CONDITION_H_INCLUDED  defined

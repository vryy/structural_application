//
//   Project Name:        Kratos
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2007-12-12 14:51:07 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_POINT_FORCE_CONDITION_H_INCLUDED )
#define  KRATOS_POINT_FORCE_CONDITION_H_INCLUDED



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
/** Detail class definition.
*/
template<int TDim, typename TNodeType>
class PointForce : public BaseCondition<TNodeType>
{
public:
    ///@name Type Definitions
    ///@{

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

    /// Counted pointer of PointForce3D
    KRATOS_CLASS_POINTER_DEFINITION(PointForce);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PointForce(IndexType NewId, typename GeometryType::Pointer pGeometry);
    PointForce(IndexType NewId, typename GeometryType::Pointer pGeometry,
                     typename PropertiesType::Pointer pProperties);

    PointForce( IndexType NewId, typename GeometryType::PointType::Pointer const& pNode, typename PropertiesType::Pointer pProperties );

    /// Destructor.
    ~PointForce() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                      typename PropertiesType::Pointer pProperties) const final;

    typename BaseType::Pointer Create(IndexType NewId, typename GeometryType::Pointer pGeom,
                                      typename PropertiesType::Pointer pProperties) const final;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType&
                              rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) final;

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                const ProcessInfo& rCurrentProcessInfo) final;

    void CalculateDampingMatrix( MatrixType& rDampMatrix,
                                 const ProcessInfo& rCurrentProcessInfo) final;

    void CalculateMassMatrix( MatrixType& rMassMatrix,
                              const ProcessInfo& rCurrentProcessInfo) final;

    void EquationIdVector(EquationIdVectorType& rResult,
                          const ProcessInfo& rCurrentProcessInfo) const final;

    void GetDofList(DofsVectorType& ConditionalDofList,
                    const ProcessInfo& CurrentProcessInfo) const final;

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
    std::string Info() const final
    {
        return "PointForce3D";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const final
    {
        rOStream << "PointForce3D #" << this->Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const final
    {
        BaseType::PrintData(rOStream);
    }

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
    PointForce() {};

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
    //PointForce& operator=(const PointForce& rOther);

    /// Copy constructor.
    //PointForce(const PointForce& rOther);

    ///@}

}; // Class PointForce

///@}

///@name Type Definitions
///@{

typedef PointForce<2, RealNode> PointForce2D;
typedef PointForce<2, ComplexNode> ComplexPointForce2D;
typedef PointForce<2, GComplexNode> GComplexPointForce2D;
typedef PointForce<3, RealNode> PointForce3D;
typedef PointForce<3, ComplexNode> ComplexPointForce3D;
typedef PointForce<3, GComplexNode> GComplexPointForce3D;

///@}
///@name Input and output
///@{


///@}

}  // namespace Kratos.

#endif // KRATOS_POINT_FORCE_CONDITION_H_INCLUDED  defined

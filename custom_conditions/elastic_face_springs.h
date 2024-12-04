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
    virtual ~ElasticFaceSprings();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                PropertiesType::Pointer pProperties) const;

    Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom,
                                PropertiesType::Pointer pProperties) const;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType&
                              rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo&
                                rCurrentProcessInfo);
    //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo);

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&
                          rCurrentProcessInfo) const;

    void GetDofList(DofsVectorType& ConditionalDofList, const ProcessInfo&
                    CurrentProcessInfo) const;

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
//      virtual String Info() const;

    /// Print information about this object.
//      virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
//      virtual void PrintData(std::ostream& rOStream) const;


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

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
    }

    virtual void load(Serializer& rSerializer)
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


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    ElasticFaceSprings& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const ElasticFaceSprings& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_ELASTIC_FACE_SPRINGS_CONDITION_H_INCLUDED  defined




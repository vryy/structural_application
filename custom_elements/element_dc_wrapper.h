//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 15 Sep 2024 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ELEMENT_DC_WRAPPER_H_INCLUDED )
#define  KRATOS_ELEMENT_DC_WRAPPER_H_INCLUDED

// Project includes
#include "includes/element.h"

namespace Kratos
{

/**
 * Wrapper to use the element in Displacement (DC) Control analysis.
 * In different to load control, the material update is performed during FinalizeNonLinearIteration.
 * Hence, the tangent in the first iteration is from the last step. It improves
 * convergence particularly for DC analysis.
 * This element must not be used with a DC version of the constitutive law.
 */
template<class TElementType>
class ElementDcWrapper : public TElementType
{
public:
    ///@name Type Definitions
    ///@{

    typedef TElementType BaseType;
    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::GeometryType GeometryType;
    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::PropertiesType PropertiesType;
    typedef typename BaseType::VectorType VectorType;
    typedef typename BaseType::MatrixType MatrixType;

    KRATOS_CLASS_POINTER_DEFINITION( ElementDcWrapper );

    ///@}
    ///@name Life Cycle
    ///@{

    // Reused constructors
    using BaseType::BaseType;

    /// Destructor.
    ~ElementDcWrapper() override
    {}

    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create( IndexType NewId, NodesArrayType const& ThisNodes, typename PropertiesType::Pointer pProperties ) const override
    {
        return Element::Pointer( new ElementDcWrapper( NewId, this->GetGeometry().Create( ThisNodes ), pProperties ) );
    }

    Element::Pointer Create( IndexType NewId, typename GeometryType::Pointer pGeom, typename PropertiesType::Pointer pProperties ) const override
    {
        return Element::Pointer( new ElementDcWrapper( NewId, pGeom, pProperties ) );
    }

#ifndef SD_APP_FORWARD_COMPATIBILITY
    Element::Pointer Create( IndexType NewId, std::vector<typename GeometryType::Pointer> pGeometries, typename PropertiesType::Pointer pProperties ) const override
    {
        return Element::Pointer( new ElementDcWrapper( NewId, pGeometries, pProperties ) );
    }
#endif

    void Initialize( const ProcessInfo& rCurrentProcessInfo ) override
    {
        BaseType::Initialize( rCurrentProcessInfo );

        typename BaseType::EquationIdVectorType EquationIds;
        this->EquationIdVector( EquationIds, rCurrentProcessInfo );

        const unsigned int mat_size = EquationIds.size();
        mRightHandSideVector.resize(mat_size, false);
        mRightHandSideVector.clear();
    }

    void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo ) override
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = true;
        bool CalculateResidualVectorFlag = false;

        this->CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                            CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );

        if (rRightHandSideVector.size() != mRightHandSideVector.size())
            rRightHandSideVector.resize(mRightHandSideVector.size());
        noalias(rRightHandSideVector) = mRightHandSideVector;
    }

    void FinalizeNonLinearIteration( const ProcessInfo& rCurrentProcessInfo ) override
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = false;
        bool CalculateResidualVectorFlag = true;

        Matrix dummy;
        this->CalculateAll( dummy, mRightHandSideVector, rCurrentProcessInfo,
                            CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );

        BaseType::FinalizeNonLinearIteration(rCurrentProcessInfo);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ElementDcWrapper<" + BaseType::Info() + ">";
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    Vector mRightHandSideVector;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
        rSerializer.save("mRightHandSideVector", mRightHandSideVector);
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
        rSerializer.load("mRightHandSideVector", mRightHandSideVector);
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    // A private default constructor necessary for serialization
    ElementDcWrapper() {}

    ///@}

}; // Class ElementDcWrapper

}  // namespace Kratos.

#endif // KRATOS_ELEMENT_DC_WRAPPER_H_INCLUDED defined

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
giang.bui@rub.de
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
//
//   Project Name:        Kratos
//   Modified by:         $Author: hbui $
//   Date:                $Date: 2 May 2025 $
//
//


#if !defined(KRATOS_KINEMATIC_LINEAR_ANTI_PLANE_INCLUDED )
#define  KRATOS_KINEMATIC_LINEAR_ANTI_PLANE_INCLUDED



// System includes


// External includes


// Project includes
#include "custom_elements/kinematic_linear.h"


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
 * Define a small strain element with strain measure as Infinitesimal strain and stress measure as Cauchy stress for anti-plane stress state
 * The anti-plane motion is characterized by only vertical displacement, hence this element only carries Uz dofs.
 * Reference:
 * +    https://chatgpt.com/share/68148590-57d8-8006-a0df-9aee9d24bb66
 */
template<class TNodeType>
class KRATOS_API(STRUCTURAL_APPLICATION) BaseKinematicLinearAntiPlane : public BaseKinematicLinear<TNodeType>
{

public:
    ///@name Type Definitions
    ///@{

    typedef BaseKinematicLinear<TNodeType> BaseType;

    typedef typename BaseType::ElementType ElementType;

    typedef typename BaseType::IndexType IndexType;

    typedef typename BaseType::IndexType SizeType;

    typedef typename BaseType::DataType DataType;

    typedef typename BaseType::ValueType ValueType;

    typedef typename BaseType::NodesArrayType NodesArrayType;

    typedef typename BaseType::GeometryType GeometryType;

    typedef typename BaseType::PropertiesType PropertiesType;

    typedef typename BaseType::EquationIdVectorType EquationIdVectorType;

    typedef typename BaseType::DofsVectorType DofsVectorType;

    typedef typename BaseType::VectorType VectorType;
    typedef typename BaseType::ZeroVectorType ZeroVectorType;
    typedef typename BaseType::MatrixType MatrixType;
    typedef typename BaseType::ZeroMatrixType ZeroMatrixType;

    KRATOS_CLASS_POINTER_DEFINITION( BaseKinematicLinearAntiPlane );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BaseKinematicLinearAntiPlane( IndexType NewId, typename GeometryType::Pointer pGeometry );
    BaseKinematicLinearAntiPlane( IndexType NewId, typename GeometryType::Pointer pGeometry, typename PropertiesType::Pointer pProperties );

    /// Destructor.
    ~BaseKinematicLinearAntiPlane() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    typename ElementType::Pointer Create( IndexType NewId, NodesArrayType const& ThisNodes, typename PropertiesType::Pointer pProperties ) const override;

    typename ElementType::Pointer Create( IndexType NewId, typename GeometryType::Pointer pGeom, typename PropertiesType::Pointer pProperties ) const override;

    void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateRightHandSide( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateMassMatrix( MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateDampingMatrix( MatrixType& rDampMatrix, const ProcessInfo& rCurrentProcessInfo ) override;

    void EquationIdVector( EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo ) const override;

    void GetDofList( DofsVectorType& ElementalDofList, const ProcessInfo& CurrentProcessInfo ) const override;

    void GetValuesVector( VectorType& values, int Step = 0 ) const override;

    void GetFirstDerivativesVector( VectorType& values, int Step = 0 ) const override;

    void GetSecondDerivativesVector( VectorType& values, int Step = 0 ) const override;

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
        return "BaseKinematicLinearAntiPlane";
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

    void CalculateBoperator( MatrixType& B_Operator, const Vector& N, const MatrixType& DN_DX ) const override;

    void CalculateBBaroperator( MatrixType& B_Operator, const MatrixType& DN_DX, const MatrixType& Bdil_bar ) const override;

    SizeType WorkingSpaceDimension() const override
    {
        return 3; // this is 2D element but working space dimension is 3
    }

    unsigned int GetStrainSize( const unsigned int dim ) const override
    {
        return 6; // although there are only two nonzero strain components. The full strain is needed to work with 3D constitutive law.
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
    ///@name Serialization
    ///@{

    friend class Serializer;

    // A private default constructor necessary for serialization
    BaseKinematicLinearAntiPlane() {}

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //BaseKinematicLinearAntiPlane& operator=(const BaseKinematicLinearAntiPlane& rOther);

    /// Copy constructor.
    //BaseKinematicLinearAntiPlane(const BaseKinematicLinearAntiPlane& rOther);

    ///@}

}; // Class BaseKinematicLinearAntiPlane

///@}
///@name Type Definitions
///@{

typedef BaseKinematicLinearAntiPlane<RealNode> KinematicLinearAntiPlane;
typedef BaseKinematicLinearAntiPlane<ComplexNode> ComplexKinematicLinearAntiPlane;
typedef BaseKinematicLinearAntiPlane<GComplexNode> GComplexKinematicLinearAntiPlane;

///@}
///@name Input and output
///@{


///@}

}  // namespace Kratos.

#endif // KRATOS_KINEMATIC_LINEAR_ANTI_PLANE_INCLUDED defined

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
//   Date:                $Date: 11 Oct 2023 $
//
//


#if !defined(KRATOS_FINITE_STRAIN_AXISYMMETRIC_INCLUDED )
#define  KRATOS_FINITE_STRAIN_AXISYMMETRIC_INCLUDED



// System includes


// External includes


// Project includes
#include "custom_elements/finite_strain.h"


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
 * Extend the finite strain element for Axisymmetric problem
 * Reference:
 * +    Ted Belytschko, Nonlinear Finite Elements for Continua and Structures, 2014, (E4.4.6)
 */
class KRATOS_API(STRUCTURAL_APPLICATION) FiniteStrainAxisymmetric : public FiniteStrain
{

public:
    ///@name Type Definitions
    ///@{

    typedef FiniteStrain BaseType;

    KRATOS_CLASS_POINTER_DEFINITION( FiniteStrainAxisymmetric );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    FiniteStrainAxisymmetric( IndexType NewId, GeometryType::Pointer pGeometry );
    FiniteStrainAxisymmetric( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    /// Destructor.
    virtual ~FiniteStrainAxisymmetric();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const override;

    Element::Pointer Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const override;

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
        return "FiniteStrainAxisymmetric";
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

    unsigned int GetStrainSize( unsigned int dim ) const override { return 4; }

    unsigned int GetFSize( unsigned int dim ) const override { return 3; }

    unsigned int GetGSize( unsigned int dim ) const override { return 5; }

    void CalculateB( Matrix& B_Operator, const VectorType& N, const Matrix& DN_DX, const Matrix& CurrentDisp ) const override;

    void CalculateG( Matrix& G_Operator, const VectorType& N, const Matrix& DN_DX ) const override;

    void CalculateG( Matrix& G_Operator, const VectorType& N, const Matrix& DN_DX, const Matrix& CurrentDisp ) const override;

    void CalculateF( Matrix& F, const Matrix& G_Operator, const Matrix& CurrentDisp ) const override;

    void CalculateStrain( const MatrixType& B, Vector& StrainVector ) const override;

    double GetIntegrationWeight( double Weight, const VectorType& N, const Matrix& CurrentDisp ) const override;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    // A private default constructor necessary for serialization
    FiniteStrainAxisymmetric() {}

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer,  BaseType );
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer,  BaseType );
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
    //FiniteStrainAxisymmetric& operator=(const FiniteStrainAxisymmetric& rOther);

    /// Copy constructor.
    //FiniteStrainAxisymmetric(const FiniteStrainAxisymmetric& rOther);


    ///@}

}; // Class FiniteStrainAxisymmetric

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

}  // namespace Kratos.

#endif // KRATOS_FINITE_STRAIN_AXISYMMETRIC_INCLUDED defined

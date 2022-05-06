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
//   Date:                $Date: 14 Feb 2022 $
//
//


#if !defined(KRATOS_KINEMATIC_LINEAR_AXISYMMETRIC_INCLUDED )
#define  KRATOS_KINEMATIC_LINEAR_AXISYMMETRIC_INCLUDED



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
Define a small strain element with strain measure as Infinitesimal strain and stress measure as Cauchy stress for Axisymmetric problem
 */
class KinematicLinearAxisymmetric : public KinematicLinear
{

public:
    ///@name Type Definitions
    ///@{

    typedef KinematicLinear BaseType;

    KRATOS_CLASS_POINTER_DEFINITION( KinematicLinearAxisymmetric );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KinematicLinearAxisymmetric( IndexType NewId, GeometryType::Pointer pGeometry );
    KinematicLinearAxisymmetric( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    /// Destructor.
    virtual ~KinematicLinearAxisymmetric();


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
        return "KinematicLinearAxisymmetric";
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

    void CalculateBoperator( Matrix& B_Operator, const Vector& N, const Matrix& DN_DX ) const override;

    void CalculateBBaroperator( Matrix& B_Operator, const Matrix& DN_DX, const Matrix& Bdil_bar ) const override;

    unsigned int GetStrainSize( const unsigned int& dim ) const override
    {
        return 4;
    }

    double GetIntegrationWeight( const GeometryType::IntegrationPointsArrayType& integration_points,
            const unsigned int& PointNumber, const Matrix& Ncontainer ) const override;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    // A private default constructor necessary for serialization
    KinematicLinearAxisymmetric() {}

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
    //KinematicLinearAxisymmetric& operator=(const KinematicLinearAxisymmetric& rOther);

    /// Copy constructor.
    //KinematicLinearAxisymmetric(const KinematicLinearAxisymmetric& rOther);


    ///@}

}; // Class KinematicLinearAxisymmetric

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
                   KinematicLinearAxisymmetric& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
                   const KinematicLinearAxisymmetric& rThis)
           {
                   rThis.PrintInfo(rOStream);
                   rOStream << std::endl;
                   rThis.PrintData(rOStream);

                   return rOStream;
}*/
///@}

}  // namespace Kratos.

#endif // KRATOS_KINEMATIC_LINEAR_AXISYMMETRIC_INCLUDED defined



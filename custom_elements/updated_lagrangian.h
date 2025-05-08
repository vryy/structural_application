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
//   Date:                $Date: 14 Sep 2023 $
//
//


#if !defined(KRATOS_UPDATED_LAGRANGIAN_INCLUDED )
#define  KRATOS_UPDATED_LAGRANGIAN_INCLUDED



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
 * This extension of the FiniteStrain element calculates the Jacobian in relative to the current position of the node.
 * Therefore, it allows for calculation with updated Jacobian. This can be controlled by updating the position
 * of the node via MoveMeshFlag.
 * Furthermore, the constitutive law must be formulated to use incremental multiplicative deformation gradient.
 */
class UpdatedLagrangian : public FiniteStrain
{

public:
    ///@name Type Definitions
    ///@{

    typedef FiniteStrain BaseType;

    KRATOS_CLASS_POINTER_DEFINITION( UpdatedLagrangian );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    UpdatedLagrangian( IndexType NewId, GeometryType::Pointer pGeometry );
    UpdatedLagrangian( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    /// Destructor.
    virtual ~UpdatedLagrangian();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const override;

    Element::Pointer Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeSolutionStep(const ProcessInfo& CurrentProcessInfo) override;

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
        return "UpdatedLagrangian";
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

    void CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                      const ProcessInfo& rCurrentProcessInfo,
                      bool CalculateStiffnessMatrixFlag,
                      bool CalculateResidualVectorFlag) override;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    // A private default constructor necessary for serialization
    UpdatedLagrangian() {}

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

    std::vector<MatrixType> mCurrentF;
    std::vector<MatrixType> mLastF;

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
    //UpdatedLagrangian& operator=(const UpdatedLagrangian& rOther);

    /// Copy constructor.
    //UpdatedLagrangian(const UpdatedLagrangian& rOther);


    ///@}

}; // Class UpdatedLagrangian

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif // KRATOS_UPDATED_LAGRANGIAN_INCLUDED defined



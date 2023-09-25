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
//   Date:                $Date: 13 Sep 2023 $
//
//


#if !defined(KRATOS_UPDATED_KINEMATIC_LINEAR_INCLUDED )
#define  KRATOS_UPDATED_KINEMATIC_LINEAR_INCLUDED



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
 * This extension of the KinematicLinear element calculates the Jacobian using the current position of the node.
 * Therefore, it allows for calculation with updated Jacobian. This can be controlled by updating the position
 * of the node via MoveMeshFlag. It is noted that, the update of the nodal coordinates should be after the step is converged.
 * Otherwise, it may lead to nonconvergence since the Jacobian is not frozen during step iteration.
 * If MoveMeshFlag is not used, this element will behave similar to KinematicLinear element.
 * Furthermore, the constitutive law must be formulated to use incremental strain instead of total strain.
 */
class UpdatedKinematicLinear : public KinematicLinear
{

public:
    ///@name Type Definitions
    ///@{

    typedef KinematicLinear BaseType;

    KRATOS_CLASS_POINTER_DEFINITION( UpdatedKinematicLinear );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    UpdatedKinematicLinear( IndexType NewId, GeometryType::Pointer pGeometry );
    UpdatedKinematicLinear( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    /// Destructor.
    virtual ~UpdatedKinematicLinear();


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
        return "UpdatedKinematicLinear";
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

    void CalculateJacobian( GeometryType::JacobiansType& J0 ) const override;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    // A private default constructor necessary for serialization
    UpdatedKinematicLinear() {}

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
    //UpdatedKinematicLinear& operator=(const UpdatedKinematicLinear& rOther);

    /// Copy constructor.
    //UpdatedKinematicLinear(const UpdatedKinematicLinear& rOther);


    ///@}

}; // Class UpdatedKinematicLinear

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
                   UpdatedKinematicLinear& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
                   const UpdatedKinematicLinear& rThis)
           {
                   rThis.PrintInfo(rOStream);
                   rOStream << std::endl;
                   rThis.PrintData(rOStream);

                   return rOStream;
}*/
///@}

}  // namespace Kratos.

#endif // KRATOS_UPDATED_KINEMATIC_LINEAR_INCLUDED defined



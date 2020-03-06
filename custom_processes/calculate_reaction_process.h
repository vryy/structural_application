/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 12 Jul 2017 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_CALCULATE_REACTION_PROCESS_H_INCLUDED )
#define  KRATOS_CALCULATE_REACTION_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "solving_strategies/schemes/scheme.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
//#include "spatial_containers/bins_static.h"
#include "spatial_containers/bins_dynamic.h"
#include "structural_application.h"

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

/// Class for calculating the internal forces
/** As the name said, this class supports for calculation of the nodal internal forces. If the d.o.f is fixed, the reaction value
 * is equivalent to the internal force value, otherwise it reflect the pure internal force of the d.o.f. This class must be called
 * before the final solve to compute the correct unbalanced forces.
 */
class CalculateReactionProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{
    typedef Element::GeometryType GeometryType;
    typedef GeometryType::PointType NodeType;
    typedef GeometryType::PointType::PointType PointType;
    typedef ModelPart::ElementsContainerType ElementsContainerType;
    typedef ModelPart::ConditionsContainerType ConditionsContainerType;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef Scheme<SparseSpaceType, LocalSpaceType> SchemeType;

    /// Pointer definition of CalculateReactionProcess
    KRATOS_CLASS_POINTER_DEFINITION(CalculateReactionProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /// This constructor will take all elements of the model_part for topology optimization
    CalculateReactionProcess(ModelPart& r_model_part, SchemeType& r_scheme)
    : mr_model_part(r_model_part), mr_scheme(r_scheme)
    {
    }

    /// Destructor.
    virtual ~CalculateReactionProcess()
    {
    }


    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{

    virtual void Execute()
    {
        Element::MatrixType LHS_Contribution = Element::MatrixType(0, 0);
        Element::VectorType RHS_Contribution = Element::VectorType(0);
        ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
        Element::EquationIdVectorType EquationId;
        Element::DofsVectorType ElementalDofList;
        std::size_t i;

        for(ElementsContainerType::ptr_iterator i_element = mr_model_part.Elements().ptr_begin();
                i_element != mr_model_part.Elements().ptr_end(); ++i_element)
        {
            if( (*i_element)->GetValue(IS_MARKED_FOR_REACTION) )
            {
                // get the list of elemental dofs
                (*i_element)->GetDofList(ElementalDofList, CurrentProcessInfo);

                // get the elemental rhs
                (*i_element)->CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);

                i = 0;
                for(typename Element::DofsVectorType::iterator i_dof = ElementalDofList.begin();
                        i_dof != ElementalDofList.end(); ++i_dof, ++i)
                {
                    (*i_dof)->GetSolutionStepReactionValue() -= RHS_Contribution[i];
                    // std::cout << "dof " << (*i_dof)->GetVariable().Name() << " of node " << (*i_dof)->Id() << " is added with " << RHS_Contribution[i] << std::endl;
                }

                // std::cout << "element " << (*i_element)->Id() << " reaction is computed, RHS_Contribution = " << RHS_Contribution << std::endl;

                // clean local elemental memory
                mr_scheme.CleanMemory(*i_element);
            }
        }

        for(ConditionsContainerType::ptr_iterator i_condition = mr_model_part.Conditions().ptr_begin();
                i_condition != mr_model_part.Conditions().ptr_end(); ++i_condition)
        {
            if( (*i_condition)->GetValue(IS_MARKED_FOR_REACTION) )
            {
                // get the list of elemental dofs
                (*i_condition)->GetDofList(ElementalDofList, CurrentProcessInfo);

                // get the elemental rhs
                (*i_condition)->CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);

                i = 0;
                for(typename Element::DofsVectorType::iterator i_dof = ElementalDofList.begin();
                        i_dof != ElementalDofList.end(); ++i_dof, ++i)
                {
                    (*i_dof)->GetSolutionStepReactionValue() -= RHS_Contribution[i];
                }

                // std::cout << "condition " << (*i_condition)->Id() << " reaction is computed, RHS_Contribution = " << RHS_Contribution << std::endl;

                // clean local elemental memory
                mr_scheme.CleanMemory(*i_condition);
            }
        }
    }

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
    virtual std::string Info() const
    {
        return "CalculateReactionProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "CalculateReactionProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
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

    ModelPart& mr_model_part;
    SchemeType& mr_scheme;

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
    CalculateReactionProcess& operator=(CalculateReactionProcess const& rOther);

    /// Copy constructor.
    //CalculateReactionProcess(FindConditionsNeighboursProcess const& rOther);


    ///@}

}; // Class CalculateReactionProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  CalculateReactionProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const CalculateReactionProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_CALCULATE_REACTION_PROCESS_H_INCLUDED  defined

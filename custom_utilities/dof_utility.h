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

/* *********************************************************
*
*   Last Modified by:    $Author: hbui $
*   Date:                $Date: 31 Mar 2016 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

#if !defined(KRATOS_DOF_UTILITY_H_INCLUDED )
#define  KRATOS_DOF_UTILITY_H_INCLUDED
//System includes

//External includes

//Project includes
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos
{

class DofUtility
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(DofUtility);

    /**
     * Constructor.
     */
    DofUtility()
    {
        std::cout << "DofUtility created" << std::endl;
    }

    /**
     * Destructor.
     */
    virtual ~DofUtility()
    {}

    /// List all the dofs in the model_part
    static void ListDofs(ModelPart::DofsArrayType& rDofSet, const std::size_t EquationSystemSize)
    {
        KRATOS_WATCH(rDofSet.size())
        KRATOS_WATCH(EquationSystemSize)

        // get the set of variables in the dof set and count them by type
        std::set<VariableData> VarSet;
        std::map<VariableData::KeyType, std::size_t> VarCount;
        for(ModelPart::DofsArrayType::iterator dof_iterator = rDofSet.begin(); dof_iterator != rDofSet.end(); ++dof_iterator)
        {
            if(dof_iterator->EquationId() < EquationSystemSize)
            {
                VarSet.insert(dof_iterator->GetVariable());
                ++VarCount[dof_iterator->GetVariable().Key()];
            }
        }

        std::cout << "List of Dofs in the system:" << std::endl;
        for(std::set<VariableData>::iterator it = VarSet.begin(); it != VarSet.end(); ++it)
        {
            std::cout << "    " << it->Name() << ": " << VarCount[it->Key()] << std::endl;
        }
    }

    /// Remove a dof from the list of dofs
    template<class TVariableType>
    static void RemoveDof(ModelPart::DofsArrayType& rDofSet, const TVariableType& rVariable)
    {
        for(ModelPart::DofsArrayType::iterator dof_iterator = rDofSet.begin(); dof_iterator != rDofSet.end();)
        {
            if(dof_iterator->GetVariable() == rVariable)
                dof_iterator = rDofSet.erase(dof_iterator);
            else
                ++dof_iterator;
        }
        std::cout << "Dof " << rVariable.Name() << " is removed from the dof set" << std::endl;
    }

    /// Print the key of a variable
    template<class TVariableType>
    static void PrintKey(const TVariableType& rThisVariable)
    {
        std::cout << "Variable " << rThisVariable.Name() << " key: " << rThisVariable.Key() << std::endl;
    }

    /// Find and print the dof containing the maximum unbalanced force
    static void PrintMaxUnbalancedForce(ModelPart::DofsArrayType& rDofSet, const Vector& forces)
    {
        double max_force = -1.0e99;
        std::size_t imax;
        for (std::size_t i = 0; i < forces.size(); ++i)
        {
            if (std::abs(forces[i]) > max_force)
            {
                max_force = std::abs(forces[i]);
                imax = i;
            }
        }

        for(ModelPart::DofsArrayType::iterator dof_iterator = rDofSet.begin(); dof_iterator != rDofSet.end(); ++dof_iterator)
        {
            if (dof_iterator->EquationId() == imax)
            {
                std::cout << "Max unbalanced force at node " << dof_iterator->Id()
                          << ", dof " << imax << " (" << dof_iterator->GetVariable().Name() << "), value = " << max_force
                          << std::endl;
            }
        }
    }

    /// Find and print the dof containing the minimum unbalanced force
    static void PrintMinUnbalancedForce(ModelPart::DofsArrayType& rDofSet, const Vector& forces)
    {
        double min_force = 1.0e99;
        std::size_t imin;
        for (std::size_t i = 0; i < forces.size(); ++i)
        {
            const auto force_abs = std::abs(forces[i]);
            if (force_abs < min_force && force_abs > 0.0)
            {
                min_force = force_abs;
                imin = i;
            }
        }

        for(ModelPart::DofsArrayType::iterator dof_iterator = rDofSet.begin(); dof_iterator != rDofSet.end(); ++dof_iterator)
        {
            if (dof_iterator->EquationId() == imin)
            {
                std::cout << "Min unbalanced force at node " << dof_iterator->Id()
                          << ", dof " << imin << " (" << dof_iterator->GetVariable().Name() << "), value = " << min_force
                          << std::endl;
            }
        }
    }

}; // class DofUtility

} //namespace Kratos.

#endif /* KRATOS_DOF_UTILITY_H_INCLUDED defined */

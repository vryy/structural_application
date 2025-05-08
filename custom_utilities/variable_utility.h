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
*   Date:                $Date: 24 May 2018 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

#if !defined(KRATOS_VARIABLE_UTILITY_INCLUDED )
#define  KRATOS_VARIABLE_UTILITY_INCLUDED

//System includes

//External includes

//Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h"

namespace Kratos
{

/**
 * Abstract class for all utility for all variable transfer utilities
 */
class VariableUtility
{
public:

    KRATOS_CLASS_POINTER_DEFINITION( VariableUtility );

    typedef ModelPart::IndexType IndexType;
    typedef ModelPart::NodesContainerType NodesContainerType;
    typedef ModelPart::ElementsContainerType ElementsContainerType;
    typedef ModelPart::DofsArrayType DofsArrayType;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> DenseSpaceType;
    typedef LinearSolver<SparseSpaceType, DenseSpaceType, ModelPart> LinearSolverType;

    struct DoubleVariableInitializer
    {
        typedef Variable<double> VariableType;
        typedef double DataType;

        static void inline Initialize(std::vector<double>& rValues)
        {
            for (std::size_t i = 0; i < rValues.size(); ++i)
                rValues[i] = 0.0;
        }

        static void inline Initialize(double& rValue)
        {
            rValue = 0.0;
        }

        static void inline Initialize(double& rA, const double rB)
        {
            rA = rB;
        }
    };

    struct Array1DVariableInitializer
    {
        typedef Variable<array_1d<double, 3> > VariableType;
        typedef array_1d<double, 3> DataType;

        static void inline Initialize(std::vector<array_1d<double, 3> >& rValues)
        {
            for (std::size_t i = 0; i < rValues.size(); ++i)
                noalias(rValues[i]) = ZeroVector(3);
        }

        static void inline Initialize(array_1d<double, 3>& rValue)
        {
            noalias(rValue) = ZeroVector(3);
        }

        static void inline Initialize(array_1d<double, 3>& rA, const array_1d<double, 3>& rB)
        {
            noalias(rA) = rB;
        }
    };

    template<int TSize>
    struct VectorVariableInitializer
    {
        typedef Variable<Vector> VariableType;
        typedef Vector DataType;

        static void inline Initialize(std::vector<Vector>& rValues)
        {
            for (std::size_t i = 0; i < rValues.size(); ++i)
            {
                rValues[i].resize(TSize, false);
                noalias(rValues[i]) = ZeroVector(TSize);
            }
        }

        static void inline Initialize(Vector& rValue)
        {
            rValue.resize(TSize, false);
            noalias(rValue) = ZeroVector(TSize);
        }

        static void inline Initialize(Vector& rA, const Vector& rB)
        {
            noalias(rA) = rB;
        }
    };

    /**
     * Constructor.
     */
    VariableUtility()
    : mEchoLevel(0)
    {
        std::cout << "VariableUtility created" << std::endl;
    }

    VariableUtility(const ElementsContainerType& pElements)
    : mEchoLevel(0), mpElements(pElements)
    {
        std::cout << "VariableUtility created, number of elements = " << mpElements.size() << std::endl;
    }

    VariableUtility(const ElementsContainerType& pElements, const int EchoLevel)
    : mEchoLevel(EchoLevel), mpElements(pElements)
    {
        std::cout << "VariableUtility created, number of elements = " << mpElements.size() << std::endl;
    }

    /// Extract the solution vector from list of dofs
    void GetSolutionVector(SparseSpaceType::VectorType& X, const DofsArrayType& rDofSet,
        const IndexType EquationSystemSize, const IndexType SolutionStepIndex = 0) const
    {
        if (SparseSpaceType::Size(X) != EquationSystemSize)
            SparseSpaceType::Resize(X, EquationSystemSize);
        SparseSpaceType::SetToZero(X);

        for (DofsArrayType::const_iterator i_dof = rDofSet.begin(); i_dof != rDofSet.end(); ++i_dof)
        {
            if (i_dof->EquationId() < EquationSystemSize)
            {
                SparseSpaceType::SetValue(X, i_dof->EquationId(), i_dof->GetSolutionStepValue(SolutionStepIndex));
            }
        }
    }

    /// Extract the solution vector from double variable
    /// The output vector is contiguous. The node id is arranged in the ascending order.
    void GetSolutionVector(SparseSpaceType::VectorType& X, const Variable<double>& rVariable,
        const ModelPart& r_model_part, const IndexType SolutionStepIndex = 0) const
    {
        // collect the number of nodes and assign unique row number
        const NodesContainerType& nodes = r_model_part.Nodes();
        std::map<IndexType, IndexType> row_indices;
        IndexType cnt = 0;
        for (NodesContainerType::const_iterator i_node = nodes.begin(); i_node != nodes.end(); ++i_node)
            row_indices[i_node->Id()] = cnt++;

        // extract the solution vector
        if (SparseSpaceType::Size(X) != cnt)
            SparseSpaceType::Resize(X, cnt);
        SparseSpaceType::SetToZero(X);

        std::size_t row;
        for (NodesContainerType::const_iterator i_node = nodes.begin(); i_node != nodes.end(); ++i_node)
        {
            row = row_indices[i_node->Id()];
            SparseSpaceType::SetValue(X, row, i_node->GetSolutionStepValue(rVariable, SolutionStepIndex));
        }
    }

    /// Extract the solution vector from double variable
    /// The output vector is contiguous. The node id is arranged in the ascending order.
    void GetSolutionVector(SparseSpaceType::VectorType& X, const Variable<array_1d<double, 3> >& rVariable,
        const ModelPart& r_model_part, const IndexType SolutionStepIndex = 0) const
    {
        // collect the number of nodes and assign unique row number
        const NodesContainerType& nodes = r_model_part.Nodes();
        std::map<IndexType, IndexType> row_indices;
        IndexType cnt = 0;
        for (NodesContainerType::const_iterator i_node = nodes.begin(); i_node != nodes.end(); ++i_node)
            row_indices[i_node->Id()] = cnt++;

        // extract the solution vector
        if (SparseSpaceType::Size(X) != 3*cnt)
            SparseSpaceType::Resize(X, 3*cnt);
        SparseSpaceType::SetToZero(X);

        std::size_t row;
        for (NodesContainerType::const_iterator i_node = nodes.begin(); i_node != nodes.end(); ++i_node)
        {
            row = row_indices[i_node->Id()];
            const array_1d<double, 3>& values = i_node->GetSolutionStepValue(rVariable, SolutionStepIndex);
            for (int i = 0; i < 3; ++i)
                SparseSpaceType::SetValue(X, 3*row+i, values[i]);
        }
    }

    /**
     * Destructor.
     */
    virtual ~VariableUtility()
    {}

    /**
     * Access
     */

    void SetEchoLevel(const int Level)
    {
        mEchoLevel = Level;
    }

    int GetEchoLevel() const
    {
        return mEchoLevel;
    }

protected:

    ElementsContainerType mpElements;

    /// Initialize the utilitey
    virtual void Initialize( const ElementsContainerType& pElements )
    {
        KRATOS_ERROR << "Error calling base class function";
    }

private:

    int mEchoLevel;

}; // Class VariableUtility

} // namespace Kratos.

#endif /* KRATOS_VARIABLE_UTILITY_INCLUDED  defined */

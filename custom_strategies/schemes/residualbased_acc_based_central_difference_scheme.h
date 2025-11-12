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
*   Date:                $Date: 18 Feb 2023 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/


#if !defined(KRATOS_STRUCTURAL_APPLICATION_RESIDUALBASED_ACC_BASED_CENTRAL_DIFFERENCE_SCHEME )
#define  KRATOS_STRUCTURAL_APPLICATION_RESIDUALBASED_ACC_BASED_CENTRAL_DIFFERENCE_SCHEME

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "includes/element.h"
#ifdef SD_APP_FORWARD_COMPATIBILITY
#include "custom_python3/legacy_structural_app_vars.h"
#else
#include "includes/legacy_structural_app_vars.h"
#endif
#include "containers/array_1d.h"
#include "solving_strategies/schemes/scheme.h"
#include "structural_application_variables.h"
#include "custom_elements/prescribed_object.h"

namespace Kratos
{
/**@name Kratos Globals */
/*@{ */
/*@} */
/**@name Type Definitions */
/*@{ */
/*@} */
/**@name  Enum's */
/*@{ */
/*@} */
/**@name  Functions */
/*@{ */
/*@} */
/**@name Kratos Classes */
/*@{ */
/**
Implementation of time integration scheme with acceleration as primal variables

In this scheme, the displacement and velocity are approximated as following:

u_{n+1} = u_{n} + dt * ud_{n} + 1/2 * dt^2 * udd_{n}
ud_{n+1/2} = ud_{n} + 1/2 * dt * udd_{n}

The structural dynamics equation to be solved:
    M*udd + D*ud + r_i(u) = r_e                     [1]

Temporal discretization of [1} is:
M*udd_{n+1} = r_ext - r_i(u_{n+1}) - D*ud_{n+1/2}     [2]

[2] can be solved trivially if M is lumped mass matrix.

The end-step velocity is updated as:

ud_{n+1} = ud_{n+1/2} + 1/2 * dt * udd_{n+1}

This scheme is second order accurate.
*/
template<class TSparseSpace,  class TDenseSpace>
class ResidualBasedAccBasedCentralDifferenceScheme: public Scheme<TSparseSpace,TDenseSpace>
{
public:
    /**@name Type Definitions */
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedAccBasedCentralDifferenceScheme );

    typedef Scheme<TSparseSpace,TDenseSpace> BaseType;

    typedef TSparseSpace SparseSpaceType;
    typedef TDenseSpace DenseSpaceType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::ElementsArrayType ElementsArrayType;

    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

    typedef typename Element::DofsVectorType DofsVectorType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    /**
     * Constructor for Central Difference Scheme
     * @ref
     */
    ResidualBasedAccBasedCentralDifferenceScheme() : BaseType()
    {
        // set all the integration flag to false. To incorporate the corresponding term, the flag shall be set by user.
        mIntegrateRotation = false;
        mIntegrateMultiplier = false;
        mIntegrateLoad = false;
        mForceLumpedMass = false;

        std::cout << "Using the Acceleration-based Central Difference scheme for structural dynamics " << std::endl;
    }

    /**
     * Constructor for Central Difference Scheme
     * @ref
     */
    ResidualBasedAccBasedCentralDifferenceScheme(const bool forceLumpedMass) : BaseType()
    {
        // set all the integration flag to false. To incorporate the corresponding term, the flag shall be set by user.
        mIntegrateRotation = false;
        mIntegrateMultiplier = false;
        mIntegrateLoad = false;
        mForceLumpedMass = forceLumpedMass;

        if (!mForceLumpedMass)
            std::cout << "Using the Acceleration-based Central Difference scheme for structural dynamics" << std::endl;
        else
            std::cout << "Using the Acceleration-based Central Difference scheme for structural dynamics with lumped mass matrix" << std::endl;
    }

    /** Destructor.*/
    virtual ~ResidualBasedAccBasedCentralDifferenceScheme()
    {}

    /**@name Operators */

    /// Copy assignment
    ResidualBasedAccBasedCentralDifferenceScheme& operator=(const ResidualBasedAccBasedCentralDifferenceScheme& rOther)
    {
        this->mIntegrateRotation = rOther.mIntegrateRotation;
        this->mIntegrateMultiplier = mIntegrateMultiplier;
        this->mIntegrateLoad = mIntegrateLoad;
        this->mForceLumpedMass = mForceLumpedMass;
        return *this;
    }

    /*@} */
    /**@name Operations */
    /*@{ */

    /// Enable integration of rotation d.o.f
    void SetIntegrateRotation(const bool& value)
    {
        mIntegrateRotation = value;
    }

    /// Enable integration of multiplier d.o.f
    void SetIntegrateMultiplier(const bool& value)
    {
        mIntegrateMultiplier = value;
    }

    /// Enable time values of load
    void SetIntegrateLoad(const bool& value)
    {
        mIntegrateLoad = value;
    }

    /** Performing the update of the solution.*/
    /**
     * incremental update within newton iteration. It updates the state variables at the end of the time step: u_{n+1}^{k+1}= u_{n+1}^{k}+ \Delta u
     * @param r_model_part
     * @param rDofSet set of all primary variables
     * @param A LHS matrix
     * @param Dx incremental update of primary variables
     * @param b RHS Vector
     */
    void UpdateAll(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b )
    {
        KRATOS_TRY

        for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ;
                i != r_model_part.NodesEnd() ; i++)
        {
            if( i->HasDofFor(DISPLACEMENT_X) )
            {
                if(i->GetDof(DISPLACEMENT_X).IsFree())
                {
                    i->GetSolutionStepValue(ACCELERATION_EINS_X)
                    =Dx[i->GetDof(DISPLACEMENT_X).EquationId()];
                }
            }
            if( i->HasDofFor(DISPLACEMENT_Y) )
            {
                if(i->GetDof(DISPLACEMENT_Y).IsFree())
                {
                    i->GetSolutionStepValue(ACCELERATION_EINS_Y)
                    =Dx[i->GetDof(DISPLACEMENT_Y).EquationId()];
                }
            }
            if( i->HasDofFor(DISPLACEMENT_Z) )
            {
                if(i->GetDof(DISPLACEMENT_Z).IsFree())
                {
                    i->GetSolutionStepValue(ACCELERATION_EINS_Z)
                    =Dx[i->GetDof(DISPLACEMENT_Z).EquationId()];
                }
            }
            if( i->HasDofFor(WATER_PRESSURE) )
            {
                if(i->GetDof(WATER_PRESSURE).IsFree())
                {
                    i->GetSolutionStepValue(WATER_PRESSURE_EINS_ACCELERATION)
                    =Dx[i->GetDof(WATER_PRESSURE).EquationId()];
                }
            }
            if( i->HasDofFor(AIR_PRESSURE) )
            {
                if(i->GetDof(AIR_PRESSURE).IsFree())
                {
                    i->GetSolutionStepValue(AIR_PRESSURE_EINS_ACCELERATION)
                    =Dx[i->GetDof(AIR_PRESSURE).EquationId()];
                }
            }

            // // update this variables to account for Lagrange multiplier
            // TODO
            // if( i->HasDofFor(LAGRANGE_DISPLACEMENT_X) )
            // {
            //     if(i->GetDof(LAGRANGE_DISPLACEMENT_X).IsFree())
            //     {
            //         i->GetSolutionStepValue(LAGRANGE_DISPLACEMENT_X_A)
            //         +=Dx[i->GetDof(LAGRANGE_DISPLACEMENT_X).EquationId()];
            //     }
            // }
            // if( i->HasDofFor(LAGRANGE_DISPLACEMENT_Y) )
            // {
            //     if(i->GetDof(LAGRANGE_DISPLACEMENT_Y).IsFree())
            //     {
            //         i->GetSolutionStepValue(LAGRANGE_DISPLACEMENT_Y)
            //         +=Dx[i->GetDof(LAGRANGE_DISPLACEMENT_Y).EquationId()];
            //     }
            // }
            // if( i->HasDofFor(LAGRANGE_DISPLACEMENT_Z) )
            // {
            //     if(i->GetDof(LAGRANGE_DISPLACEMENT_Z).IsFree())
            //     {
            //         i->GetSolutionStepValue(LAGRANGE_DISPLACEMENT_Z)
            //         +=Dx[i->GetDof(LAGRANGE_DISPLACEMENT_Z).EquationId()];
            //     }
            // }
            // if( i->HasDofFor(LAGRANGE_WATER_PRESSURE) )
            // {
            //     if(i->GetDof(LAGRANGE_WATER_PRESSURE).IsFree())
            //     {
            //         i->GetSolutionStepValue(LAGRANGE_WATER_PRESSURE)
            //         +=Dx[i->GetDof(LAGRANGE_WATER_PRESSURE).EquationId()];
            //     }
            // }
            // if( i->HasDofFor(LAGRANGE_AIR_PRESSURE) )
            // {
            //     if(i->GetDof(LAGRANGE_AIR_PRESSURE).IsFree())
            //     {
            //         i->GetSolutionStepValue(LAGRANGE_AIR_PRESSURE)
            //         +=Dx[i->GetDof(LAGRANGE_AIR_PRESSURE).EquationId()];
            //     }
            // }
        }

        KRATOS_CATCH("")
    }

    /** Performing the update of the solution.*/
    /**
     * incremental update within newton iteration. It updates the state variables at the end of the time step: u_{n+1}^{k+1}= u_{n+1}^{k}+ \Delta u
     * @param r_model_part
     * @param rDofSet set of all primary variables
     * @param A LHS matrix
     * @param Dx incremental update of primary variables
     * @param b RHS Vector
     */
    void Update(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b ) final
    {
        KRATOS_TRY

        for (typename DofsArrayType::iterator dof_iterator = rDofSet.begin(); dof_iterator != rDofSet.end(); ++dof_iterator)
        {
            ModelPart::NodeType& rNode = r_model_part.GetNode(dof_iterator->Id());

            if (dof_iterator->GetVariable() == DISPLACEMENT_X)
            {
                if (dof_iterator->IsFree())
                {
                    rNode.GetSolutionStepValue(ACCELERATION_EINS_X)
                    = Dx[dof_iterator->EquationId()];
                }
            }
            else if (dof_iterator->GetVariable() == DISPLACEMENT_Y)
            {
                if (dof_iterator->IsFree())
                {
                    rNode.GetSolutionStepValue(ACCELERATION_EINS_Y)
                    = Dx[dof_iterator->EquationId()];
                }
            }
            else if (dof_iterator->GetVariable() == DISPLACEMENT_Z)
            {
                if (dof_iterator->IsFree())
                {
                    rNode.GetSolutionStepValue(ACCELERATION_EINS_Z)
                    = Dx[dof_iterator->EquationId()];
                }
            }
            else if (dof_iterator->GetVariable() == WATER_PRESSURE)
            {
                if (dof_iterator->IsFree())
                {
                    rNode.GetSolutionStepValue(WATER_PRESSURE_EINS_ACCELERATION)
                    = Dx[dof_iterator->EquationId()];
                }
            }
            else if (dof_iterator->GetVariable() == AIR_PRESSURE)
            {
                if (dof_iterator->IsFree())
                {
                    rNode.GetSolutionStepValue(AIR_PRESSURE_EINS_ACCELERATION)
                    = Dx[dof_iterator->EquationId()];
                }
            }
            else
            {
                if (mIntegrateRotation)
                {
                    if (dof_iterator->GetVariable() == ROTATION_X)
                    {
                        if (dof_iterator->IsFree())
                        {
                            rNode.GetSolutionStepValue(ANGULAR_ACCELERATION_EINS_X)
                            = Dx[dof_iterator->EquationId()];
                        }
                    }
                    else if (dof_iterator->GetVariable() == ROTATION_Y)
                    {
                        if (dof_iterator->IsFree())
                        {
                            rNode.GetSolutionStepValue(ANGULAR_ACCELERATION_EINS_Y)
                            = Dx[dof_iterator->EquationId()];
                        }
                    }
                    else if (dof_iterator->GetVariable() == ROTATION_Z)
                    {
                        if (dof_iterator->IsFree())
                        {
                            rNode.GetSolutionStepValue(ANGULAR_ACCELERATION_EINS_Z)
                            = Dx[dof_iterator->EquationId()];
                        }
                    }
                }

                if (mIntegrateMultiplier)
                {
                    if (dof_iterator->GetVariable() == LAMBDA)
                    {
                        if (dof_iterator->IsFree())
                        {
                            rNode.GetSolutionStepValue(LAMBDA_EINS_DT2)
                            = Dx[dof_iterator->EquationId()];
                        }
                    }
                }

                if (!mIntegrateMultiplier && !mIntegrateRotation)
                {
                    // update for non-dynamics variable
                    if (dof_iterator->IsFree())
                    {
                        dof_iterator->GetSolutionStepValue()
                        = Dx[dof_iterator->EquationId()];
                    }
                }
            }
        }

        KRATOS_CATCH("")
    }

    /**
     * initializes next iteration step. In the explicit scheme, this is called once.
     * @param r_model_part
     * @param A LHS matrix
     * @param Dx incremental update of primary variables
     * @param b RHS Vector
     */
    void InitializeNonLinIteration(ModelPart& r_model_part,
                                   TSystemMatrixType& A,
                                   TSystemVectorType& Dx,
                                   TSystemVectorType& b ) final
    {
        KRATOS_TRY

        const ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
        const double Dt = CurrentProcessInfo[DELTA_TIME];

        ElementsArrayType& pElements = r_model_part.Elements();
        bool element_is_active;
        for (typename ElementsArrayType::iterator it = pElements.begin(); it != pElements.end(); ++it)
        {
            element_is_active = true;
            if(it->IsDefined(ACTIVE))
                element_is_active = it->Is(ACTIVE);
            if (it->Has(IS_INACTIVE))
                element_is_active |= (!it->GetValue(IS_INACTIVE));

            if (element_is_active)
                it->InitializeNonLinearIteration(CurrentProcessInfo);
        }

        bool condition_is_active;
        ConditionsArrayType& pConditions = r_model_part.Conditions();
        for (typename ConditionsArrayType::iterator it = pConditions.begin(); it != pConditions.end(); ++it)
        {
            condition_is_active = true;
            if( it->IsDefined( ACTIVE ) )
                condition_is_active = it->Is(ACTIVE);
            if (it->Has(IS_INACTIVE))
                condition_is_active |= (!it->GetValue(IS_INACTIVE));

            if (condition_is_active)
                it->InitializeNonLinearIteration(CurrentProcessInfo);
        }

        // Update nodal values and nodal velocities
        for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; i != r_model_part.NodesEnd() ; i++)
        {
            if( i->HasDofFor(DISPLACEMENT_X) )
            {
                i->GetSolutionStepValue(DISPLACEMENT_EINS_X) = i->GetSolutionStepValue(DISPLACEMENT_NULL_X)
                    + Dt * i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_X)
                    + 0.5 * Dt * Dt * i->GetSolutionStepValue(ACCELERATION_NULL_X);

                i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_X) = i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_X)
                    + 0.5 * Dt * i->GetSolutionStepValue(ACCELERATION_NULL_X);

                i->GetSolutionStepValue(DISPLACEMENT_X) = i->GetSolutionStepValue(DISPLACEMENT_EINS_X);
                i->GetSolutionStepValue(DISPLACEMENT_DT_X) = i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_X);
            }
            if( i->HasDofFor(DISPLACEMENT_Y) )
            {
                i->GetSolutionStepValue(DISPLACEMENT_EINS_Y) = i->GetSolutionStepValue(DISPLACEMENT_NULL_Y)
                    + Dt * i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_Y)
                    + 0.5 * Dt * Dt * i->GetSolutionStepValue(ACCELERATION_NULL_Y);

                i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_Y) = i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_Y)
                    + 0.5 * Dt * i->GetSolutionStepValue(ACCELERATION_NULL_Y);

                i->GetSolutionStepValue(DISPLACEMENT_Y) = i->GetSolutionStepValue(DISPLACEMENT_EINS_Y);
                i->GetSolutionStepValue(DISPLACEMENT_DT_Y) = i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_Y);
            }
            if( i->HasDofFor(DISPLACEMENT_Z) )
            {
                i->GetSolutionStepValue(DISPLACEMENT_EINS_Z) = i->GetSolutionStepValue(DISPLACEMENT_NULL_Z)
                    + Dt * i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_Z)
                    + 0.5 * Dt * Dt * i->GetSolutionStepValue(ACCELERATION_NULL_Z);

                i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_Z) = i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_Z)
                    + 0.5 * Dt * i->GetSolutionStepValue(ACCELERATION_NULL_Z);

                i->GetSolutionStepValue(DISPLACEMENT_Z) = i->GetSolutionStepValue(DISPLACEMENT_EINS_Z);
                i->GetSolutionStepValue(DISPLACEMENT_DT_Z) = i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_Z);
            }
            if( i->HasDofFor(WATER_PRESSURE) )
            {
                i->GetSolutionStepValue(WATER_PRESSURE_EINS) = i->GetSolutionStepValue(WATER_PRESSURE_NULL)
                    + Dt * i->GetSolutionStepValue(WATER_PRESSURE_NULL_DT)
                    + 0.5 * Dt * Dt * i->GetSolutionStepValue(WATER_PRESSURE_NULL_ACCELERATION);

                i->GetSolutionStepValue(WATER_PRESSURE_EINS_DT) = i->GetSolutionStepValue(WATER_PRESSURE_NULL_DT)
                    + 0.5 * Dt * i->GetSolutionStepValue(WATER_PRESSURE_NULL_ACCELERATION);

                i->GetSolutionStepValue(WATER_PRESSURE) = i->GetSolutionStepValue(WATER_PRESSURE_EINS);
                i->GetSolutionStepValue(WATER_PRESSURE_DT) = i->GetSolutionStepValue(WATER_PRESSURE_EINS_DT);
            }
            if( i->HasDofFor(AIR_PRESSURE) )
            {
                i->GetSolutionStepValue(AIR_PRESSURE_EINS) = i->GetSolutionStepValue(AIR_PRESSURE_NULL)
                    + Dt * i->GetSolutionStepValue(AIR_PRESSURE_NULL_DT)
                    + 0.5 * Dt * Dt * i->GetSolutionStepValue(AIR_PRESSURE_NULL_ACCELERATION);

                i->GetSolutionStepValue(AIR_PRESSURE_EINS_DT) = i->GetSolutionStepValue(AIR_PRESSURE_NULL_DT)
                    + 0.5 * Dt * i->GetSolutionStepValue(AIR_PRESSURE_NULL_ACCELERATION);

                i->GetSolutionStepValue(AIR_PRESSURE) = i->GetSolutionStepValue(AIR_PRESSURE_EINS);
                i->GetSolutionStepValue(AIR_PRESSURE_DT) = i->GetSolutionStepValue(AIR_PRESSURE_EINS_DT);
            }

            if (mIntegrateRotation)
            {
                if( i->HasDofFor(ROTATION_X) )
                {
                    i->GetSolutionStepValue(ROTATION_EINS_X) = i->GetSolutionStepValue(ROTATION_NULL_X)
                        + Dt * i->GetSolutionStepValue(ROTATION_NULL_DT_X)
                        + 0.5 * Dt * Dt * i->GetSolutionStepValue(ANGULAR_ACCELERATION_NULL_X);

                    i->GetSolutionStepValue(ROTATION_EINS_DT_X) = i->GetSolutionStepValue(ROTATION_NULL_DT_X)
                        + 0.5 * Dt * i->GetSolutionStepValue(ANGULAR_ACCELERATION_NULL_X);

                    i->GetSolutionStepValue(ROTATION_X) = i->GetSolutionStepValue(ROTATION_EINS_X);
                    i->GetSolutionStepValue(ROTATION_DT_X) = i->GetSolutionStepValue(ROTATION_EINS_DT_X);
                }
                if( i->HasDofFor(ROTATION_Y) )
                {
                    i->GetSolutionStepValue(ROTATION_EINS_Y) = i->GetSolutionStepValue(ROTATION_NULL_Y)
                        + Dt * i->GetSolutionStepValue(ROTATION_NULL_DT_Y)
                        + 0.5 * Dt * Dt * i->GetSolutionStepValue(ANGULAR_ACCELERATION_NULL_Y);

                    i->GetSolutionStepValue(ROTATION_EINS_DT_Y) = i->GetSolutionStepValue(ROTATION_NULL_DT_Y)
                        + 0.5 * Dt * i->GetSolutionStepValue(ANGULAR_ACCELERATION_NULL_Y);

                    i->GetSolutionStepValue(ROTATION_Y) = i->GetSolutionStepValue(ROTATION_EINS_Y);
                    i->GetSolutionStepValue(ROTATION_DT_Y) = i->GetSolutionStepValue(ROTATION_EINS_DT_Y);
                }
                if( i->HasDofFor(ROTATION_Z) )
                {
                    i->GetSolutionStepValue(ROTATION_EINS_Z) = i->GetSolutionStepValue(ROTATION_NULL_Z)
                        + Dt * i->GetSolutionStepValue(ROTATION_NULL_DT_Z)
                        + 0.5 * Dt * Dt * i->GetSolutionStepValue(ANGULAR_ACCELERATION_NULL_Z);

                    i->GetSolutionStepValue(ROTATION_EINS_DT_Z) = i->GetSolutionStepValue(ROTATION_NULL_DT_Z)
                        + 0.5 * Dt * i->GetSolutionStepValue(ANGULAR_ACCELERATION_NULL_Z);

                    i->GetSolutionStepValue(ROTATION_Z) = i->GetSolutionStepValue(ROTATION_EINS_Z);
                    i->GetSolutionStepValue(ROTATION_DT_Z) = i->GetSolutionStepValue(ROTATION_EINS_DT_Z);
                }
            }

            if (mIntegrateMultiplier)
            {
                if( i->HasDofFor(LAMBDA) )
                {
                    i->GetSolutionStepValue(LAMBDA_EINS) = i->GetSolutionStepValue(LAMBDA_NULL)
                        + Dt * i->GetSolutionStepValue(LAMBDA_NULL_DT)
                        + 0.5 * Dt * Dt * i->GetSolutionStepValue(LAMBDA_NULL_DT2);

                    i->GetSolutionStepValue(LAMBDA_EINS_DT) = i->GetSolutionStepValue(LAMBDA_NULL_DT)
                        + 0.5 * Dt * i->GetSolutionStepValue(LAMBDA_NULL_DT2);

                    i->GetSolutionStepValue(LAMBDA) = i->GetSolutionStepValue(LAMBDA_EINS);
                    i->GetSolutionStepValue(LAMBDA_DT) = i->GetSolutionStepValue(LAMBDA_EINS_DT);
                }
            }

            if (mIntegrateLoad)
            {
                if( i->SolutionStepsDataHas(FACE_LOAD) )
                {
                    noalias( i->GetSolutionStepValue(FACE_LOAD) )
                    = i->GetSolutionStepValue(FACE_LOAD_NULL);
                }

                if( i->SolutionStepsDataHas(FORCE) )
                {
                    noalias( i->GetSolutionStepValue(FORCE) )
                    = i->GetSolutionStepValue(FORCE_NULL);
                }
            }
        }

        //For total Lagrangian
        for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ;
                i != r_model_part.NodesEnd() ; ++i)
        {
            (i)->X() = (i)->X0() + i->GetSolutionStepValue(DISPLACEMENT_X);
            (i)->Y() = (i)->Y0() + i->GetSolutionStepValue(DISPLACEMENT_Y);
            (i)->Z() = (i)->Z0() + i->GetSolutionStepValue(DISPLACEMENT_Z);
        }

        KRATOS_CATCH("")
    }

    /**
     */
    void FinalizeNonLinIteration(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) final
    {
        KRATOS_TRY

        // to account for prescribed displacement, the displacement at prescribed nodes need to be updated
        double curr_disp, delta_disp;
        for (ModelPart::NodesContainerType::iterator it_node = r_model_part.Nodes().begin(); it_node != r_model_part.Nodes().end(); ++it_node)
        {
            if (it_node->IsFixed(DISPLACEMENT_X))
            {
                curr_disp = it_node->GetSolutionStepValue(DISPLACEMENT_X);
                delta_disp = it_node->GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X);
                it_node->GetSolutionStepValue(DISPLACEMENT_X) = curr_disp + delta_disp;
                it_node->GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X) = 0.0; // set the prescribed displacement to zero to avoid update in the second step
            }

            if (it_node->IsFixed(DISPLACEMENT_Y))
            {
                curr_disp = it_node->GetSolutionStepValue(DISPLACEMENT_Y);
                delta_disp = it_node->GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y);
                it_node->GetSolutionStepValue(DISPLACEMENT_Y) = curr_disp + delta_disp;
                it_node->GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y) = 0.0; // set the prescribed displacement to zero to avoid update in the second step
            }

            if (it_node->IsFixed(DISPLACEMENT_Z))
            {
                curr_disp = it_node->GetSolutionStepValue(DISPLACEMENT_Z);
                delta_disp = it_node->GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Z);
                it_node->GetSolutionStepValue(DISPLACEMENT_Z) = curr_disp + delta_disp;
                it_node->GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Z) = 0.0; // set the prescribed displacement to zero to avoid update in the second step
            }
        }

        // invoking the element and condition finalization after an iteration

        const ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
        const double Dt = CurrentProcessInfo[DELTA_TIME];

        ElementsArrayType& pElements = r_model_part.Elements();
        bool element_is_active;
        for (typename ElementsArrayType::iterator it = pElements.begin(); it != pElements.end(); ++it)
        {
            element_is_active = true;
            if(it->IsDefined(ACTIVE))
                element_is_active = it->Is(ACTIVE);
            if (it->Has(IS_INACTIVE))
                element_is_active |= (!it->GetValue(IS_INACTIVE));

            if (element_is_active)
                it->FinalizeNonLinearIteration(CurrentProcessInfo);
        }

        bool condition_is_active;
        ConditionsArrayType& pConditions = r_model_part.Conditions();
        for (typename ConditionsArrayType::iterator it = pConditions.begin(); it != pConditions.end(); ++it)
        {
            condition_is_active = true;
            if( it->IsDefined( ACTIVE ) )
                condition_is_active = it->Is(ACTIVE);
            if (it->Has(IS_INACTIVE))
                condition_is_active |= (!it->GetValue(IS_INACTIVE));

            if (condition_is_active)
                it->FinalizeNonLinearIteration(CurrentProcessInfo);
        }

        // Update nodal values and nodal velocities
        for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; i != r_model_part.NodesEnd() ; i++)
        {
            if( i->HasDofFor(DISPLACEMENT_X) )
            {
                i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_X) += 0.5 * Dt * i->GetSolutionStepValue(ACCELERATION_EINS_X);
            }
            if( i->HasDofFor(DISPLACEMENT_Y) )
            {
                i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_Y) += 0.5 * Dt * i->GetSolutionStepValue(ACCELERATION_EINS_Y);
            }
            if( i->HasDofFor(DISPLACEMENT_Z) )
            {
                i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_Z) += 0.5 * Dt * i->GetSolutionStepValue(ACCELERATION_EINS_Z);
            }
            if( i->HasDofFor(WATER_PRESSURE) )
            {
                i->GetSolutionStepValue(WATER_PRESSURE_EINS_DT) += 0.5 * Dt * i->GetSolutionStepValue(WATER_PRESSURE_EINS_ACCELERATION);
            }
            if( i->HasDofFor(AIR_PRESSURE) )
            {
                i->GetSolutionStepValue(AIR_PRESSURE_EINS_DT) += 0.5 * Dt * i->GetSolutionStepValue(AIR_PRESSURE_EINS_ACCELERATION);
            }

            if (mIntegrateRotation)
            {
                if( i->HasDofFor(ROTATION_X) )
                {
                    i->GetSolutionStepValue(ROTATION_EINS_DT_X) += 0.5 * Dt * i->GetSolutionStepValue(ANGULAR_ACCELERATION_EINS_X);
                }
                if( i->HasDofFor(ROTATION_Y) )
                {
                    i->GetSolutionStepValue(ROTATION_EINS_DT_Y) += 0.5 * Dt * i->GetSolutionStepValue(ANGULAR_ACCELERATION_EINS_Y);
                }
                if( i->HasDofFor(ROTATION_Z) )
                {
                    i->GetSolutionStepValue(ROTATION_EINS_DT_Z) += 0.5 * Dt * i->GetSolutionStepValue(ANGULAR_ACCELERATION_EINS_Z);
                }
            }

            if (mIntegrateMultiplier)
            {
                if( i->HasDofFor(LAMBDA) )
                {
                    i->GetSolutionStepValue(LAMBDA_EINS_DT) += 0.5 * Dt * i->GetSolutionStepValue(LAMBDA_EINS_DT2);
                }
            }
        }

        KRATOS_CATCH("")
    }

    /**
     * initializes time step solution
     * only for reasons if the time step solution is restarted
     * @param r_model_part
     * @param A LHS matrix
     * @param Dx incremental update of primary variables
     * @param b RHS Vector
     */
    void InitializeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) final
    {
        KRATOS_TRY

        const ProcessInfo& CurrentProcessInfo= r_model_part.GetProcessInfo();

        BaseType::InitializeSolutionStep(r_model_part,A,Dx,b);
        //Update nodal values and nodal velocities
        for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ;
                i != r_model_part.NodesEnd() ; i++)
        {
            if( i->HasDofFor(DISPLACEMENT_X) &&  i->GetDof(DISPLACEMENT_X).IsFree())
            {
                i->GetSolutionStepValue(ACCELERATION_EINS_X)=
                    i->GetSolutionStepValue(ACCELERATION_NULL_X);
                i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_X)=
                    i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_X);
                i->GetSolutionStepValue(DISPLACEMENT_EINS_X )=
                    i->GetSolutionStepValue(DISPLACEMENT_NULL_X);
            }
            if( i->HasDofFor(DISPLACEMENT_Y) &&  i->GetDof(DISPLACEMENT_Y).IsFree())
            {
                i->GetSolutionStepValue(ACCELERATION_EINS_Y)=
                    i->GetSolutionStepValue(ACCELERATION_NULL_Y);
                i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_Y)=
                    i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_Y);
                i->GetSolutionStepValue(DISPLACEMENT_EINS_Y)=
                    i->GetSolutionStepValue(DISPLACEMENT_NULL_Y);
            }
            if( i->HasDofFor(DISPLACEMENT_Z) && i->GetDof(DISPLACEMENT_Z).IsFree() )
            {
                i->GetSolutionStepValue(ACCELERATION_EINS_Y)=
                    i->GetSolutionStepValue(ACCELERATION_NULL_Z);
                i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_Z)=
                    i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_Z);
                i->GetSolutionStepValue(DISPLACEMENT_EINS_Z)=
                    i->GetSolutionStepValue(DISPLACEMENT_NULL_Z);
            }
            if( i->HasDofFor(WATER_PRESSURE) && i->GetDof(WATER_PRESSURE).IsFree())
            {
                i->GetSolutionStepValue(WATER_PRESSURE_EINS_DT)=
                    i->GetSolutionStepValue(WATER_PRESSURE_NULL_DT);
                i->GetSolutionStepValue(WATER_PRESSURE_EINS_ACCELERATION)=
                    i->GetSolutionStepValue(WATER_PRESSURE_NULL_ACCELERATION);
                i->GetSolutionStepValue(WATER_PRESSURE_EINS)=
                    i->GetSolutionStepValue(WATER_PRESSURE_NULL);
            }
            if( i->HasDofFor(AIR_PRESSURE) && i->GetDof(AIR_PRESSURE).IsFree())
            {
                i->GetSolutionStepValue(AIR_PRESSURE_EINS_DT)=
                    i->GetSolutionStepValue(AIR_PRESSURE_NULL_DT);
                i->GetSolutionStepValue(AIR_PRESSURE_EINS_ACCELERATION)=
                    i->GetSolutionStepValue(AIR_PRESSURE_NULL_ACCELERATION);
                i->GetSolutionStepValue(AIR_PRESSURE_EINS)=
                    i->GetSolutionStepValue(AIR_PRESSURE_NULL);
            }

            if (mIntegrateRotation)
            {
                if( i->HasDofFor(ROTATION_X) &&  i->GetDof(ROTATION_X).IsFree())
                {
                    i->GetSolutionStepValue(ANGULAR_ACCELERATION_EINS_X)=
                        i->GetSolutionStepValue(ANGULAR_ACCELERATION_NULL_X);
                    i->GetSolutionStepValue(ROTATION_EINS_DT_X)=
                        i->GetSolutionStepValue(ROTATION_NULL_DT_X);
                    i->GetSolutionStepValue(ROTATION_EINS_X)=
                        i->GetSolutionStepValue(ROTATION_NULL_X);
                }
                if( i->HasDofFor(ROTATION_Y) &&  i->GetDof(ROTATION_Y).IsFree())
                {
                    i->GetSolutionStepValue(ANGULAR_ACCELERATION_EINS_Y)=
                        i->GetSolutionStepValue(ANGULAR_ACCELERATION_NULL_Y);
                    i->GetSolutionStepValue(ROTATION_EINS_DT_Y)=
                        i->GetSolutionStepValue(ROTATION_NULL_DT_Y);
                    i->GetSolutionStepValue(ROTATION_EINS_Y)=
                        i->GetSolutionStepValue(ROTATION_NULL_Y);
                }
                if( i->HasDofFor(ROTATION_Z) &&  i->GetDof(ROTATION_X).IsFree())
                {
                    i->GetSolutionStepValue(ANGULAR_ACCELERATION_EINS_Z)=
                        i->GetSolutionStepValue(ANGULAR_ACCELERATION_NULL_Z);
                    i->GetSolutionStepValue(ROTATION_EINS_DT_Z)=
                        i->GetSolutionStepValue(ROTATION_NULL_DT_Z);
                    i->GetSolutionStepValue(ROTATION_EINS_Z)=
                        i->GetSolutionStepValue(ROTATION_NULL_Z);
                }
            }

            if (mIntegrateMultiplier)
            {
                if( i->HasDofFor(LAMBDA) &&  i->GetDof(LAMBDA).IsFree())
                {
                    i->GetSolutionStepValue(LAMBDA_EINS_DT2)=
                        i->GetSolutionStepValue(LAMBDA_NULL_DT2);
                    i->GetSolutionStepValue(LAMBDA_EINS_DT)=
                        i->GetSolutionStepValue(LAMBDA_NULL_DT);
                    i->GetSolutionStepValue(LAMBDA_EINS )=
                        i->GetSolutionStepValue(LAMBDA_NULL);
                }
            }

            // if (mIntegrateLoad)
            // {
            //     if( i->SolutionStepsDataHas(FACE_LOAD) )
            //     {
            //         noalias( i->GetSolutionStepValue(FACE_LOAD_EINS) )=
            //             i->GetSolutionStepValue(FACE_LOAD_NULL);
            //     }
            //     if( i->SolutionStepsDataHas(FORCE) )
            //     {
            //         noalias( i->GetSolutionStepValue(FORCE_EINS) )=
            //             i->GetSolutionStepValue(FORCE_NULL);
            //     }
            // }
        }

        KRATOS_CATCH("")
    }

    /**
     * finalizes time step solution
     * by setting u_n= u_n+1^k etc.
     * u_{n+1-alpha}= u_n*alpha+U_n+1^k+1*(1-alpha)
     * @param r_model_part
     * @param A LHS matrix
     * @param Dx incremental update of primary variables
     * @param b RHS Vector
     */
    void FinalizeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) final
    {
        const ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
        const double Dt = CurrentProcessInfo[DELTA_TIME];

        KRATOS_WATCH(CurrentProcessInfo[FIRST_TIME_STEP])

        // we manually FinalizeSolutionStep for each entities because the parent function is multithreaded
        ElementsArrayType& pElements = r_model_part.Elements();
        for (typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)
        {
            if ( (*it)->Is(ACTIVE) )
                (*it)->FinalizeSolutionStep(CurrentProcessInfo);
        }

        ConditionsArrayType& pConditions = r_model_part.Conditions();
        for (typename ConditionsArrayType::ptr_iterator it = pConditions.ptr_begin(); it != pConditions.ptr_end(); ++it)
        {
            if ( (*it)->Is(ACTIVE) )
                (*it)->FinalizeSolutionStep(CurrentProcessInfo);
        }

        if(CurrentProcessInfo[CALCULATE_INSITU_STRESS])
        {
            for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; i != r_model_part.NodesEnd() ; ++i)
            {
                if( i->HasDofFor(DISPLACEMENT_X))
                {
                    i->GetSolutionStepValue(DISPLACEMENT_OLD_X) = 0.0;
                    i->GetSolutionStepValue(DISPLACEMENT_NULL_X) = 0.0;
                    i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_X) = 0.0;
                    i->GetSolutionStepValue(ACCELERATION_NULL_X) = 0.0;
                }
                if( i->HasDofFor(DISPLACEMENT_Y) )
                {
                    i->GetSolutionStepValue(DISPLACEMENT_OLD_Y) = 0.0;
                    i->GetSolutionStepValue(DISPLACEMENT_NULL_Y) = 0.0;
                    i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_Y) = 0.0;
                    i->GetSolutionStepValue(ACCELERATION_NULL_Y) = 0.0;
                }
                if( i->HasDofFor(DISPLACEMENT_Z) )
                {
                    i->GetSolutionStepValue(DISPLACEMENT_OLD_Z) = 0.0;
                    i->GetSolutionStepValue(DISPLACEMENT_NULL_Z) = 0.0;
                    i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_Z) = 0.0;
                    i->GetSolutionStepValue(ACCELERATION_NULL_Z) = 0.0;
                }
            }
        }
        else
        {
            for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; i != r_model_part.NodesEnd() ; ++i)
            {
                if( i->HasDofFor(DISPLACEMENT_X))
                {
                    i->GetSolutionStepValue(DISPLACEMENT_OLD_X) = i->GetSolutionStepValue(DISPLACEMENT_X);

                    if(CurrentProcessInfo[FIRST_TIME_STEP])
                    {
                        i->GetSolutionStepValue(ACCELERATION_NULL_X) = i->GetSolutionStepValue(ACCELERATION_X);
                        i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_X) = i->GetSolutionStepValue(DISPLACEMENT_DT_X);
                        i->GetSolutionStepValue(DISPLACEMENT_NULL_X) = i->GetSolutionStepValue(DISPLACEMENT_X);
                    }
                    else
                    {
                        i->GetSolutionStepValue(DISPLACEMENT_NULL_X) = i->GetSolutionStepValue(DISPLACEMENT_EINS_X);
                        i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_X) = i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_X);
                        i->GetSolutionStepValue(ACCELERATION_NULL_X) = i->GetSolutionStepValue(ACCELERATION_EINS_X);
                    }

                    // // here we update the current values at the end of time step
                    // i->GetSolutionStepValue(DISPLACEMENT_X) = i->GetSolutionStepValue(DISPLACEMENT_EINS_X);
                    // i->GetSolutionStepValue(DISPLACEMENT_DT_X) = i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_X);
                    // i->GetSolutionStepValue(ACCELERATION_X) = i->GetSolutionStepValue(ACCELERATION_EINS_X);
                }
                if( i->HasDofFor(DISPLACEMENT_Y) )
                {
                    i->GetSolutionStepValue(DISPLACEMENT_OLD_Y) = i->GetSolutionStepValue(DISPLACEMENT_Y);

                    // KRATOS_WATCH(CurrentProcessInfo[FIRST_TIME_STEP])
                    // std::cout << "Node " << i->Id() << " DISPLACEMENT_NULL_Y: " << i->GetSolutionStepValue(DISPLACEMENT_NULL_Y) << std::endl;
                    // std::cout << "Node " << i->Id() << " DISPLACEMENT_EINS_Y: " << i->GetSolutionStepValue(DISPLACEMENT_EINS_Y) << std::endl;
                    // std::cout << "Node " << i->Id() << " DISPLACEMENT_Y: " << i->GetSolutionStepValue(DISPLACEMENT_Y) << std::endl;

                    if(CurrentProcessInfo[FIRST_TIME_STEP])
                    {
                        i->GetSolutionStepValue(ACCELERATION_NULL_Y) = i->GetSolutionStepValue(ACCELERATION_Y);
                        i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_Y) = i->GetSolutionStepValue(DISPLACEMENT_DT_Y);
                        i->GetSolutionStepValue(DISPLACEMENT_NULL_Y) = i->GetSolutionStepValue(DISPLACEMENT_Y);
                    }
                    else
                    {
                        i->GetSolutionStepValue(DISPLACEMENT_NULL_Y) = i->GetSolutionStepValue(DISPLACEMENT_EINS_Y);
                        i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_Y) = i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_Y);
                        i->GetSolutionStepValue(ACCELERATION_NULL_Y) = i->GetSolutionStepValue(ACCELERATION_EINS_Y);
                    }

                    // // here we update the current values at the end of time step
                    // i->GetSolutionStepValue(DISPLACEMENT_Y) = i->GetSolutionStepValue(DISPLACEMENT_EINS_Y);
                    // i->GetSolutionStepValue(DISPLACEMENT_DT_Y) = i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_Y);
                    // i->GetSolutionStepValue(ACCELERATION_Y) = i->GetSolutionStepValue(ACCELERATION_EINS_Y);
                }
                if( i->HasDofFor(DISPLACEMENT_Z) )
                {
                    i->GetSolutionStepValue(DISPLACEMENT_OLD_Z) = i->GetSolutionStepValue(DISPLACEMENT_Z);

                    if(CurrentProcessInfo[FIRST_TIME_STEP])
                    {
                        i->GetSolutionStepValue(ACCELERATION_NULL_Z) = i->GetSolutionStepValue(ACCELERATION_Z);
                        i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_Z) = i->GetSolutionStepValue(DISPLACEMENT_DT_Z);
                        i->GetSolutionStepValue(DISPLACEMENT_NULL_Z) = i->GetSolutionStepValue(DISPLACEMENT_Z);
                    }
                    else
                    {
                        i->GetSolutionStepValue(DISPLACEMENT_NULL_Z) = i->GetSolutionStepValue(DISPLACEMENT_EINS_Z);
                        i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_Z) = i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_Z);
                        i->GetSolutionStepValue(ACCELERATION_NULL_Z) = i->GetSolutionStepValue(ACCELERATION_EINS_Z);
                    }

                    // // here we update the current values at the end of time step
                    // i->GetSolutionStepValue(DISPLACEMENT_Z) = i->GetSolutionStepValue(DISPLACEMENT_EINS_Z);
                    // i->GetSolutionStepValue(DISPLACEMENT_DT_Z) = i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_Z);
                    // i->GetSolutionStepValue(ACCELERATION_Z) = i->GetSolutionStepValue(ACCELERATION_EINS_Z);
                }
                if( i->HasDofFor(WATER_PRESSURE))
                {
                    if(CurrentProcessInfo[FIRST_TIME_STEP])
                    {
                        i->GetSolutionStepValue(WATER_PRESSURE_NULL_ACCELERATION) = i->GetSolutionStepValue(WATER_PRESSURE_ACCELERATION);
                        i->GetSolutionStepValue(WATER_PRESSURE_NULL_DT) = i->GetSolutionStepValue(WATER_PRESSURE_DT);
                        i->GetSolutionStepValue(WATER_PRESSURE_NULL) = i->GetSolutionStepValue(WATER_PRESSURE);
                    }
                    else
                    {
                        i->GetSolutionStepValue(WATER_PRESSURE_NULL_DT) = i->GetSolutionStepValue(WATER_PRESSURE_EINS_DT);
                        i->GetSolutionStepValue(WATER_PRESSURE_NULL) = i->GetSolutionStepValue(WATER_PRESSURE_EINS);
                        i->GetSolutionStepValue(WATER_PRESSURE_NULL_ACCELERATION) = i->GetSolutionStepValue(WATER_PRESSURE_EINS_ACCELERATION);
                    }

                    // // here we update the current values at the end of time step
                    // i->GetSolutionStepValue(WATER_PRESSURE) = i->GetSolutionStepValue(WATER_PRESSURE_EINS);
                    // i->GetSolutionStepValue(WATER_PRESSURE_DT) = i->GetSolutionStepValue(WATER_PRESSURE_EINS_DT);
                    // i->GetSolutionStepValue(WATER_PRESSURE_ACCELERATION) = i->GetSolutionStepValue(WATER_PRESSURE_EINS_ACCELERATION);
                }
                if( i->HasDofFor(AIR_PRESSURE) )
                {
                    if(CurrentProcessInfo[FIRST_TIME_STEP])
                    {
                        i->GetSolutionStepValue(AIR_PRESSURE_NULL_ACCELERATION) = i->GetSolutionStepValue(AIR_PRESSURE_ACCELERATION);
                        i->GetSolutionStepValue(AIR_PRESSURE_NULL_DT) = i->GetSolutionStepValue(AIR_PRESSURE_DT);
                        i->GetSolutionStepValue(AIR_PRESSURE_NULL) = i->GetSolutionStepValue(AIR_PRESSURE);
                    }
                    else
                    {
                        i->GetSolutionStepValue(AIR_PRESSURE_NULL_DT) = i->GetSolutionStepValue(AIR_PRESSURE_EINS_DT);
                        i->GetSolutionStepValue(AIR_PRESSURE_NULL) = i->GetSolutionStepValue(AIR_PRESSURE_EINS);
                        i->GetSolutionStepValue(AIR_PRESSURE_NULL_ACCELERATION) = i->GetSolutionStepValue(AIR_PRESSURE_EINS_ACCELERATION);
                    }

                    // // here we update the current values at the end of time step
                    // i->GetSolutionStepValue(AIR_PRESSURE) = i->GetSolutionStepValue(AIR_PRESSURE_EINS);
                    // i->GetSolutionStepValue(AIR_PRESSURE_DT) = i->GetSolutionStepValue(AIR_PRESSURE_EINS_DT);
                    // i->GetSolutionStepValue(AIR_PRESSURE_ACCELERATION) = i->GetSolutionStepValue(AIR_PRESSURE_EINS_ACCELERATION);
                }

                if (mIntegrateRotation)
                {
                    if( i->HasDofFor(ROTATION_X) )
                    {
                        if(CurrentProcessInfo[FIRST_TIME_STEP])
                        {
                            i->GetSolutionStepValue(ANGULAR_ACCELERATION_NULL_X) = i->GetSolutionStepValue(ANGULAR_ACCELERATION_X);
                            i->GetSolutionStepValue(ROTATION_NULL_DT_X) = i->GetSolutionStepValue(ROTATION_DT_X);
                            i->GetSolutionStepValue(ROTATION_NULL_X) = i->GetSolutionStepValue(ROTATION_X);
                        }
                        else
                        {
                            i->GetSolutionStepValue(ROTATION_NULL_X) = i->GetSolutionStepValue(ROTATION_EINS_X);
                            i->GetSolutionStepValue(ROTATION_NULL_DT_X) = i->GetSolutionStepValue(ROTATION_EINS_DT_X);
                            i->GetSolutionStepValue(ANGULAR_ACCELERATION_NULL_X) = i->GetSolutionStepValue(ANGULAR_ACCELERATION_EINS_X);
                        }

                        // // here we update the current values at the end of time step
                        // i->GetSolutionStepValue(ROTATION_X) = i->GetSolutionStepValue(ROTATION_EINS_X);
                        // i->GetSolutionStepValue(ROTATION_DT_X) = i->GetSolutionStepValue(ROTATION_EINS_DT_X);
                        // i->GetSolutionStepValue(ANGULAR_ACCELERATION_X) = i->GetSolutionStepValue(ANGULAR_ACCELERATION_EINS_X);
                    }
                    if( i->HasDofFor(ROTATION_Y) )
                    {
                        if(CurrentProcessInfo[FIRST_TIME_STEP])
                        {
                            i->GetSolutionStepValue(ANGULAR_ACCELERATION_NULL_Y) = i->GetSolutionStepValue(ANGULAR_ACCELERATION_Y);
                            i->GetSolutionStepValue(ROTATION_NULL_DT_Y) = i->GetSolutionStepValue(ROTATION_DT_Y);
                            i->GetSolutionStepValue(ROTATION_NULL_Y) = i->GetSolutionStepValue(ROTATION_Y);
                        }
                        else
                        {
                            i->GetSolutionStepValue(ROTATION_NULL_Y) = i->GetSolutionStepValue(ROTATION_EINS_Y);
                            i->GetSolutionStepValue(ROTATION_NULL_DT_Y) = i->GetSolutionStepValue(ROTATION_EINS_DT_Y);
                            i->GetSolutionStepValue(ANGULAR_ACCELERATION_NULL_Y) = i->GetSolutionStepValue(ANGULAR_ACCELERATION_EINS_Y);
                        }

                        // // here we update the current values at the end of time step
                        // i->GetSolutionStepValue(ROTATION_Y) = i->GetSolutionStepValue(ROTATION_EINS_Y);
                        // i->GetSolutionStepValue(ROTATION_DT_Y) = i->GetSolutionStepValue(ROTATION_EINS_DT_Y);
                        // i->GetSolutionStepValue(ANGULAR_ACCELERATION_Y) = i->GetSolutionStepValue(ANGULAR_ACCELERATION_EINS_Y);
                    }
                    if( i->HasDofFor(ROTATION_Z) )
                    {
                        if(CurrentProcessInfo[FIRST_TIME_STEP])
                        {
                            i->GetSolutionStepValue(ANGULAR_ACCELERATION_NULL_Z) = i->GetSolutionStepValue(ANGULAR_ACCELERATION_Z);
                            i->GetSolutionStepValue(ROTATION_NULL_DT_Z) = i->GetSolutionStepValue(ROTATION_DT_Z);
                            i->GetSolutionStepValue(ROTATION_NULL_Z) = i->GetSolutionStepValue(ROTATION_Z);
                        }
                        else
                        {
                            i->GetSolutionStepValue(ROTATION_NULL_Z) = i->GetSolutionStepValue(ROTATION_EINS_Z);
                            i->GetSolutionStepValue(ROTATION_NULL_DT_Z) = i->GetSolutionStepValue(ROTATION_EINS_DT_Z);
                            i->GetSolutionStepValue(ANGULAR_ACCELERATION_NULL_Z) = i->GetSolutionStepValue(ANGULAR_ACCELERATION_EINS_Z);
                        }

                        // // here we update the current values at the end of time step
                        // i->GetSolutionStepValue(ROTATION_Z) = i->GetSolutionStepValue(ROTATION_EINS_Z);
                        // i->GetSolutionStepValue(ROTATION_DT_Z) = i->GetSolutionStepValue(ROTATION_EINS_DT_Z);
                        // i->GetSolutionStepValue(ANGULAR_ACCELERATION_Z) = i->GetSolutionStepValue(ANGULAR_ACCELERATION_EINS_Z);
                    }
                }

                if (mIntegrateMultiplier)
                {
                    if( i->HasDofFor(LAMBDA) )
                    {
                        i->GetSolutionStepValue(LAMBDA_OLD) = i->GetSolutionStepValue(LAMBDA);

                        if(CurrentProcessInfo[FIRST_TIME_STEP])
                        {
                            i->GetSolutionStepValue(LAMBDA_NULL_DT2) = i->GetSolutionStepValue(LAMBDA_DT2);
                            i->GetSolutionStepValue(LAMBDA_NULL_DT) = i->GetSolutionStepValue(LAMBDA_DT);
                            i->GetSolutionStepValue(LAMBDA_NULL) = i->GetSolutionStepValue(LAMBDA);
                        }
                        else
                        {
                            i->GetSolutionStepValue(LAMBDA_NULL) = i->GetSolutionStepValue(LAMBDA_EINS);
                            i->GetSolutionStepValue(LAMBDA_NULL_DT) = i->GetSolutionStepValue(LAMBDA_EINS_DT);
                            i->GetSolutionStepValue(LAMBDA_NULL_DT2) =i->GetSolutionStepValue(LAMBDA_EINS_DT2);
                        }

                        // // here we update the current values at the end of time step
                        // i->GetSolutionStepValue(LAMBDA) = i->GetSolutionStepValue(LAMBDA_EINS);
                        // i->GetSolutionStepValue(LAMBDA_DT) = i->GetSolutionStepValue(LAMBDA_EINS_DT);
                        // i->GetSolutionStepValue(LAMBDA_DT2) = i->GetSolutionStepValue(LAMBDA_EINS_DT2);
                    }
                }

                if (mIntegrateLoad)
                {
                    if( i->SolutionStepsDataHas(FACE_LOAD) )
                    {
                        if(CurrentProcessInfo[FIRST_TIME_STEP])
                        {
                            i->GetSolutionStepValue(FACE_LOAD_NULL_X) = i->GetSolutionStepValue(FACE_LOAD_X);
                            i->GetSolutionStepValue(FACE_LOAD_NULL_Y) = i->GetSolutionStepValue(FACE_LOAD_Y);
                            i->GetSolutionStepValue(FACE_LOAD_NULL_Z) = i->GetSolutionStepValue(FACE_LOAD_Z);
                        }
                        else
                        {
                            i->GetSolutionStepValue(FACE_LOAD_NULL_X) = i->GetSolutionStepValue(FACE_LOAD_EINS_X);
                            i->GetSolutionStepValue(FACE_LOAD_NULL_Y) = i->GetSolutionStepValue(FACE_LOAD_EINS_Y);
                            i->GetSolutionStepValue(FACE_LOAD_NULL_Z) = i->GetSolutionStepValue(FACE_LOAD_EINS_Z);
                        }
                    }
                    if( i->SolutionStepsDataHas(FORCE) )
                    {
                        if(CurrentProcessInfo[FIRST_TIME_STEP])
                        {
                            i->GetSolutionStepValue(FORCE_NULL_X) = i->GetSolutionStepValue(FORCE_X);
                            i->GetSolutionStepValue(FORCE_NULL_Y) = i->GetSolutionStepValue(FORCE_Y);
                            i->GetSolutionStepValue(FORCE_NULL_Z) = i->GetSolutionStepValue(FORCE_Z);
                        }
                        else
                        {
                            i->GetSolutionStepValue(FORCE_NULL_X) = i->GetSolutionStepValue(FORCE_EINS_X);
                            i->GetSolutionStepValue(FORCE_NULL_Y) = i->GetSolutionStepValue(FORCE_EINS_Y);
                            i->GetSolutionStepValue(FORCE_NULL_Z) = i->GetSolutionStepValue(FORCE_EINS_Z);
                        }
                    }
                }
            }
        }
    }
    //***************************************************************************
    //***************************************************************************

    void CalculateSystemContributions(
        Element& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) final
    {
        KRATOS_TRY

        rCurrentElement.CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);
        rCurrentElement.EquationIdVector(EquationId,CurrentProcessInfo);

        Matrix DampingMatrix;

        Matrix MassMatrix;

        rCurrentElement.CalculateDampingMatrix(DampingMatrix, CurrentProcessInfo);

        rCurrentElement.CalculateMassMatrix(MassMatrix, CurrentProcessInfo);

        // KRATOS_WATCH(rCurrentElement.Id())
        // KRATOS_WATCH(MassMatrix)
        // KRATOS_WATCH(DampingMatrix)
        // KRATOS_WATCH(LHS_Contribution)

        if ((norm_frobenius(DampingMatrix) == 0.0)
        && (((DampingMatrix.size1() != LHS_Contribution.size1())
          || (DampingMatrix.size2() != LHS_Contribution.size2()))))
        {
            if ((DampingMatrix.size1() != LHS_Contribution.size1()) || (DampingMatrix.size2() != LHS_Contribution.size2()))
                DampingMatrix.resize(LHS_Contribution.size1(), LHS_Contribution.size2(), false);
            noalias(DampingMatrix) = ZeroMatrix(LHS_Contribution.size1(), LHS_Contribution.size2());
        }

        if ((norm_frobenius(MassMatrix) == 0.0)
        && (((MassMatrix.size1() != LHS_Contribution.size1())
          || (MassMatrix.size2() != LHS_Contribution.size2()))))
        {
            if ((MassMatrix.size1() != LHS_Contribution.size1()) || (MassMatrix.size2() != LHS_Contribution.size2()))
                MassMatrix.resize(LHS_Contribution.size1(), LHS_Contribution.size2(), false);
            noalias(MassMatrix) = ZeroMatrix(LHS_Contribution.size1(), LHS_Contribution.size2());
        }

        AddDampingToRHS(rCurrentElement, RHS_Contribution, DampingMatrix, CurrentProcessInfo);

        if (mForceLumpedMass)
            this->LumpMatrix(MassMatrix);

        AssembleTimeSpaceLHS_Dynamics(LHS_Contribution, DampingMatrix, MassMatrix, CurrentProcessInfo);

        // if (rCurrentElement.Id() == 743)
        // {
        //     Vector u;
        //     rCurrentElement.GetValuesVector(u, 0);
        //     KRATOS_WATCH(u)
        //     Vector RHS(u.size());
        //     noalias(RHS) = prod( LHS_Contribution, u );
        //     KRATOS_WATCH(RHS)
        //     KRATOS_WATCH(RHS_Contribution)
        // }

        // account for prescription of dofs
        try
        {
            const PrescribedObject& rObject = dynamic_cast<PrescribedObject&>(rCurrentElement);
            rObject.ApplyPrescribedDofs(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);
        }
        catch (std::bad_cast& bc)
        {
            // DO NOTHING
        }

        KRATOS_CATCH("")
    }

    void CalculateRHSContribution(
        Element& rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) final
    {
        KRATOS_TRY

        rCurrentElement.CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId,CurrentProcessInfo);

        Matrix DampingMatrix;

        Matrix MassMatrix;

        rCurrentElement.CalculateDampingMatrix(DampingMatrix, CurrentProcessInfo);

        rCurrentElement.CalculateMassMatrix(MassMatrix, CurrentProcessInfo);

        if ((norm_frobenius(DampingMatrix) == 0.0)
        && (((DampingMatrix.size1() != RHS_Contribution.size())
          || (DampingMatrix.size2() != RHS_Contribution.size()))))
        {
            if ((DampingMatrix.size1() != RHS_Contribution.size()) || (DampingMatrix.size2() != RHS_Contribution.size()))
                DampingMatrix.resize(RHS_Contribution.size(), RHS_Contribution.size(), false);
            noalias(DampingMatrix) = ZeroMatrix(RHS_Contribution.size(), RHS_Contribution.size());
        }

        if ((norm_frobenius(MassMatrix) == 0.0)
        && (((MassMatrix.size1() != RHS_Contribution.size())
          || (MassMatrix.size2() != RHS_Contribution.size()))))
        {
            if ((MassMatrix.size1() != RHS_Contribution.size()) || (MassMatrix.size2() != RHS_Contribution.size()))
                MassMatrix.resize(RHS_Contribution.size(), RHS_Contribution.size(), false);
            noalias(MassMatrix) = ZeroMatrix(RHS_Contribution.size(), RHS_Contribution.size());
        }

        if (mForceLumpedMass)
            this->LumpMatrix(MassMatrix);

        AddDampingToRHS(rCurrentElement, RHS_Contribution, DampingMatrix, CurrentProcessInfo);

        KRATOS_CATCH("")
    }


    void CalculateLHSContribution(
        Element& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) final
    {
        KRATOS_TRY

        rCurrentElement.CalculateLeftHandSide(LHS_Contribution, CurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId,CurrentProcessInfo);

        Matrix DampingMatrix;

        Matrix MassMatrix;

        rCurrentElement.CalculateDampingMatrix(DampingMatrix, CurrentProcessInfo);

        rCurrentElement.CalculateMassMatrix(MassMatrix, CurrentProcessInfo);

        if ((norm_frobenius(DampingMatrix) == 0.0)
        && (((DampingMatrix.size1() != LHS_Contribution.size1())
          || (DampingMatrix.size2() != LHS_Contribution.size2()))))
        {
            if ((DampingMatrix.size1() != LHS_Contribution.size1()) || (DampingMatrix.size2() != LHS_Contribution.size2()))
                DampingMatrix.resize(LHS_Contribution.size1(), LHS_Contribution.size2(), false);
            noalias(DampingMatrix) = ZeroMatrix(LHS_Contribution.size1(), LHS_Contribution.size2());
        }

        if ((norm_frobenius(MassMatrix) == 0.0)
        && (((MassMatrix.size1() != LHS_Contribution.size1())
          || (MassMatrix.size2() != LHS_Contribution.size2()))))
        {
            if ((MassMatrix.size1() != LHS_Contribution.size1()) || (MassMatrix.size2() != LHS_Contribution.size2()))
                MassMatrix.resize(LHS_Contribution.size1(), LHS_Contribution.size2(), false);
            noalias(MassMatrix) = ZeroMatrix(LHS_Contribution.size1(), LHS_Contribution.size2());
        }

        if (mForceLumpedMass)
            this->LumpMatrix(MassMatrix);

        AssembleTimeSpaceLHS_Dynamics(LHS_Contribution, DampingMatrix, MassMatrix, CurrentProcessInfo);

        KRATOS_CATCH("")
    }


    void CalculateSystemContributions(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) final
    {
        KRATOS_TRY

        Matrix DampingMatrix;

        rCurrentCondition.CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

        rCurrentCondition.EquationIdVector(EquationId, CurrentProcessInfo);

        rCurrentCondition.CalculateDampingMatrix(DampingMatrix, CurrentProcessInfo);

        Matrix MassMatrix = ZeroMatrix(LHS_Contribution.size1(), LHS_Contribution.size2());

        if (mForceLumpedMass)
            this->LumpMatrix(MassMatrix);

        AssembleTimeSpaceLHS_Dynamics(LHS_Contribution, DampingMatrix, MassMatrix, CurrentProcessInfo);

        KRATOS_CATCH("")
    }


    void CalculateLHSContribution(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) final
    {
        KRATOS_TRY

        rCurrentCondition.CalculateLeftHandSide(LHS_Contribution,CurrentProcessInfo);
        rCurrentCondition.EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH("")
    }


    void CalculateRHSContribution(
        Condition& rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) final
    {
        KRATOS_TRY

        rCurrentCondition.CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);
        rCurrentCondition.EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH("")
    }


    void UpdateForces(ModelPart& r_model_part )
    {
        KRATOS_TRY

        const ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        const double Dt = CurrentProcessInfo[DELTA_TIME];

        //Update nodal values and nodal velocities
        for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; i != r_model_part.NodesEnd() ; i++)
        {
            if( i->SolutionStepsDataHas(FACE_LOAD) )
            {
                noalias( i->GetSolutionStepValue(FACE_LOAD) )
                = i->GetSolutionStepValue(FACE_LOAD_NULL);
            }

            if( i->SolutionStepsDataHas(FORCE) )
            {
                noalias( i->GetSolutionStepValue(FORCE) )
                = i->GetSolutionStepValue(FORCE_NULL);
            }
        }

        KRATOS_CATCH("")
    }

    /*@} */
    /**@name Access */
    /*@{ */


    /*@} */
    /**@name Inquiry */
    /*@{ */


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ResidualBasedAccBasedCentralDifferenceScheme";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    /*@} */
    /**@name Friends */
    /*@{ */

protected:

    void AddDampingToRHS(
        Element& rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        LocalSystemMatrixType& D,
        const ProcessInfo& CurrentProcessInfo)
    {
        const double Dt = CurrentProcessInfo[DELTA_TIME];

        // adding damping contribution
        Vector ud;
        rCurrentElement.GetFirstDerivativesVector(ud, 0);
        noalias(RHS_Contribution) -= prod(D, ud);
    }

    void AssembleTimeSpaceLHS_Dynamics(
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemMatrixType& DampingMatrix,
        LocalSystemMatrixType& MassMatrix,
        const ProcessInfo& CurrentProcessInfo)
    {
        // set LHS as mass matrix
        noalias(LHS_Contribution) = MassMatrix;
    }

    /*@} */
    /**@name Protected member Variables */
    /*@{ */
    /*@} */
    /**@name Protected Operators*/
    /*@{ */
    /*@} */
    /**@name Protected Operations*/
    /*@{ */
    /*@} */
    /**@name Protected  Access */
    /*@{ */
    /*@} */
    /**@name Protected Inquiry */
    /*@{ */
    /*@} */
    /**@name Protected LifeCycle */
    /*@{ */
private:
    /**@name Static Member Variables */
    /*@{ */
    /*@} */
    /**@name Member Variables */
    /*@{ */

    bool mIntegrateRotation;
    bool mIntegrateMultiplier;
    bool mIntegrateLoad;
    bool mForceLumpedMass;

    void LumpMatrix(Matrix& MassMatrix) const
    {
        for (std::size_t i = 0; i < MassMatrix.size1(); ++i)
        {
            double row_sum = 0.0;
            for (std::size_t j = 0; j < MassMatrix.size1(); ++j)
            {
                row_sum += MassMatrix(i, j);
                MassMatrix(i, j) = 0.0;
            }

            MassMatrix(i, i) = row_sum;
        }
    }

    /*@} */
    /**@name Private Operators*/
    /*@{ */
    /*@} */
    /**@name Private Operations*/
    /*@{ */
    /*@} */
    /**@name Private  Access */
    /*@{ */
    /*@} */
    /**@name Private Inquiry */
    /*@{ */
    /*@} */
    /**@name Un accessible methods */
    /*@{ */
}; /* Class Scheme */
}  /* namespace Kratos.*/

#endif /* KRATOS_STRUCTURAL_APPLICATION_RESIDUALBASED_ACC_BASED_CENTRAL_DIFFERENCE_SCHEME  defined */

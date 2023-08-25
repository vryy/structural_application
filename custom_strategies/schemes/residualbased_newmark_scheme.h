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
*   Last Modified by:    $Author: nagel $
*   Date:                $Date: 2009-03-20 08:56:48 $
*   Revision:            $Revision: 1.7 $
*   Date:                $Date: 23/07/2023 $
*   Revision:            $Revision: 1.8 $
*
* ***********************************************************/


#if !defined(KRATOS_STRUCTURAL_APPLICATION_RESIDUALBASED_NEWMARK_SCHEME_H_INCLUDED )
#define  KRATOS_STRUCTURAL_APPLICATION_RESIDUALBASED_NEWMARK_SCHEME_H_INCLUDED

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "includes/element.h"
#include "includes/fnv_1a_hash.h"
#ifdef SD_APP_FORWARD_COMPATIBILITY
#include "custom_python3/legacy_structural_app_vars.h"
#else
#include "includes/legacy_structural_app_vars.h"
#endif
#include "containers/array_1d.h"
#include "solving_strategies/schemes/scheme.h"
#include "structural_application_variables.h"
#include "custom_elements/prescribed_object.h"
#include "residualbased_newmark_helper.h"

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
 * Implementation of Time Integration scheme based on Generalized-Alpha Method
 */
template<class TSparseSpace, class TDenseSpace, int TMassDampingType = 0>
class ResidualBasedNewmarkScheme : public Scheme<TSparseSpace,TDenseSpace>
{
public:
    /**@name Type Definitions */
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedNewmarkScheme );

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

    typedef ResidualBasedNewmarkHelper<TMassDampingType> ResidualBasedNewmarkHelperType;

    /**
     * Constructor for PURE Newmark scheme
     * @ref Chung&Hulbert: "A time integration algorithm for structural dynamics with improved
     *  numerical dissipation: The generalized alpha method" Journal of applied mechanics, 60, 371-375
     */
    ResidualBasedNewmarkScheme() : BaseType()
    {
        // pure Newmark scheme
        mAlpha_f = 0.0;
        mAlpha_m = 0.0;
        mBeta = 1.0 / 4.0;
        mGamma = 0.5;

        // set all the integration flag to false. To incorporate the corresponding term, the flag shall be set by user.
        mIntegrateRotation = false;
        mIntegrateMultiplier = false;
        mIntegrateLoad = false;

        std::cout << "PURE Newmark Time Integration Scheme!!!!!!!!!!!!!!!!!!!!!" << " alpha_f= " << mAlpha_f << " alpha_m= " << mAlpha_m << " beta= " << mBeta << " gamma= " << mGamma << std::endl;
        std::cout << "Nonlinear mass damping: " << TMassDampingType << std::endl;
        //(...)_NULL DOF at the begin of time step, (...)_EINS DOF at the end of time step, (...)
        // DOF at the midpoint of time step, Please recognize that the name of the DOF is (...)
        // while the iteration is done towards the (...)_EINS DOF value at end of time step
    }

    /**
     * Constructor.
     * @param mDissipationRadius if the scheme is numerically energy conserving or not
     *                              == 1.0 energy conserving
     *                              < 1.0 numerically discipating
     * @ref Chung&Hulbert: "A time integration algorithm for structural dynamics with improved
     *  numerical dissipation: The generalized alpha method" Journal of applied mechanics, 60, 371-375
     */
    ResidualBasedNewmarkScheme(double mDissipationRadius ) : BaseType()
    {
        // standard Newmark scheme
        mAlpha_f = mDissipationRadius / (1.0 + mDissipationRadius);
        mAlpha_m = (2.0 * mDissipationRadius - 1.0) / (mDissipationRadius + 1.0);
        mBeta = (1.0 + mAlpha_f - mAlpha_m) * (1.0 + mAlpha_f - mAlpha_m)/4.0;
        mGamma = 0.5 + mAlpha_f - mAlpha_m;

        // set all the integration flag to false. To incorporate the corresponding term, the flag shall be set by user.
        mIntegrateRotation = false;
        mIntegrateMultiplier = false;
        mIntegrateLoad = false;

        std::cout << "Using the Generalized alpha Time Integration Scheme with radius= "<< mDissipationRadius << " alpha_f= " << mAlpha_f << " alpha_m= " << mAlpha_m << " beta= " << mBeta << " gamma= " << mGamma << std::endl;
        std::cout << "Nonlinear mass damping: " << TMassDampingType << std::endl;
    }

    /**
     * Constructor with full options.
     * @param mDissipationRadius if the scheme is numerically energy conserving or not
     *                              == 1.0 energy conserving
     *                              < 1.0 numerically discipating
     * @ref Chung&Hulbert: "A time integration algorithm for structural dynamics with improved
     *  numerical dissipation: The generalized alpha method" Journal of applied mechanics, 60, 371-375
     */
    ResidualBasedNewmarkScheme(int option, double mDissipationRadius ) : BaseType()
    {
        // set all the integration flag to false. To incorporate the corresponding term, the flag shall be set by user.
        mIntegrateRotation = false;
        mIntegrateMultiplier = false;
        mIntegrateLoad = false;

        if (option == 0)
        {
            mAlpha_f = 0.0;
            mAlpha_m = 0.0;
            mBeta = 1.0 / 4.0;
            mGamma = 0.5;

            std::cout << "PURE Newmark Time Integration Scheme!!!!!!!!!!!!!!!!!!!!!" << " alpha_f= " << mAlpha_f << " alpha_m= " << mAlpha_m << " beta= " << mBeta << " gamma= " << mGamma << std::endl;
        }
        else if (option == 1)
        {
            mAlpha_f = mDissipationRadius / (1.0 + mDissipationRadius);
            mAlpha_m = (2.0 * mDissipationRadius - 1.0) / (mDissipationRadius + 1.0);
            mBeta = (1.0 + mAlpha_f - mAlpha_m) * (1.0 + mAlpha_f - mAlpha_m)/4.0;
            mGamma = 0.5 + mAlpha_f - mAlpha_m;

            std::cout << "using the Generalized-alpha Time Integration Scheme with radius= "<< mDissipationRadius << " alpha_f= " << mAlpha_f << " alpha_m= " << mAlpha_m << " beta= " << mBeta << " gamma= " << mGamma << std::endl;
        }
        else if (option == 2)
        {
            mAlpha_f = 0.0;
            mAlpha_m = (2.0 * mDissipationRadius - 1.0) / (mDissipationRadius + 1.0);
            mBeta = (1.0 - mAlpha_m) * (1.0 - mAlpha_m)/4.0;
            mGamma = 0.5 - mAlpha_m;

            std::cout << "using the Bossak-alpha Time Integration Scheme with radius= "<< mDissipationRadius << " alpha_m= " << mAlpha_m << " beta= " << mBeta << " gamma= " << mGamma << std::endl;
        }
        else if (option == 3)
        {
            mAlpha_f = mDissipationRadius / (1.0 + mDissipationRadius);
            mAlpha_m = 0.0;
            mBeta = (1.0 + mAlpha_f) * (1.0 + mAlpha_f)/4.0;
            mGamma = 0.5 + mAlpha_f;

            std::cout << "using the Hilber-alpha Time Integration Scheme with radius= "<< mDissipationRadius << " alpha_f= " << mAlpha_f << " beta= " << mBeta << " gamma= " << mGamma << std::endl;
        }
        else
        {
            KRATOS_THROW_ERROR(std::logic_error, "Invalid Time integration option", option)
        }

        std::cout << "Nonlinear mass damping: " << TMassDampingType << std::endl;
    }

    /** Destructor.*/
    virtual ~ResidualBasedNewmarkScheme()
    {}

    /**@name Operators */

    /// Copy assignment
    ResidualBasedNewmarkScheme& operator=(const ResidualBasedNewmarkScheme& rOther)
    {
        this->mAlpha_f = rOther.mAlpha_f;
        this->mAlpha_m = rOther.mAlpha_m;
        this->mBeta = rOther.mBeta;
        this->mGamma = rOther.mGamma;
        this->mIntegrateRotation = rOther.mIntegrateRotation;
        this->mIntegrateMultiplier = mIntegrateMultiplier;
        this->mIntegrateLoad = mIntegrateLoad;
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

    /// Initialize the scheme
    void Initialize(ModelPart& r_model_part) override
    {
        BaseType::Initialize(r_model_part);

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        CurrentProcessInfo[TIME_INTEGRATION_SCHEME] = FNV1a32Hash::CalculateHash(Info().c_str());
        CurrentProcessInfo[NEWMARK_ALPHAF] = mAlpha_f;
        CurrentProcessInfo[NEWMARK_ALPHAM] = mAlpha_m;
        CurrentProcessInfo[NEWMARK_BETA]   = mBeta;
        CurrentProcessInfo[NEWMARK_GAMMA]  = mGamma;

        std::cout << "ModelPart " << r_model_part.Name() << " is initialized by " << Info() << std::endl;
        KRATOS_WATCH(CurrentProcessInfo[TIME_INTEGRATION_SCHEME])
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
                    i->GetSolutionStepValue(DISPLACEMENT_EINS_X)
                    +=Dx[i->GetDof(DISPLACEMENT_X).EquationId()];
                }
            }
            if( i->HasDofFor(DISPLACEMENT_Y) )
            {
                if(i->GetDof(DISPLACEMENT_Y).IsFree())
                {
                    i->GetSolutionStepValue(DISPLACEMENT_EINS_Y)
                    +=Dx[i->GetDof(DISPLACEMENT_Y).EquationId()];
                }
            }
            if( i->HasDofFor(DISPLACEMENT_Z) )
            {
                if(i->GetDof(DISPLACEMENT_Z).IsFree())
                {
                    i->GetSolutionStepValue(DISPLACEMENT_EINS_Z)
                    +=Dx[i->GetDof(DISPLACEMENT_Z).EquationId()];
                }
            }
            if( i->HasDofFor(WATER_PRESSURE) )
            {
                if(i->GetDof(WATER_PRESSURE).IsFree())
                {
                    i->GetSolutionStepValue(WATER_PRESSURE_EINS)
                    +=Dx[i->GetDof(WATER_PRESSURE).EquationId()];
                }
            }
            if( i->HasDofFor(AIR_PRESSURE) )
            {
                if(i->GetDof(AIR_PRESSURE).IsFree())
                {
                    i->GetSolutionStepValue(AIR_PRESSURE_EINS)
                    +=Dx[i->GetDof(AIR_PRESSURE).EquationId()];
                }
            }
            if( i->HasDofFor(TEMPERATURE) )
            {
                if(i->GetDof(TEMPERATURE).IsFree())
                {
                    i->GetSolutionStepValue(TEMPERATURE_EINS)
                    +=Dx[i->GetDof(TEMPERATURE).EquationId()];
                }
            }
            if( i->HasDofFor(DENSITY) )
            {
                if(i->GetDof(DENSITY).IsFree())
                {
                    i->GetSolutionStepValue(DENSITY_EINS)
                    +=Dx[i->GetDof(DENSITY).EquationId()];
                }
            }

            // update this variables to account for Lagrange multiplier
            if( i->HasDofFor(LAGRANGE_DISPLACEMENT_X) )
            {
                if(i->GetDof(LAGRANGE_DISPLACEMENT_X).IsFree())
                {
                    i->GetSolutionStepValue(LAGRANGE_DISPLACEMENT_X)
                    +=Dx[i->GetDof(LAGRANGE_DISPLACEMENT_X).EquationId()];
                }
            }
            if( i->HasDofFor(LAGRANGE_DISPLACEMENT_Y) )
            {
                if(i->GetDof(LAGRANGE_DISPLACEMENT_Y).IsFree())
                {
                    i->GetSolutionStepValue(LAGRANGE_DISPLACEMENT_Y)
                    +=Dx[i->GetDof(LAGRANGE_DISPLACEMENT_Y).EquationId()];
                }
            }
            if( i->HasDofFor(LAGRANGE_DISPLACEMENT_Z) )
            {
                if(i->GetDof(LAGRANGE_DISPLACEMENT_Z).IsFree())
                {
                    i->GetSolutionStepValue(LAGRANGE_DISPLACEMENT_Z)
                    +=Dx[i->GetDof(LAGRANGE_DISPLACEMENT_Z).EquationId()];
                }
            }
            if( i->HasDofFor(LAGRANGE_WATER_PRESSURE) )
            {
                if(i->GetDof(LAGRANGE_WATER_PRESSURE).IsFree())
                {
                    i->GetSolutionStepValue(LAGRANGE_WATER_PRESSURE)
                    +=Dx[i->GetDof(LAGRANGE_WATER_PRESSURE).EquationId()];
                }
            }
            if( i->HasDofFor(LAGRANGE_AIR_PRESSURE) )
            {
                if(i->GetDof(LAGRANGE_AIR_PRESSURE).IsFree())
                {
                    i->GetSolutionStepValue(LAGRANGE_AIR_PRESSURE)
                    +=Dx[i->GetDof(LAGRANGE_AIR_PRESSURE).EquationId()];
                }
            }
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
                    rNode.GetSolutionStepValue(DISPLACEMENT_EINS_X)
                    += Dx[dof_iterator->EquationId()];
                }
            }
            else if (dof_iterator->GetVariable() == DISPLACEMENT_Y)
            {
                if (dof_iterator->IsFree())
                {
                    rNode.GetSolutionStepValue(DISPLACEMENT_EINS_Y)
                    += Dx[dof_iterator->EquationId()];
                }
            }
            else if (dof_iterator->GetVariable() == DISPLACEMENT_Z)
            {
                if (dof_iterator->IsFree())
                {
                    rNode.GetSolutionStepValue(DISPLACEMENT_EINS_Z)
                    += Dx[dof_iterator->EquationId()];
                }
            }
            else if (dof_iterator->GetVariable() == WATER_PRESSURE)
            {
                if (dof_iterator->IsFree())
                {
                    rNode.GetSolutionStepValue(WATER_PRESSURE_EINS)
                    += Dx[dof_iterator->EquationId()];
                }
            }
            else if (dof_iterator->GetVariable() == AIR_PRESSURE)
            {
                if (dof_iterator->IsFree())
                {
                    rNode.GetSolutionStepValue(AIR_PRESSURE_EINS)
                    += Dx[dof_iterator->EquationId()];
                }
            }
            else if (dof_iterator->GetVariable() == TEMPERATURE)
            {
                if (dof_iterator->IsFree())
                {
                    rNode.GetSolutionStepValue(TEMPERATURE_EINS)
                    += Dx[dof_iterator->EquationId()];
                }
            }
            else if (dof_iterator->GetVariable() == DENSITY)
            {
                if (dof_iterator->IsFree())
                {
                    rNode.GetSolutionStepValue(DENSITY_EINS)
                    += Dx[dof_iterator->EquationId()];
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
                            rNode.GetSolutionStepValue(ROTATION_EINS_X)
                            += Dx[dof_iterator->EquationId()];
                        }
                    }
                    else if (dof_iterator->GetVariable() == ROTATION_Y)
                    {
                        if (dof_iterator->IsFree())
                        {
                            rNode.GetSolutionStepValue(ROTATION_EINS_Y)
                            += Dx[dof_iterator->EquationId()];
                        }
                    }
                    else if (dof_iterator->GetVariable() == ROTATION_Z)
                    {
                        if (dof_iterator->IsFree())
                        {
                            rNode.GetSolutionStepValue(ROTATION_EINS_Z)
                            += Dx[dof_iterator->EquationId()];
                        }
                    }
                }

                if (mIntegrateMultiplier)
                {
                    if (dof_iterator->GetVariable() == LAMBDA)
                    {
                        if (dof_iterator->IsFree())
                        {
                            rNode.GetSolutionStepValue(LAMBDA_EINS)
                            += Dx[dof_iterator->EquationId()];
                        }
                    }
                }

                if (!mIntegrateMultiplier && !mIntegrateRotation)
                {
                    // update for non-dynamics variable
                    if (dof_iterator->IsFree())
                    {
                        dof_iterator->GetSolutionStepValue()
                        += Dx[dof_iterator->EquationId()];
                    }
                }
            }
        }

        KRATOS_CATCH("")
    }

    /**
     * initializes next newton step by calculating
     * u_{n+1-alpha}= u_n*alpha+U_n+1^k+1*(1-alpha)
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

        const double Dt = CurrentProcessInfo[DELTA_TIME];

        //Update nodal values and nodal velocities at mAlpha_f
        for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; i != r_model_part.NodesEnd() ; i++)
        {
            if( i->HasDofFor(DISPLACEMENT_X) )
            {
                i->GetSolutionStepValue(ACCELERATION_EINS_X)
                = 1.0/(mBeta*Dt*Dt)
                  * (i->GetSolutionStepValue(DISPLACEMENT_EINS_X)
                     -i->GetSolutionStepValue(DISPLACEMENT_NULL_X))
                  -1.0/(mBeta*Dt)
                  *i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_X)
                  -(1.0-2.0*mBeta)/(2.0*mBeta)*i->GetSolutionStepValue(ACCELERATION_NULL_X);

                i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_X)
                = (i->GetSolutionStepValue(DISPLACEMENT_EINS_X)
                   -i->GetSolutionStepValue(DISPLACEMENT_NULL_X))
                  *mGamma/(mBeta*Dt)
                  -(mGamma-mBeta)/mBeta*(i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_X))
                  -(mGamma-2.0*mBeta)/(2.0*mBeta)*Dt
                  *(i->GetSolutionStepValue(ACCELERATION_NULL_X));

                i->GetSolutionStepValue(ACCELERATION_X)
                = mAlpha_m*i->GetSolutionStepValue(ACCELERATION_NULL_X)
                  +(1.0-mAlpha_m)*i->GetSolutionStepValue(ACCELERATION_EINS_X);

                i->GetSolutionStepValue(DISPLACEMENT_DT_X)
                = mAlpha_f*i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_X)
                  +(1.0-mAlpha_f)*i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_X);

                i->GetSolutionStepValue(DISPLACEMENT_X)
                = mAlpha_f*i->GetSolutionStepValue(DISPLACEMENT_NULL_X)
                  +(1.0-mAlpha_f)*i->GetSolutionStepValue(DISPLACEMENT_EINS_X);
            }
            if( i->HasDofFor(DISPLACEMENT_Y) )
            {
                i->GetSolutionStepValue(ACCELERATION_EINS_Y)
                = 1.0/(mBeta*Dt*Dt)
                  * (i->GetSolutionStepValue(DISPLACEMENT_EINS_Y)
                    -i->GetSolutionStepValue(DISPLACEMENT_NULL_Y))
                 -1.0/(mBeta*Dt)
                 *i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_Y)
                 -(1.0-2.0*mBeta)/(2.0*mBeta)*
                 i->GetSolutionStepValue(ACCELERATION_NULL_Y);

                i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_Y)
                =(i->GetSolutionStepValue(DISPLACEMENT_EINS_Y)
                  -i->GetSolutionStepValue(DISPLACEMENT_NULL_Y))
                 *mGamma/(mBeta*Dt)
                 -(mGamma-mBeta)/mBeta*(i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_Y))
                 -(mGamma-2.0*mBeta)/(2.0*mBeta)*Dt
                 *(i->GetSolutionStepValue(ACCELERATION_NULL_Y));

                i->GetSolutionStepValue(ACCELERATION_Y)
                =mAlpha_m*i->GetSolutionStepValue(ACCELERATION_NULL_Y)
                 +(1.0-mAlpha_m)*i->GetSolutionStepValue(ACCELERATION_EINS_Y);

                i->GetSolutionStepValue(DISPLACEMENT_DT_Y)
                = mAlpha_f*i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_Y)
                  +(1.0-mAlpha_f)*i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_Y);

                i->GetSolutionStepValue(DISPLACEMENT_Y)
                = mAlpha_f*i->GetSolutionStepValue(DISPLACEMENT_NULL_Y)+(1.0-mAlpha_f)*
                  i->GetSolutionStepValue(DISPLACEMENT_EINS_Y);
            }
            if( i->HasDofFor(DISPLACEMENT_Z) )
            {
                i->GetSolutionStepValue(ACCELERATION_EINS_Z)
                = 1.0/(mBeta*Dt*Dt)
                  * (i->GetSolutionStepValue(DISPLACEMENT_EINS_Z)
                     -i->GetSolutionStepValue(DISPLACEMENT_NULL_Z))
                  -1.0/(mBeta*Dt)
                  *i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_Z)
                  -(1.0-2.0*mBeta)/(2.0*mBeta)*i->GetSolutionStepValue(ACCELERATION_NULL_Z);

                i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_Z)
                = (i->GetSolutionStepValue(DISPLACEMENT_EINS_Z)
                   -i->GetSolutionStepValue(DISPLACEMENT_NULL_Z))
                  *mGamma/(mBeta*Dt)-(mGamma-mBeta)/mBeta*
                  (i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_Z))
                  -(mGamma-2.0*mBeta)/(2.0*mBeta)*Dt
                  *(i->GetSolutionStepValue(ACCELERATION_NULL_Z));

                i->GetSolutionStepValue(ACCELERATION_Z)
                = mAlpha_m*i->GetSolutionStepValue(ACCELERATION_NULL_Z)
                  +(1.0-mAlpha_m)*i->GetSolutionStepValue(ACCELERATION_EINS_Z);

                i->GetSolutionStepValue(DISPLACEMENT_DT_Z)
                = mAlpha_f*i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_Z)
                  +(1.0-mAlpha_f)*i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_Z);

                i->GetSolutionStepValue(DISPLACEMENT_Z)
                = mAlpha_f*i->GetSolutionStepValue(DISPLACEMENT_NULL_Z)
                  +(1.0-mAlpha_f)*i->GetSolutionStepValue(DISPLACEMENT_EINS_Z);
            }
            if( i->HasDofFor(WATER_PRESSURE) )
            {
                i->GetSolutionStepValue(WATER_PRESSURE_EINS_ACCELERATION)
                = 1.0/(mBeta*Dt*Dt)
                  * (i->GetSolutionStepValue(WATER_PRESSURE_EINS)
                     -i->GetSolutionStepValue(WATER_PRESSURE_NULL))
                  -1.0/(mBeta*Dt)
                  *i->GetSolutionStepValue(WATER_PRESSURE_NULL_DT)
                  -(1.0-2.0*mBeta)/(2.0*mBeta)*i->GetSolutionStepValue(WATER_PRESSURE_NULL_ACCELERATION);

                i->GetSolutionStepValue(WATER_PRESSURE_EINS_DT)
                = (i->GetSolutionStepValue(WATER_PRESSURE_EINS)
                   -i->GetSolutionStepValue(WATER_PRESSURE_NULL))
                  *mGamma/(mBeta*Dt)
                  -(mGamma-mBeta)/mBeta*(i->GetSolutionStepValue(WATER_PRESSURE_NULL_DT))
                  // TODO check if the inertial term is also needed, see the update for TEMPERATURE below
                  -(mGamma-2.0*mBeta)/(2.0*mBeta)*Dt
                  *(i->GetSolutionStepValue(WATER_PRESSURE_NULL_ACCELERATION));

                i->GetSolutionStepValue(WATER_PRESSURE_ACCELERATION)
                = mAlpha_m*i->GetSolutionStepValue(WATER_PRESSURE_NULL_ACCELERATION)
                  +(1.0-mAlpha_m)*i->GetSolutionStepValue(WATER_PRESSURE_EINS_ACCELERATION);

                i->GetSolutionStepValue(WATER_PRESSURE_DT)
                = mAlpha_f*i->GetSolutionStepValue(WATER_PRESSURE_NULL_DT)
                  +(1.0-mAlpha_f)*i->GetSolutionStepValue(WATER_PRESSURE_EINS_DT);

                i->GetSolutionStepValue(WATER_PRESSURE)
                = mAlpha_f*i->GetSolutionStepValue(WATER_PRESSURE_NULL)
                  +(1.0-mAlpha_f)*i->GetSolutionStepValue(WATER_PRESSURE_EINS);
            }
            if( i->HasDofFor(AIR_PRESSURE) )
            {
                i->GetSolutionStepValue(AIR_PRESSURE_EINS_ACCELERATION)
                = 1.0/(mBeta*Dt*Dt)
                  * (i->GetSolutionStepValue(AIR_PRESSURE_EINS)
                     -i->GetSolutionStepValue(AIR_PRESSURE_NULL))
                  -1.0/(mBeta*Dt)
                  *i->GetSolutionStepValue(AIR_PRESSURE_NULL_DT)
                  -(1.0-2.0*mBeta)/(2.0*mBeta)*i->GetSolutionStepValue(AIR_PRESSURE_NULL_ACCELERATION);

                i->GetSolutionStepValue(AIR_PRESSURE_EINS_DT)
                = (i->GetSolutionStepValue(AIR_PRESSURE_EINS)
                   -i->GetSolutionStepValue(AIR_PRESSURE_NULL))
                  *mGamma/(mBeta*Dt)
                  -(mGamma-mBeta)/mBeta*
                  (i->GetSolutionStepValue(AIR_PRESSURE_NULL_DT))
                  // TODO check if the inertial term is also needed, see the update for TEMPERATURE below
                  -(mGamma-2.0*mBeta)/(2.0*mBeta)*Dt
                  *(i->GetSolutionStepValue(AIR_PRESSURE_NULL_ACCELERATION));

                i->GetSolutionStepValue(AIR_PRESSURE_ACCELERATION)
                = mAlpha_m*i->GetSolutionStepValue(AIR_PRESSURE_NULL_ACCELERATION)
                  +(1.0-mAlpha_m)*i->GetSolutionStepValue(AIR_PRESSURE_EINS_ACCELERATION);

                i->GetSolutionStepValue(AIR_PRESSURE_DT)
                = mAlpha_f*i->GetSolutionStepValue(AIR_PRESSURE_NULL_DT)
                  +(1.0-mAlpha_f)*i->GetSolutionStepValue(AIR_PRESSURE_EINS_DT);

                i->GetSolutionStepValue(AIR_PRESSURE)
                = mAlpha_f*i->GetSolutionStepValue(AIR_PRESSURE_NULL)
                  +(1.0-mAlpha_f)*i->GetSolutionStepValue(AIR_PRESSURE_EINS);
            }
            if( i->HasDofFor(TEMPERATURE) )
            {
                i->GetSolutionStepValue(TEMPERATURE_EINS_ACCELERATION)
                = 1.0/(mBeta*Dt*Dt)
                  * (i->GetSolutionStepValue(TEMPERATURE_EINS)
                     -i->GetSolutionStepValue(TEMPERATURE_NULL))
                  -1.0/(mBeta*Dt)
                  *i->GetSolutionStepValue(TEMPERATURE_NULL_DT)
                  -(1.0-2.0*mBeta)/(2.0*mBeta)*i->GetSolutionStepValue(TEMPERATURE_NULL_ACCELERATION);

                i->GetSolutionStepValue(TEMPERATURE_EINS_DT)
                = (i->GetSolutionStepValue(TEMPERATURE_EINS)
                   -i->GetSolutionStepValue(TEMPERATURE_NULL))
                  *mGamma/(mBeta*Dt)
                  -(mGamma-mBeta)/mBeta*
                  (i->GetSolutionStepValue(TEMPERATURE_NULL_DT));
                // temperature equation is a parabolic equation. Therefore we
                // don't need to add the inertial term to the approximation of
                // the rate of temperature
                  // -(mGamma-2.0*mBeta)/(2.0*mBeta)*Dt
                  // *(i->GetSolutionStepValue(TEMPERATURE_NULL_ACCELERATION));

                i->GetSolutionStepValue(TEMPERATURE_ACCELERATION)
                = mAlpha_m*i->GetSolutionStepValue(TEMPERATURE_NULL_ACCELERATION)
                  +(1.0-mAlpha_m)*i->GetSolutionStepValue(TEMPERATURE_EINS_ACCELERATION);

                i->GetSolutionStepValue(TEMPERATURE_DT)
                = mAlpha_f*i->GetSolutionStepValue(TEMPERATURE_NULL_DT)
                  +(1.0-mAlpha_f)*i->GetSolutionStepValue(TEMPERATURE_EINS_DT);

                i->GetSolutionStepValue(TEMPERATURE)
                = mAlpha_f*i->GetSolutionStepValue(TEMPERATURE_NULL)
                  +(1.0-mAlpha_f)*i->GetSolutionStepValue(TEMPERATURE_EINS);
            }
            if( i->HasDofFor(DENSITY) )
            {
                i->GetSolutionStepValue(DENSITY_EINS_ACCELERATION)
                = 1.0/(mBeta*Dt*Dt)
                  * (i->GetSolutionStepValue(DENSITY_EINS)
                     -i->GetSolutionStepValue(DENSITY_NULL))
                  -1.0/(mBeta*Dt)
                  *i->GetSolutionStepValue(DENSITY_NULL_DT)
                  -(1.0-2.0*mBeta)/(2.0*mBeta)*i->GetSolutionStepValue(DENSITY_NULL_ACCELERATION);

                i->GetSolutionStepValue(DENSITY_EINS_DT)
                = (i->GetSolutionStepValue(DENSITY_EINS)
                   -i->GetSolutionStepValue(DENSITY_NULL))
                  *mGamma/(mBeta*Dt)
                  -(mGamma-mBeta)/mBeta*
                  (i->GetSolutionStepValue(DENSITY_NULL_DT));
                // DENSITY equation is a parabolic equation. Therefore we
                // don't need to add the inertial term to the approximation of
                // the rate of DENSITY
                  // -(mGamma-2.0*mBeta)/(2.0*mBeta)*Dt
                  // *(i->GetSolutionStepValue(DENSITY_NULL_ACCELERATION));

                i->GetSolutionStepValue(DENSITY_ACCELERATION)
                = mAlpha_m*i->GetSolutionStepValue(DENSITY_NULL_ACCELERATION)
                  +(1.0-mAlpha_m)*i->GetSolutionStepValue(DENSITY_EINS_ACCELERATION);

                i->GetSolutionStepValue(DENSITY_DT)
                = mAlpha_f*i->GetSolutionStepValue(DENSITY_NULL_DT)
                  +(1.0-mAlpha_f)*i->GetSolutionStepValue(DENSITY_EINS_DT);

                i->GetSolutionStepValue(DENSITY)
                = mAlpha_f*i->GetSolutionStepValue(DENSITY_NULL)
                  +(1.0-mAlpha_f)*i->GetSolutionStepValue(DENSITY_EINS);
            }

            if (mIntegrateRotation)
            {
                if( i->HasDofFor(ROTATION_X) )
                {
                    i->GetSolutionStepValue(ANGULAR_ACCELERATION_EINS_X)
                    = 1.0/(mBeta*Dt*Dt)
                      * (i->GetSolutionStepValue(ROTATION_EINS_X)
                         -i->GetSolutionStepValue(ROTATION_NULL_X))
                      -1.0/(mBeta*Dt)
                      *i->GetSolutionStepValue(ROTATION_NULL_DT_X)
                      -(1.0-2.0*mBeta)/(2.0*mBeta)*i->GetSolutionStepValue(ANGULAR_ACCELERATION_NULL_X);

                    i->GetSolutionStepValue(ROTATION_EINS_DT_X)
                    = (i->GetSolutionStepValue(ROTATION_EINS_X)
                       -i->GetSolutionStepValue(ROTATION_NULL_X))
                      *mGamma/(mBeta*Dt)
                      -(mGamma-mBeta)/mBeta*(i->GetSolutionStepValue(ROTATION_NULL_DT_X))
                      -(mGamma-2.0*mBeta)/(2.0*mBeta)*Dt
                      *(i->GetSolutionStepValue(ANGULAR_ACCELERATION_NULL_X));

                    i->GetSolutionStepValue(ANGULAR_ACCELERATION_X)
                    = mAlpha_m*i->GetSolutionStepValue(ANGULAR_ACCELERATION_NULL_X)
                      +(1.0-mAlpha_m)*i->GetSolutionStepValue(ANGULAR_ACCELERATION_EINS_X);

                    i->GetSolutionStepValue(ROTATION_DT_X)
                    = mAlpha_f*i->GetSolutionStepValue(ROTATION_NULL_DT_X)
                      +(1.0-mAlpha_f)*i->GetSolutionStepValue(ROTATION_EINS_DT_X);

                    i->GetSolutionStepValue(ROTATION_X)
                    = mAlpha_f*i->GetSolutionStepValue(ROTATION_NULL_X)
                      +(1.0-mAlpha_f)*i->GetSolutionStepValue(ROTATION_EINS_X);

                }
                if( i->HasDofFor(ROTATION_Y) )
                {
                    i->GetSolutionStepValue(ANGULAR_ACCELERATION_EINS_Y)
                    =1.0/(mBeta*Dt*Dt)
                     * (i->GetSolutionStepValue(ROTATION_EINS_Y)
                        -i->GetSolutionStepValue(ROTATION_NULL_Y))
                     -1.0/(mBeta*Dt)
                     *i->GetSolutionStepValue(ROTATION_NULL_DT_Y)
                     -(1.0-2.0*mBeta)/(2.0*mBeta)*
                     i->GetSolutionStepValue(ANGULAR_ACCELERATION_NULL_Y);

                    i->GetSolutionStepValue(ROTATION_EINS_DT_Y)
                    =(i->GetSolutionStepValue(ROTATION_EINS_Y)
                      -i->GetSolutionStepValue(ROTATION_NULL_Y))
                     *mGamma/(mBeta*Dt)
                     -(mGamma-mBeta)/mBeta*(i->GetSolutionStepValue(ROTATION_NULL_DT_Y))
                     -(mGamma-2.0*mBeta)/(2.0*mBeta)*Dt
                     *(i->GetSolutionStepValue(ANGULAR_ACCELERATION_NULL_Y));

                    i->GetSolutionStepValue(ANGULAR_ACCELERATION_Y)
                    =mAlpha_m*i->GetSolutionStepValue(ANGULAR_ACCELERATION_NULL_Y)
                     +(1.0-mAlpha_m)*i->GetSolutionStepValue(ANGULAR_ACCELERATION_EINS_Y);

                    i->GetSolutionStepValue(ROTATION_DT_Y)
                    = mAlpha_f*i->GetSolutionStepValue(ROTATION_NULL_DT_Y)
                      +(1.0-mAlpha_f)*i->GetSolutionStepValue(ROTATION_EINS_DT_Y);

                    i->GetSolutionStepValue(ROTATION_Y)
                    = mAlpha_f*i->GetSolutionStepValue(ROTATION_NULL_Y)+(1.0-mAlpha_f)*
                      i->GetSolutionStepValue(ROTATION_EINS_Y);
                }
                if( i->HasDofFor(ROTATION_Z) )
                {
                    i->GetSolutionStepValue(ANGULAR_ACCELERATION_EINS_Z)
                    = 1.0/(mBeta*Dt*Dt)
                      * (i->GetSolutionStepValue(ROTATION_EINS_Z)
                         -i->GetSolutionStepValue(ROTATION_NULL_Z))
                      -1.0/(mBeta*Dt)
                      *i->GetSolutionStepValue(ROTATION_NULL_DT_Z)
                      -(1.0-2.0*mBeta)/(2.0*mBeta)*i->GetSolutionStepValue(ANGULAR_ACCELERATION_NULL_Z);

                    i->GetSolutionStepValue(ROTATION_EINS_DT_Z)
                    = (i->GetSolutionStepValue(ROTATION_EINS_Z)
                       -i->GetSolutionStepValue(ROTATION_NULL_Z))
                      *mGamma/(mBeta*Dt)-(mGamma-mBeta)/mBeta*
                      (i->GetSolutionStepValue(ROTATION_NULL_DT_Z))
                      -(mGamma-2.0*mBeta)/(2.0*mBeta)*Dt
                      *(i->GetSolutionStepValue(ANGULAR_ACCELERATION_NULL_Z));

                    i->GetSolutionStepValue(ANGULAR_ACCELERATION_Z)
                    = mAlpha_m*i->GetSolutionStepValue(ANGULAR_ACCELERATION_NULL_Z)
                      +(1.0-mAlpha_m)*i->GetSolutionStepValue(ANGULAR_ACCELERATION_EINS_Z);

                    i->GetSolutionStepValue(ROTATION_DT_Z)
                    = mAlpha_f*i->GetSolutionStepValue(ROTATION_NULL_DT_Z)
                      +(1.0-mAlpha_f)*i->GetSolutionStepValue(ROTATION_EINS_DT_Z);

                    i->GetSolutionStepValue(ROTATION_Z)
                    = mAlpha_f*i->GetSolutionStepValue(ROTATION_NULL_Z)
                      +(1.0-mAlpha_f)*i->GetSolutionStepValue(ROTATION_EINS_Z);
                }
            }

            if (mIntegrateMultiplier)
            {
                if( i->HasDofFor(LAMBDA) )
                {
                    i->GetSolutionStepValue(LAMBDA_EINS_DT2)
                    = 1.0/(mBeta*Dt*Dt)
                      * (i->GetSolutionStepValue(LAMBDA_EINS)
                         -i->GetSolutionStepValue(LAMBDA_NULL))
                      -1.0/(mBeta*Dt)
                      *i->GetSolutionStepValue(LAMBDA_NULL_DT)
                      -(1.0-2.0*mBeta)/(2.0*mBeta)*i->GetSolutionStepValue(LAMBDA_NULL_DT2);

                    i->GetSolutionStepValue(LAMBDA_EINS_DT)
                    = (i->GetSolutionStepValue(LAMBDA_EINS)
                       -i->GetSolutionStepValue(LAMBDA_NULL))
                      *mGamma/(mBeta*Dt)
                      -(mGamma-mBeta)/mBeta*(i->GetSolutionStepValue(LAMBDA_NULL_DT))
                      -(mGamma-2.0*mBeta)/(2.0*mBeta)*Dt
                      *(i->GetSolutionStepValue(LAMBDA_NULL_DT2));

                    i->GetSolutionStepValue(LAMBDA_DT2)
                    = mAlpha_m*i->GetSolutionStepValue(LAMBDA_NULL_DT2)
                      +(1.0-mAlpha_m)*i->GetSolutionStepValue(LAMBDA_EINS_DT2);

                    i->GetSolutionStepValue(LAMBDA_DT)
                    = mAlpha_f*i->GetSolutionStepValue(LAMBDA_NULL_DT)
                      +(1.0-mAlpha_f)*i->GetSolutionStepValue(LAMBDA_EINS_DT);

                    i->GetSolutionStepValue(LAMBDA)
                    = mAlpha_f*i->GetSolutionStepValue(LAMBDA_NULL)
                      +(1.0-mAlpha_f)*i->GetSolutionStepValue(LAMBDA_EINS);
                }
            }

            if (mIntegrateLoad)
            {
                if( i->SolutionStepsDataHas(FACE_LOAD) )
                {
                    noalias( i->GetSolutionStepValue(FACE_LOAD) )
                    = mAlpha_f*i->GetSolutionStepValue(FACE_LOAD_NULL)
                      +(1.0-mAlpha_f)*i->GetSolutionStepValue(FACE_LOAD_EINS);
                }

                if( i->SolutionStepsDataHas(FORCE) )
                {
                    noalias( i->GetSolutionStepValue(FORCE) )
                    = mAlpha_f*i->GetSolutionStepValue(FORCE_NULL)
                      +(1.0-mAlpha_f)*i->GetSolutionStepValue(FORCE_EINS);
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
        //Update nodal values and nodal velocities at mAlpha_f
        for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ;
                i != r_model_part.NodesEnd() ; i++)
        {
            if( i->HasDofFor(DISPLACEMENT_X) && i->GetDof(DISPLACEMENT_X).IsFree())
            {
                i->GetSolutionStepValue(ACCELERATION_EINS_X)=
                    i->GetSolutionStepValue(ACCELERATION_NULL_X);
                i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_X)=
                    i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_X);
                i->GetSolutionStepValue(DISPLACEMENT_EINS_X )=
                    i->GetSolutionStepValue(DISPLACEMENT_NULL_X);
            }
            if( i->HasDofFor(DISPLACEMENT_Y) && i->GetDof(DISPLACEMENT_Y).IsFree())
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
            if( i->HasDofFor(TEMPERATURE) && i->GetDof(TEMPERATURE).IsFree())
            {
                i->GetSolutionStepValue(TEMPERATURE_EINS_DT)=
                    i->GetSolutionStepValue(TEMPERATURE_NULL_DT);
                i->GetSolutionStepValue(TEMPERATURE_EINS_ACCELERATION)=
                    i->GetSolutionStepValue(TEMPERATURE_NULL_ACCELERATION);
                i->GetSolutionStepValue(TEMPERATURE_EINS)=
                    i->GetSolutionStepValue(TEMPERATURE_NULL);
            }
            if( i->HasDofFor(DENSITY) && i->GetDof(DENSITY).IsFree())
            {
                i->GetSolutionStepValue(DENSITY_EINS_DT)=
                    i->GetSolutionStepValue(DENSITY_NULL_DT);
                i->GetSolutionStepValue(DENSITY_EINS_ACCELERATION)=
                    i->GetSolutionStepValue(DENSITY_NULL_ACCELERATION);
                i->GetSolutionStepValue(DENSITY_EINS)=
                    i->GetSolutionStepValue(DENSITY_NULL);
            }

            if (mIntegrateRotation)
            {
                if( i->HasDofFor(ROTATION_X) && i->GetDof(ROTATION_X).IsFree())
                {
                    i->GetSolutionStepValue(ANGULAR_ACCELERATION_EINS_X)=
                        i->GetSolutionStepValue(ANGULAR_ACCELERATION_NULL_X);
                    i->GetSolutionStepValue(ROTATION_EINS_DT_X)=
                        i->GetSolutionStepValue(ROTATION_NULL_DT_X);
                    i->GetSolutionStepValue(ROTATION_EINS_X)=
                        i->GetSolutionStepValue(ROTATION_NULL_X);
                }
                if( i->HasDofFor(ROTATION_Y) && i->GetDof(ROTATION_Y).IsFree())
                {
                    i->GetSolutionStepValue(ANGULAR_ACCELERATION_EINS_Y)=
                        i->GetSolutionStepValue(ANGULAR_ACCELERATION_NULL_Y);
                    i->GetSolutionStepValue(ROTATION_EINS_DT_Y)=
                        i->GetSolutionStepValue(ROTATION_NULL_DT_Y);
                    i->GetSolutionStepValue(ROTATION_EINS_Y)=
                        i->GetSolutionStepValue(ROTATION_NULL_Y);
                }
                if( i->HasDofFor(ROTATION_Z) && i->GetDof(ROTATION_Z).IsFree())
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
                if( i->HasDofFor(LAMBDA) && i->GetDof(LAMBDA).IsFree())
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
                        i->GetSolutionStepValue(ACCELERATION_NULL_X) =i->GetSolutionStepValue(ACCELERATION_X);
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
                if( i->HasDofFor(TEMPERATURE) )
                {
                    if(CurrentProcessInfo[FIRST_TIME_STEP])
                    {
                        i->GetSolutionStepValue(TEMPERATURE_NULL_ACCELERATION) = i->GetSolutionStepValue(TEMPERATURE_ACCELERATION);
                        i->GetSolutionStepValue(TEMPERATURE_NULL_DT) = i->GetSolutionStepValue(TEMPERATURE_DT);
                        i->GetSolutionStepValue(TEMPERATURE_NULL) = i->GetSolutionStepValue(TEMPERATURE);
                    }
                    else
                    {
                        i->GetSolutionStepValue(TEMPERATURE_NULL_DT) = i->GetSolutionStepValue(TEMPERATURE_EINS_DT);
                        i->GetSolutionStepValue(TEMPERATURE_NULL) = i->GetSolutionStepValue(TEMPERATURE_EINS);
                        i->GetSolutionStepValue(TEMPERATURE_NULL_ACCELERATION) = i->GetSolutionStepValue(TEMPERATURE_EINS_ACCELERATION);
                    }

                    // // here we update the current values at the end of time step
                    // i->GetSolutionStepValue(TEMPERATURE) = i->GetSolutionStepValue(TEMPERATURE_EINS);
                    // i->GetSolutionStepValue(TEMPERATURE_DT) = i->GetSolutionStepValue(TEMPERATURE_EINS_DT);
                    // i->GetSolutionStepValue(TEMPERATURE_ACCELERATION) = i->GetSolutionStepValue(TEMPERATURE_EINS_ACCELERATION);
                }
                if( i->HasDofFor(DENSITY) )
                {
                    if(CurrentProcessInfo[FIRST_TIME_STEP])
                    {
                        i->GetSolutionStepValue(DENSITY_NULL_ACCELERATION) = i->GetSolutionStepValue(DENSITY_ACCELERATION);
                        i->GetSolutionStepValue(DENSITY_NULL_DT) = i->GetSolutionStepValue(DENSITY_DT);
                        i->GetSolutionStepValue(DENSITY_NULL) = i->GetSolutionStepValue(DENSITY);
                    }
                    else
                    {
                        i->GetSolutionStepValue(DENSITY_NULL_DT) = i->GetSolutionStepValue(DENSITY_EINS_DT);
                        i->GetSolutionStepValue(DENSITY_NULL) = i->GetSolutionStepValue(DENSITY_EINS);
                        i->GetSolutionStepValue(DENSITY_NULL_ACCELERATION) = i->GetSolutionStepValue(DENSITY_EINS_ACCELERATION);
                    }

                    // // here we update the current values at the end of time step
                    // i->GetSolutionStepValue(DENSITY) = i->GetSolutionStepValue(DENSITY_EINS);
                    // i->GetSolutionStepValue(DENSITY_DT) = i->GetSolutionStepValue(DENSITY_EINS_DT);
                    // i->GetSolutionStepValue(DENSITY_ACCELERATION) = i->GetSolutionStepValue(DENSITY_EINS_ACCELERATION);
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
                        // if(CurrentProcessInfo[FIRST_TIME_STEP])
                        // {
                        //     i->GetSolutionStepValue(FACE_LOAD_NULL_X) = i->GetSolutionStepValue(FACE_LOAD_X);
                        //     i->GetSolutionStepValue(FACE_LOAD_NULL_Y) = i->GetSolutionStepValue(FACE_LOAD_Y);
                        //     i->GetSolutionStepValue(FACE_LOAD_NULL_Z) = i->GetSolutionStepValue(FACE_LOAD_Z);
                        // }
                        // else
                        // {
                            i->GetSolutionStepValue(FACE_LOAD_NULL_X) = i->GetSolutionStepValue(FACE_LOAD_EINS_X);
                            i->GetSolutionStepValue(FACE_LOAD_NULL_Y) = i->GetSolutionStepValue(FACE_LOAD_EINS_Y);
                            i->GetSolutionStepValue(FACE_LOAD_NULL_Z) = i->GetSolutionStepValue(FACE_LOAD_EINS_Z);
                        // }
                    }
                    if( i->SolutionStepsDataHas(FORCE) )
                    {
                        // if(CurrentProcessInfo[FIRST_TIME_STEP])
                        // {
                        //     i->GetSolutionStepValue(FORCE_NULL_X) = i->GetSolutionStepValue(FORCE_X);
                        //     i->GetSolutionStepValue(FORCE_NULL_Y) = i->GetSolutionStepValue(FORCE_Y);
                        //     i->GetSolutionStepValue(FORCE_NULL_Z) = i->GetSolutionStepValue(FORCE_Z);
                        // }
                        // else
                        // {
                            i->GetSolutionStepValue(FORCE_NULL_X) = i->GetSolutionStepValue(FORCE_EINS_X);
                            i->GetSolutionStepValue(FORCE_NULL_Y) = i->GetSolutionStepValue(FORCE_EINS_Y);
                            i->GetSolutionStepValue(FORCE_NULL_Z) = i->GetSolutionStepValue(FORCE_EINS_Z);
                        // }
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
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        ResidualBasedNewmarkHelperType::CalculateSystemContributions( rCurrentElement,
                LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo );

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
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        ResidualBasedNewmarkHelperType::CalculateRHSContribution( rCurrentElement,
                RHS_Contribution, EquationId, CurrentProcessInfo );

        KRATOS_CATCH("")
    }


    void CalculateLHSContribution(
        Element& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        ResidualBasedNewmarkHelperType::CalculateLHSContribution( rCurrentElement,
                LHS_Contribution, EquationId, CurrentProcessInfo );

        KRATOS_CATCH("")
    }


    void CalculateSystemContributions(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        ResidualBasedNewmarkHelperType::CalculateSystemContributions( rCurrentCondition,
                LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo );

        KRATOS_CATCH("")
    }


    void CalculateLHSContribution(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        ResidualBasedNewmarkHelperType::CalculateLHSContribution( rCurrentCondition,
                LHS_Contribution, EquationId, CurrentProcessInfo );

        KRATOS_CATCH("")
    }


    void CalculateRHSContribution(
        Condition& rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        ResidualBasedNewmarkHelperType::CalculateRHSContribution( rCurrentCondition,
                RHS_Contribution, EquationId, CurrentProcessInfo );

        KRATOS_CATCH("")
    }


    void UpdateForces(ModelPart& r_model_part )
    {
        KRATOS_TRY

        const ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        const double Dt = CurrentProcessInfo[DELTA_TIME];

        //Update nodal values and nodal velocities at mAlpha_f
        for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; i != r_model_part.NodesEnd() ; i++)
        {
            if( i->SolutionStepsDataHas(FACE_LOAD) )
            {
                noalias( i->GetSolutionStepValue(FACE_LOAD) )
                = mAlpha_f*i->GetSolutionStepValue(FACE_LOAD_NULL)
                  +(1.0-mAlpha_f)*i->GetSolutionStepValue(FACE_LOAD_EINS);
            }

            if( i->SolutionStepsDataHas(FORCE) )
            {
                noalias( i->GetSolutionStepValue(FORCE) )
                = mAlpha_f*i->GetSolutionStepValue(FORCE_NULL)
                  +(1.0-mAlpha_f)*i->GetSolutionStepValue(FORCE_EINS);
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
        std::stringstream ss;
        ss << "ResidualBasedNewmarkScheme<" << TMassDampingType << ">";
        return ss.str();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << " alpha_f: " << mAlpha_f
                 << " alpha_m: " << mAlpha_m
                 << " beta: " << mBeta
                 << " gamma: " << mGamma;
    }

    /*@} */
    /**@name Friends */
    /*@{ */

protected:
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

    double mAlpha_f;
    double mAlpha_m;
    double mBeta;
    double mGamma;

    bool mIntegrateRotation;
    bool mIntegrateMultiplier;
    bool mIntegrateLoad;

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
}; /* class ResidualBasedNewmarkScheme */

}  /* namespace Kratos.*/

#endif /* KRATOS_STRUCTURAL_APPLICATION_RESIDUALBASED_NEWMARK_SCHEME_H_INCLUDED  defined */

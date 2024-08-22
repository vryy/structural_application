/*
==============================================================================
see soil_mechanics_appplication/LICENSE.txt
==============================================================================
*/
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Giang Bui-Hoang $
//   Date:                $Date: 16 Feb 2020 $
//   Revision:            $Revision: 1.0 $
//
//



#if !defined(KRATOS_SOIL_MECHANICS_GENERAL_ELASTOPLASTICITY_LAW_H_INCLUDED )
#define  KRATOS_SOIL_MECHANICS_GENERAL_ELASTOPLASTICITY_LAW_H_INCLUDED

/* System includes */

/* External includes */

/* Project includes */
#include "plasticity_laws/general_plasticity_law.h"

// #define DEBUG_GENERAL_PLASTICITY_LAW

namespace Kratos
{

/**
 * Asbtract class for the general elastoplasticity law.
 * It is noted that q denotes the stress-like internal variables, not deviatoric pressure
 * REF: gen_plas.pdf
 */
class GeneralElastoplasticityLaw : virtual public GeneralPlasticityLaw
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(GeneralElastoplasticityLaw);
    typedef GeneralPlasticityLaw BaseType;
    typedef BaseType::Third_Order_Tensor Third_Order_Tensor;
    typedef BaseType::Fourth_Order_Tensor Fourth_Order_Tensor;

    /**
     * Constructor.
     */
    GeneralElastoplasticityLaw()
    {}

    /**
     * Destructor.
     */
    virtual ~GeneralElastoplasticityLaw()
    {}

    /**
     * Operations
     */


    ///////////////////////////////////////////////////////////
    ///////////// PLASTIC CALCULATION SUBROUTINES /////////////
    ///////////////////////////////////////////////////////////

    /// Perform the plastic integration with sub-stepping for elastoplastic material law
    /// stress, q and alpha should be stress_n, q_n and alpha_n. On output, the values at n+1 will be returned.
    /// On return, the vector of converged delta factor is returned
    std::vector<double> PlasticIntegration_Substepping(Matrix& stress, Vector& q, Vector& alpha,
        const Matrix& incremental_strain, const Fourth_Order_Tensor& Ce,
        const double FTOL, const int max_iters,
        const ProcessInfo& CurrentProcessInfo, const Properties& props,
        const int debug_level=0) const;

    /// Perform the plastic integration with sub-stepping on a specific loading profile for elastoplastic material law
    /// stress, q and alpha should be stress_n, q_n and alpha_n. On output, the values at n+1 will be returned.
    /// On return, the vector of converged delta factor is returned
    int PlasticIntegration_Substepping(Matrix& stress, Vector& q, Vector& alpha,
        const std::vector<double>& loads, const Matrix& incremental_strain, const Fourth_Order_Tensor& Ce,
        const double FTOL, const int max_iters,
        const ProcessInfo& CurrentProcessInfo, const Properties& props,
        const int debug_level=0) const;

    /// Perform the plastic integration for elastoplastic constitutive law
    /// q and alpha should be q_n and alpha_n. On output, the values at n+1 will be returned.
    int PlasticIntegration(Matrix& stress, Vector& q, Vector& alpha, double& dlambda,
        const Matrix& stress_trial, const Fourth_Order_Tensor& Ce,
        const double FTOL, const int max_iters,
        const ProcessInfo& CurrentProcessInfo, const Properties& props,
        const int debug_level=0) const;

    /// Perform the plastic integration for elastoplastic constitutive law using cutting plane algorithm
    /// q and alpha should be q_n and alpha_n. On output, the values at n+1 will be returned.
    int PlasticIntegration_CuttingPlane(Matrix& stress, Vector& q, Vector& alpha, double& dlambda,
        const Matrix& stress_trial, const Fourth_Order_Tensor& Ce,
        const double FTOL, const int max_iters,
        const ProcessInfo& CurrentProcessInfo, const Properties& props,
        const int debug_level=0) const;

    /// Compute the (continuum) plastic tangent
    void ComputeContinuumPlasticTangent(Fourth_Order_Tensor& Cep,
        const Fourth_Order_Tensor& Ce, const Matrix& stress, const Vector& q, const Vector& alpha,
        const ProcessInfo& CurrentProcessInfo, const Properties& props) const;

    /// Compute the (consistent) plastic tangent according to the procedure in plastic integration
    void ComputeConsistentPlasticTangent(Fourth_Order_Tensor& Cep,
        const Fourth_Order_Tensor& Ce, const Matrix& stress, const Vector& q, const Vector& alpha,
        const double dlambda,
        const ProcessInfo& CurrentProcessInfo, const Properties& props) const;

    /// Compute the (consistent) plastic tangent according to the procedure in plastic integration with sub-stepping
    void ComputeConsistentPlasticTangent_Substepping(Fourth_Order_Tensor& Cep,
        const Fourth_Order_Tensor& Ce, const Matrix& stressn, const Vector& qn, const Vector& alphan,
        const std::vector<double>& loads, const Matrix& incremental_strain,
        const double FTOL, const int max_iters,
        const ProcessInfo& CurrentProcessInfo, const Properties& props) const;

    /// Compute the (numerical) plastic tangent according to the procedure in plastic integration with sub-stepping
    void ComputeNumericalPlasticTangent_Substepping(Fourth_Order_Tensor& Cep,
        const Fourth_Order_Tensor& Ce, const Matrix& stressn, const Vector& qn, const Vector& alphan,
        const std::vector<double>& loads, const Matrix& incremental_strain, const double epsilon,
        const double FTOL, const int max_iters,
        const ProcessInfo& CurrentProcessInfo, const Properties& props) const;

    ///////////////////////////////////////////////////////////

    /**
     * Turn back information as a string.
     */
    std::string Name() const override
    {
        return "GeneralElastoplasticityLaw";
    }

private:

    /// Compute the stress providing dlambda
    int PlasticIntegration_ComputeStress(Matrix& stress, const Vector& q, const Vector& alpha,
        const double dlambda, const Matrix& stress_trial, const Fourth_Order_Tensor& Ce,
        const double FTOL, const int max_iters,
        const ProcessInfo& CurrentProcessInfo, const Properties& props) const;

    double PlasticIntegration_ComputeRHS(const Matrix& stress, const Vector& q, const Vector& alpha,
        const double lambda,
        const Matrix& stress_trial, const Fourth_Order_Tensor& Ce,
        const ProcessInfo& CurrentProcessInfo, const Properties& props) const;

    double PlasticIntegration_ComputeNumLHS(const Matrix& stress, const Vector& q, const Vector& alpha,
        const double dlambda, const double ddlambda,
        const Matrix& stress_trial, const Fourth_Order_Tensor& Ce,
        const double FTOL, const int max_iters,
        const ProcessInfo& CurrentProcessInfo, const Properties& props) const;

}; /* Class GeneralElastoplasticityLaw */

} /* namespace Kratos.*/

#endif /* KRATOS_SOIL_MECHANICS_GENERAL_ELASTOPLASTICITY_LAW_H_INCLUDED  defined */

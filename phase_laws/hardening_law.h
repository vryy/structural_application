//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 5 Dec 2018 $
//   Revision:            $Revision: 1.0 $
//
//



#if !defined(KRATOS_STRUCTURAL_APPLICATION_HARDENING_LAW_H_INCLUDED )
#define  KRATOS_STRUCTURAL_APPLICATION_HARDENING_LAW_H_INCLUDED

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/variables.h"
#include "includes/properties.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"
#include "containers/flags.h"


namespace Kratos
{

/**
 * Asbtract class for the hardening law.
 */
class HardeningLaw : public Flags
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(HardeningLaw);

    /**
     * Constructor.
     */
    HardeningLaw()
    {}

    /**
     * Destructor.
     */
    virtual ~HardeningLaw()
    {}

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     * NOTE: implementation scheme:
     *      HardeningLaw::Pointer p_clone(new HardeningLaw());
     *      return p_clone;
     */
    virtual HardeningLaw::Pointer Clone() const
    {
        return HardeningLaw::Pointer(new HardeningLaw());
    }

    /**
     * Operations
     */

    virtual bool Has( const Variable<int>& rThisVariable ) {}
    virtual bool Has( const Variable<double>& rThisVariable ) {}
    virtual bool Has( const Variable<Vector>& rThisVariable ) {}
    virtual bool Has( const Variable<Matrix>& rThisVariable ) {}

    virtual int& GetValue( const Variable<int>& rThisVariable, int& rValue ) {}
    virtual double& GetValue( const Variable<double>& rThisVariable, double& rValue ) {}
    virtual Vector& GetValue( const Variable<Vector>& rThisVariable, Vector& rValue ) {}
    virtual Matrix& GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue ) {}

    virtual void SetValue( const Variable<int>& rThisVariable, const int rValue,
                   const ProcessInfo& rCurrentProcessInfo ) {}
    virtual void SetValue( const Variable<double>& rThisVariable, const double rValue,
                   const ProcessInfo& rCurrentProcessInfo ) {}
    virtual void SetValue( const Variable<array_1d<double, 3> >& rThisVariable,
                   const array_1d<double, 3>& rValue, const ProcessInfo& rCurrentProcessInfo ) {}
    virtual void SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) {}
    virtual void SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) {}

    /// Get the value of the hardening function w.r.t constitutive variable
    virtual double GetValue(const double phi) const {return 0.0;}

    /// Get the value of the hardening function w.r.t constitutive variables
    virtual double GetValue(const std::vector<double>& values) const {return 0.0;}

    /// Get the derivative of the hardening function w.r.t constitutive variable
    virtual double GetDerivative(const double phi) const {return 0.0;}

    /// Get the partial derivative of the hardening function w.r.t constitutive variable
    virtual double GetDerivative(const unsigned int i, const std::vector<double>& values) const {return 0.0;}

    /// Utility function to set the hardening_law to the properties
    static void Assign( const Variable<HardeningLaw::Pointer>& rThisVariable, const HardeningLaw::Pointer pValue,
                const Properties::Pointer pProperties )
    {
        pProperties->SetValue(rThisVariable, pValue);
    }

    /// Fit a specific hardening curve
    virtual void Fit(const Vector& x, const Vector& y) const
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "HardeningLaw";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

private:

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags );
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags );
    }

    ///@}

}; /* Class HardeningLaw */

/// input stream function
inline std::istream & operator >>(std::istream& rIStream,
                                  HardeningLaw& rThis);

/// output stream function

inline std::ostream & operator <<(std::ostream& rOStream,
                                  const HardeningLaw& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

//template class KRATOS_API(KRATOS_CORE) KratosComponents<HardeningLaw>;

//void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, HardeningLaw const& ThisComponent);

/**
 * Definition of HardeningLaw variable
 */

// #undef  KRATOS_EXPORT_MACRO
// #define KRATOS_EXPORT_MACRO KRATOS_API

// KRATOS_DEFINE_VARIABLE(HardeningLaw::Pointer, HARDENING_LAW)
// KRATOS_DEFINE_VARIABLE(HardeningLaw::Pointer, HARDENING_LAW_1)
// KRATOS_DEFINE_VARIABLE(HardeningLaw::Pointer, HARDENING_LAW_2)
// KRATOS_DEFINE_VARIABLE(HardeningLaw::Pointer, ISOTROPIC_HARDENING_LAW)
// KRATOS_DEFINE_VARIABLE(HardeningLaw::Pointer, KINEMATIC_HARDENING_LAW)

// #undef  KRATOS_EXPORT_MACRO
// #define KRATOS_EXPORT_MACRO KRATOS_NO_EXPORT

} /* namespace Kratos.*/

#endif /* KRATOS_STRUCTURAL_APPLICATION_HARDENING_LAW_H_INCLUDED  defined */

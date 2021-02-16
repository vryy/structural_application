//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 22 Jun 20 $
//   Revision:            $Revision: 1.0 $
//
//
#if !defined(KRATOS_TIMOSHENKO_LINEAR_BEAM_ELEMENT_H_INCLUDED )
#define  KRATOS_TIMOSHENKO_LINEAR_BEAM_ELEMENT_H_INCLUDED


// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

namespace Kratos
{

/**
 */
class TimoshenkoLinearBeamElement : public Element
{
public:
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef GeometryType::IntegrationPointType IntegrationPointType;

    typedef ConstitutiveLaw ConstitutiveLawType;

    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;

    // Counted pointer of TimoshenkoLinearBeamElement
    KRATOS_CLASS_POINTER_DEFINITION(TimoshenkoLinearBeamElement);

    /**
     * Default constructor.
     */
    TimoshenkoLinearBeamElement();
    TimoshenkoLinearBeamElement( IndexType NewId, GeometryType::Pointer pGeometry);
    TimoshenkoLinearBeamElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /**
     * Destructor.
     */
    virtual ~TimoshenkoLinearBeamElement();

    /**
     * Operations.
     */

    virtual Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const;

    virtual Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom,
                            PropertiesType::Pointer pProperties) const;

    void Initialize(const ProcessInfo& rCurrentProcessInfo);

    void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                               VectorType& rRightHandSideVector,
                               const ProcessInfo& rCurrentProcessInfo);

    void CalculateMassMatrix( MatrixType& rMassMatrix,
                              const ProcessInfo& rCurrentProcessInfo);

    void CalculateDampingMatrix( MatrixType& rDampMatrix,
                                 const ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSide( VectorType& rRightHandSideVector,
                                 const ProcessInfo& rCurrentProcessInfo);

    void EquationIdVector( EquationIdVectorType& rResult,
                           const ProcessInfo& rCurrentProcessInfo) const;

    void GetDofList( DofsVectorType& ElementalDofList,
                     const ProcessInfo& CurrentProcessInfo) const;

    void CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& Output, const ProcessInfo& rCurrentProcessInfo );

    void CalculateOnIntegrationPoints( const Variable<array_1d<double, 3> >& rVariable, std::vector<array_1d<double, 3> >& Output, const ProcessInfo& rCurrentProcessInfo );

    int Check(const ProcessInfo& rCurrentProcessInfo) const;

    /**
     * Turn back information as a string.
     * (DEACTIVATED)
     */
    //std::string Info();

    /**
     * Print information about this object.
     * (DEACTIVATED)
     */
    //virtual void PrintInfo(std::ostream& rOStream) const;

    /**
     * Print object's data.
     * (DEACTIVATED)
     */
    //virtual void PrintData(std::ostream& rOStream) const;

protected:


private:

    // double mE, mG, mA, mNU, mI, mk;

    Matrix mInitialDisp;
    Matrix mInitialRot;

    friend class Serializer;

    virtual void save ( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, Element )
    }

    virtual void load ( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, Element )
    }

    void CalculateAll( MatrixType& rLeftHandSideMatrix,
                       VectorType& rRightHandSideVector,
                       const ProcessInfo& rCurrentProcessInfo,
                       const bool& CalculateStiffnessMatrixFlag,
                       const bool& CalculateResidualVectorFlag);

    void GetCurrentDisplacement( Vector& Disp ) const;

    void CalculateInitialLocalCS(Matrix& transformation_matrix) const;

    void CreateElementStiffnessMatrix_Material(Matrix& local_stiffness_matrix) const;

}; // Class TimoshenkoLinearBeamElement

}  // namespace Kratos.


#endif // KRATOS_TIMOSHENKO_LINEAR_BEAM_ELEMENT_H_INCLUDED defined


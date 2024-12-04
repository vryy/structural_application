//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 28 Mar 17 $
//   Revision:            $Revision: 1.0 $
//
//
#if !defined(KRATOS_DUMMY_ELEMENT_H_INCLUDED )
#define  KRATOS_DUMMY_ELEMENT_H_INCLUDED


// External includes

// Project includes
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

namespace Kratos
{

/**
 * A zero contribution element to support for assigning condition
 * It also supports for post processing at integration points using interpolation value from node
 */
class KRATOS_API(STRUCTURAL_APPLICATION) DummyElement : public Element
{
    public:
        // Counted pointer of DummyElement
        KRATOS_CLASS_POINTER_DEFINITION(DummyElement);

        /**
         * Default constructor.
         */
        DummyElement();
        DummyElement( IndexType NewId, GeometryType::Pointer pGeometry);
        DummyElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

        /**
         * Destructor.
         */
        virtual ~DummyElement();

        /**
         * Operations.
         */

        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                PropertiesType::Pointer pProperties) const final;

        Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom,
                                PropertiesType::Pointer pProperties) const final;

        void Initialize(const ProcessInfo& rCurrentProcessInfo) final;

        void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                   VectorType& rRightHandSideVector,
                                   const ProcessInfo& rCurrentProcessInfo) final;

        void CalculateRightHandSide( VectorType& rRightHandSideVector,
                                     const ProcessInfo& rCurrentProcessInfo) final;

        void EquationIdVector( EquationIdVectorType& rResult,
                               const ProcessInfo& rCurrentProcessInfo) const final;

        void GetDofList( DofsVectorType& ElementalDofList,
                         const ProcessInfo& CurrentProcessInfo) const final;

        void CalculateOnIntegrationPoints( const Variable<MatrixType>& rVariable, std::vector<MatrixType>& rValues, const ProcessInfo& rCurrentProcessInfo ) override;

        void CalculateOnIntegrationPoints( const Variable<VectorType>& rVariable, std::vector<VectorType>& rValues, const ProcessInfo& rCurrentProcessInfo ) override;

        void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3> >& rVariable, std::vector<array_1d<double, 3> >& rValues, const ProcessInfo& rCurrentProcessInfo) override;

        void CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo ) override;

        int Check(const ProcessInfo& rCurrentProcessInfo) const final;

        /**
         * Turn back information as a string.
         */
        std::string Info() const override
        {
            return "DummyElement";
        }

        /**
         * Print information about this object.
         */
        void PrintInfo(std::ostream& rOStream) const override
        {
            rOStream << Info();
        }

        /**
         * Print object's data.
         */
        void PrintData(std::ostream& rOStream) const override
        {
        }


    private:

        IntegrationMethod mThisIntegrationMethod;

        friend class Serializer;

        void save ( Serializer& rSerializer ) const final
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, Element )
        }

        void load ( Serializer& rSerializer ) final
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, Element )
        }

        void CalculateAll( MatrixType& rLeftHandSideMatrix,
                           VectorType& rRightHandSideVector,
                           const ProcessInfo& rCurrentProcessInfo,
                           bool CalculateStiffnessMatrixFlag,
                           bool CalculateResidualVectorFlag);

}; // Class DummyElement

}  // namespace Kratos.


#endif // KRATOS_DUMMY_ELEMENT_H_INCLUDED defined


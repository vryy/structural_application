/*
see license.txt
*/
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 6 Jul 2016 $
//   Revision:            $Revision: 0.0 $
//
#if !defined(KRATOS_EMBEDDED_NODE_PENALTY_TYING_CONDITION_H_INCLUDED )
#define  KRATOS_EMBEDDED_NODE_PENALTY_TYING_CONDITION_H_INCLUDED

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "utilities/math_utils.h"


namespace Kratos
{
    /**
     * Tying link from a node to an element using penalty method
     */
    class EmbeddedNodePenaltyTyingCondition : public Condition
    {
        public:
            typedef Condition                       BaseType;
            typedef BaseType::NodeType              NodeType;
            typedef NodeType::PointType             PointType;

            // Counted pointer of TipCondition
            KRATOS_CLASS_POINTER_DEFINITION(EmbeddedNodePenaltyTyingCondition);

            /**
             * Default constructor.
             */
            EmbeddedNodePenaltyTyingCondition( );
            EmbeddedNodePenaltyTyingCondition( IndexType NewId, GeometryType::Pointer pGeometry );
            EmbeddedNodePenaltyTyingCondition( IndexType NewId, GeometryType::Pointer pGeometry,
                        NodeType::Pointer& pSlaveNode,      // slave node that will be tied to the parent element
                        Element::Pointer& pParentElement,   // parent element
                        PointType& rLocalPoint );      // local coordinates of the slave node in the parent element

            /**
             * Destructor.
             */
            virtual ~EmbeddedNodePenaltyTyingCondition();

            /**
             * Operations.
             */
            Condition::Pointer Create( IndexType NewId, GeometryType::Pointer pGeometry,
                        NodeType::Pointer& pSlaveNode, Element::Pointer& pParentElement, PointType& rSolidLocalPoint ) const;

            void Initialize(const ProcessInfo& rCurrentProcessInfo);

            void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo);

            void CalculateRightHandSide( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo);

            void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                               const ProcessInfo& rCurrentProcessInfo,
                               bool CalculateStiffnessMatrixFlag, bool CalculateResidualVectorFlag);

            void EquationIdVector( EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const;

            void GetDofList( DofsVectorType& ConditionalDofList, const ProcessInfo& CurrentProcessInfo) const;

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

        private:

            friend class Serializer;

            virtual void save ( Serializer& rSerializer ) const
            {
                KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, Condition )
            }

            virtual void load ( Serializer& rSerializer )
            {
                KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, Condition )
            }

            NodeType::Pointer mpSlaveNode;
            Element::Pointer mpMasterElement;
            PointType mLocalPoint;

    }; // Class EmbeddedNodePenaltyTyingCondition
}  // namespace Kratos.

#endif // KRATOS_EMBEDDED_NODE_PENALTY_TYING_CONDITION_H_INCLUDED defined


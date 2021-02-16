/*
==============================================================================
KratosR1StructuralApplication
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
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: jelena $
//   Date:                $Date: 2013-10-16 $
//   Revision:            $Revision: 1.1 $
//
//
#if !defined(KRATOS_PILE_KINEMATIC_LINEAR_H_INCLUDED )
#define  KRATOS_PILE_KINEMATIC_LINEAR_H_INCLUDED


// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
// #include "utilities/math_utils.h"
#include "includes/serializer.h"

namespace Kratos
{
    /**
     * Contact link element for 3D contact problems.
     * This condition links two contact surfaces (one master and
     * one slave element) to a condition which may be assembled
     * at once by the builder.
     * Within this condition all system contribution with regard to
     * the gap function, the contact stress, the normal vectors
     * and their linearizations are calculated.
     */

    class Pile_Kinematic_Linear : public Condition
    {
        public:
            // Counted pointer of Pile_Kinematic_Linear
            KRATOS_CLASS_POINTER_DEFINITION( Pile_Kinematic_Linear );

            #ifdef SD_APP_FORWARD_COMPATIBILITY
            typedef Point PointType;
            #else
            typedef Point<3> PointType;
            #endif

            /**
             * Default constructor.
             */
            Pile_Kinematic_Linear();
            Pile_Kinematic_Linear( IndexType NewId, GeometryType::Pointer pGeometry );
            Pile_Kinematic_Linear( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

            Pile_Kinematic_Linear( IndexType NewId,
                                   GeometryType::Pointer pGeometry,
                                   PropertiesType::Pointer pProperties,
                                   Element::Pointer soilElement,
                                   Element::Pointer pileElement,
                                   PointType& rSoilLocalPoint,
                                   PointType& rPileLocalPoint,
                                   int SoilIntegrationPointIndex );

            /**
             * Destructor.
             */
            virtual ~Pile_Kinematic_Linear();

            /**
             * Operations.
             */

            void Initialize(const ProcessInfo& rCurrentProcessInfo);

            void InitializeSolutionStep(const ProcessInfo& CurrentProcessInfo);

            /**
             * Calculates the local system contributions for this contact element
             */
            void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                       VectorType& rRightHandSideVector,
                                       const ProcessInfo& rCurrentProcessInfo );

            void CalculateRightHandSide( VectorType& rRightHandSideVector,
                                         const ProcessInfo& rCurrentProcessInfo );

            void DampMatrix( MatrixType& rDampMatrix, const ProcessInfo& rCurrentProcessInfo );

            void EquationIdVector( EquationIdVectorType& rResult,
                                   const ProcessInfo& rCurrentProcessInfo ) const;

            void GetDofList( DofsVectorType& ConditionalDofList,
                             const ProcessInfo& CurrentProcessInfo ) const;

            /**
             * Turn back information as a string.
             * (DEACTIVATED)
             */
            //std::string Info();

            /**
             * Print information about this object.
             */
            virtual void PrintInfo( std::ostream& rOStream ) const;

            /**
             * Print object's data.
             */
            virtual void PrintData( std::ostream& rOStream ) const;

        protected:


        private:

            friend class Serializer;

            virtual void save( Serializer& rSerializer ) const
            {
                rSerializer.save( "name", "Pile_Kinematic_Linear" );
                KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition )
            }

            virtual void load( Serializer& rSerializer )
            {
                KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition )
            }

            Vector mvPile;
            Matrix mTPile;

            PointType mPileLocalPoint;
            PointType mSoilLocalPoint;
            Element::Pointer mpPileElement;
            Element::Pointer mpSoilElement;
            PointType mPileGlobalPoint;
            PointType mSoilGlobalPoint;
//            Vector mTPileGlobalVector;
            int mPileIntegrationPointIndex;

            void CalculateAll( MatrixType& rLeftHandSideMatrix,
                               VectorType& rRightHandSideVector,
                               const ProcessInfo& rCurrentProcessInfo,
                               bool CalculateStiffnessMatrixFlag,
                               bool CalculateResidualVectorFlag );

            void CalculateAndAdd_RHS( Vector& rRightHandSideVector,
                                      const Vector& NSoil,
                                      const Vector& NPile,
                                      const Vector& vPileNorm,
                                      const Matrix& T,
                                      const Vector& tangentialStresses,
                                      double Gap,
                                      double normalStress,
                                      double SoilIntegrationWeight,
                                      double dASoil,
                                      double friction_coeff );

            Vector GetRelativTangentialVelocity( Matrix& T );

            Vector GetRelativVelocity();

            void UpdatePileLocalPoint( );

            Vector NormalVector( Element::Pointer rElement, const GeometryType::CoordinatesArrayType& LocalPoint );

            PointType& GetGlobalCoordinates( Element::Pointer rElement, PointType& rResult, const PointType& LocalCoordinates );

            Matrix CalculateTangentVectors( const Element::Pointer rElement, const Matrix& DN );

            Matrix TangentialVectors( Element::Pointer rElement,
                                      const GeometryType::CoordinatesArrayType& LocalPoint );

            Matrix TangentialVectorsTotal( Element::Pointer rElement,
                                           const GeometryType::CoordinatesArrayType& LocalPoint );

            Matrix TangentialVectorsGlobal( Element::Pointer rElement,
                                            const GeometryType::CoordinatesArrayType& LocalPoint );

            Matrix TangentialVectors_inOrigin( Element::Pointer rElement,
                                               const GeometryType::CoordinatesArrayType& rPoint );
            //   Vector& GlobalCoordinatesTPileVector(Element::Pointer PileElements, Vector& rResult, PointType const& LocalCoordinates);
    }; // Class Pile_Kinematic_Linear
}  // namespace Kratos.


#endif // KRATOS_PILE_KINEMATIC_LINEAR_H_INCLUDED defined


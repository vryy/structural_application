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
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: nagel $
//   Date:                $Date: 2009-02-24 08:06:20 $
//   Revision:            $Revision: 1.2 $
//
//
#if !defined(KRATOS_CONTACT_LINK_3D_LAGRANGE_TYING_CONDITION_H_INCLUDED )
#define  KRATOS_CONTACT_LINK_3D_LAGRANGE_TYING_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

#include "custom_conditions/master_contact_face_3D.h"
#include "custom_conditions/slave_contact_face_3D.h"


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
    class Contact_Link_3D_Lagrange_Tying : public Condition
    {
    public:
        // Counted pointer of Contact_Link_3D_Lagrange_Tying
        KRATOS_CLASS_POINTER_DEFINITION(Contact_Link_3D_Lagrange_Tying);

        #ifdef SD_APP_FORWARD_COMPATIBILITY
        typedef Point PointType;
        #else
        typedef Point<3> PointType;
        #endif

        /**
         * Default constructor.
         */
        Contact_Link_3D_Lagrange_Tying( IndexType NewId, GeometryType::Pointer pGeometry);

        Contact_Link_3D_Lagrange_Tying( IndexType NewId, GeometryType::Pointer pGeometry,
                                       PropertiesType::Pointer pProperties
                                       );


        Contact_Link_3D_Lagrange_Tying( IndexType NewId, GeometryType::Pointer pGeometry,
                                       PropertiesType::Pointer pProperties,
                                       Condition::Pointer Master,
                                       Condition::Pointer Slave,
                                       PointType& MasterContactLocalPoint,
                                       PointType& SlaveContactLocalPoint,
                                       int SlaveIntegrationPointIndex
                                       );
        /**
         * Destructor.
         */
        virtual ~Contact_Link_3D_Lagrange_Tying();

        /**
         * Operations.
         */



        Condition::Pointer Create( IndexType NewId,
                                  NodesArrayType const& ThisNodes,
                                  PropertiesType::Pointer pProperties) const;

        Condition::Pointer Create( IndexType NewId,
                                  GeometryType::Pointer pGeom,
                                  PropertiesType::Pointer pProperties) const;

        void InitializeSolutionStep(const ProcessInfo& CurrentProcessInfo);
        /**
         * Calculates the local system contributions for this contact element
         */
        void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                  VectorType& rRightHandSideVector,
                                  const ProcessInfo& rCurrentProcessInfo);

        void CalculateRightHandSide( VectorType& rRightHandSideVector,
                                    const ProcessInfo& rCurrentProcessInfo);

        void CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo);

        void EquationIdVector( EquationIdVectorType& rResult,
                              const ProcessInfo& rCurrentProcessInfo) const;

        void GetDofList( DofsVectorType& ConditionalDofList,
                        const ProcessInfo& CurrentProcessInfo) const;

        /**
         * Turn back information as a string.
         * (DEACTIVATED)
         */
        //std::string Info();

        /**
         * Print information about this object.
         */
        virtual void PrintInfo(std::ostream& rOStream) const;

        /**
         * Print object's data.
         */
        virtual void PrintData(std::ostream& rOStream) const;

    protected:


    private:
        void CalculateAll( MatrixType& rLeftHandSideMatrix,
                          VectorType& rRightHandSideVector,
                          const ProcessInfo& rCurrentProcessInfo,
                          bool CalculateStiffnessMatrixFlag,
                          bool CalculateResidualVectorFlag);


        /////////////////////////
        /////
        //////////////////////////
        /*void CalculateAndAdd_RHS(
         Vector& residualvector,
         const Vector& NMaster,
         const Vector& NSlave,
         const Vector& vMaster,
         double normalStress,
         double SlaveIntegrationWeight,
         double dASlave ); */

        void CalculateAndAdd_RHS( Vector& residualvector,
                                 const Vector& NMaster,
                                 const Vector& NSlave,
                                 const Vector& vMaster,
                                 const Matrix& T,
                                 const Vector& tangentialStresses,
                                 double Gap,
                                 double normalStress,
                                 double SlaveIntegrationWeight,
                                 double dASlave);




        /**
         * This function calculates updates the local and global coordinates
         * of the master contact partner in order to follow the movement of
         * the slave surface along the master surface
         */
        void UpdateMasterLocalPoint( );


        Vector NormalVector( Condition::Pointer Surface,
                            const GeometryType::CoordinatesArrayType& LocalPoint );

        Matrix TangentialVectors( Condition::Pointer Surface,
                                 const GeometryType::CoordinatesArrayType& LocalPoint );

        Matrix TangentialVectorsGlobal( Condition::Pointer Surface,
                                       const GeometryType::CoordinatesArrayType& LocalPoint );

        Matrix TangentialVectors_inOrigin( Condition::Pointer Surface,
                                          const GeometryType::CoordinatesArrayType& rPoint );

        PointType& GlobalCoordinates(Condition::Pointer Surface, PointType& rResult, PointType const& LocalCoordinates);


        Vector GetRelativTangentialVelocity(Matrix& T);

        Vector GetRelativVelocity();
        /**
         * Assignment operator.
         * (DEACTIVATED)
         */
        //Contact_Link_3D_Lagrange_Tying& operator=(const Contact_Link_3D_Lagrange_Tying& rOther);

        /**
         * Copy constructor.
         * (DEACTIVATED)
         */
        //Contact_Link_3D_Lagrange_Tying(const Contact_Link_3D_Lagrange_Tying& rOther);


        /**
         * private members
         */

        ///@}
        ///@name Serialization
        ///@{
        friend class Serializer;

        // A private default constructor necessary for serialization
        Contact_Link_3D_Lagrange_Tying() {};

        virtual void save(Serializer& rSerializer) const
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
        }

        virtual void load(Serializer& rSerializer)
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
        }

        Vector mvMaster;
        Matrix mTMaster;
        //             Condition::Pointer mpSlave;
        //             Condition::Pointer mpMaster;
        //             PointType mMasterContactLocalPoint;
        //             PointType mSlaveContactLocalPoint;
        //             PointType mMasterContactGlobalPoint;
        //             PointType mSlaveContactGlobalPoint;
    }; // Class Contact_Link_3D_Lagrange_Tying
}  // namespace Kratos.

#endif // KRATOS_CONTACT_LINK_3D_LAGRANGE_TYING_CONDITION_H_INCLUDED  defined

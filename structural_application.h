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
//   Date:                $Date: 2009-03-20 08:55:34 $
//   Revision:            $Revision: 1.20 $
//
//


#if !defined(KRATOS_STRUCTURAL_APPLICATION_H_INCLUDED )
#define  KRATOS_STRUCTURAL_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/kratos_application.h"

#ifdef SD_APP_FORWARD_COMPATIBILITY
    #include "custom_elements/truss_element.h"
    #include "custom_elements/crisfield_truss_element.h"
    #include "custom_elements/timoshenko_beam_element.h"
    #include "custom_elements/timoshenko_linear_beam_element.h"
    #include "custom_elements/beam_element.h"
    #include "custom_elements/kinematic_linear.h"
    #include "custom_elements/kinematic_linear_axisymmetric.h"
    #include "custom_elements/total_lagrangian.h"
    #include "custom_conditions/point_force.h"
    #include "custom_conditions/pointmoment3D.h"
    #include "custom_conditions/face2D.h"
    #include "custom_conditions/face3D.h"
    #include "custom_conditions/face_pressure3D.h"
    #include "custom_conditions/face_pressure3D_total_lagrangian.h"
    #include "custom_conditions/face_traction3D.h"
    #include "custom_conditions/faceforce3D.h"
    #include "custom_conditions/line_force.h"
    #include "custom_conditions/line_pressure.h"
    #include "custom_conditions/line_traction.h"
#else
    #include "custom_elements/total_lagrangian.h"
    #include "custom_elements/total_lagrangian_axisymmetric.h"
    #include "custom_elements/updated_lagrangian.h"
    #include "custom_elements/finite_strain.h"
    #include "custom_elements/finite_strain_axisymmetric.h"
    #include "custom_elements/beam_element.h"
    #include "custom_elements/timoshenko_beam_element.h"
    #include "custom_elements/timoshenko_linear_beam_element.h"
    #include "custom_elements/corotational_linear_beam_element.h"
    #include "custom_elements/kinematic_linear.h"
    #include "custom_elements/kinematic_linear_axisymmetric.h"
    #include "custom_elements/kinematic_linear_anti_plane.h"
    #include "custom_elements/updated_kinematic_linear.h"
    #include "custom_elements/shell_isotropic.h"
    #include "custom_elements/shell_anisotropic.h"
    #include "custom_elements/crisfield_truss_element.h"
    #include "custom_elements/truss_element.h"
    #include "custom_elements/ebst.h"
    #include "custom_elements/ebst_vel.h"
    #include "custom_elements/eas_element_q4e4.h"
    #include "custom_elements/dummy_element.h"

    #include "custom_conditions/face2D.h"
    #include "custom_conditions/face3D.h"
    #include "custom_conditions/face_pressure3D.h"
    #include "custom_conditions/face_pressure3D_total_lagrangian.h"
    #include "custom_conditions/face_traction3D.h"
    #include "custom_conditions/faceforce3D.h"
    #include "custom_conditions/line_force.h"
    #include "custom_conditions/line_pressure.h"
    #include "custom_conditions/line_pressure_distributed.h"
    #include "custom_conditions/line_traction.h"
    #include "custom_conditions/point_force.h"
    #include "custom_conditions/point_point_joint_condition.h"
    #include "custom_conditions/point_point_lagrange_condition.h"
    #include "custom_conditions/elastic_constraint.h"
    #include "custom_conditions/elastic_line_springs.h"
    #include "custom_conditions/elastic_face_springs.h"
    #include "custom_conditions/nitsche_isotropic_constraint.h"
    #include "custom_conditions/roller_constraint.h"
    #include "custom_conditions/mean_displacement_constraint.h"
    #include "custom_conditions/dummy_condition.h"

    #include "constitutive_laws/isotropic_3d.h"
    #include "constitutive_laws/dummy_constitutive_law.h"
    #include "constitutive_laws/cam_clay_3d.h"
#endif

#define STRUCTURAL_APPLICATION_DEFINE_ELEMENT_ALL_GEOMETRIES(element_type) \
    const element_type m##element_type##2D3N; \
    const element_type m##element_type##2D6N; \
    const element_type m##element_type##2D4N; \
    const element_type m##element_type##2D8N; \
    const element_type m##element_type##2D9N; \
    const element_type m##element_type##3D4N; \
    const element_type m##element_type##3D10N; \
    const element_type m##element_type##3D8N; \
    const element_type m##element_type##3D20N; \
    const element_type m##element_type##3D27N; \
    const element_type m##element_type##3D6N; \
    const element_type m##element_type##3D15N;

#define STRUCTURAL_APPLICATION_DEFINE_ELEMENT_ALL_2D_GEOMETRIES(element_type) \
    const element_type m##element_type##3N; \
    const element_type m##element_type##6N; \
    const element_type m##element_type##4N; \
    const element_type m##element_type##8N; \
    const element_type m##element_type##9N; \

#define STRUCTURAL_APPLICATION_DEFINE_CONDITION_ALL_GEOMETRIES(condition_type) \
    const condition_type m##condition_type##2D2N; \
    const condition_type m##condition_type##2D3N; \
    const condition_type m##condition_type##3D3N; \
    const condition_type m##condition_type##3D6N; \
    const condition_type m##condition_type##3D4N; \
    const condition_type m##condition_type##3D8N; \
    const condition_type m##condition_type##3D9N;

namespace Kratos
{
///@name Kratos Globals
///@{
///@}
///@name Type Definitions
///@{
///@}
///@name  Enum's
///@{
///@}
///@name  Functions
///@{
///@}
///@name Kratos Classes
///@{

/// Structural Application for KRATOS.
/**
 * This application features Elements, Conditions, Constitutive laws and Utilities
 * for structural analysis problems
 */
class KRATOS_API(STRUCTURAL_APPLICATION) KratosStructuralApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosStructuralApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosStructuralApplication);
    ///@}
    ///@name Life Cycle
    ///@{
    /// Default constructor.
    KratosStructuralApplication();

    /// Destructor.
    ~KratosStructuralApplication() override {}

    ///@}
    ///@name Operators
    ///@{
    ///@}
    ///@name Operations
    ///@{

    /**
     * Registers the structural application in the KRATOS kernel
     */
    void Register() final;

    /**
     * Registers the structural application variables in the KRATOS kernel
     */
    void RegisterVariables() final;

    ///@}
    ///@name Access
    ///@{
    ///@}
    ///@name Inquiry
    ///@{

    /// A size map for vector variable which can be used for visualization
    static std::map<Variable<Vector>, std::size_t> StandardVectorVariableSizeMap();

    ///@}
    ///@name Input and output
    ///@{
    /// Turn back information as a string.
    std::string Info() const final
    {
        return "KratosStructuralApplication";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const final
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const final
    {
        KRATOS_WATCH("KratosStructuralApplication application contains following components:");
        KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );
        rOStream << "Variables:" << std::endl;
        KratosComponents<VariableData>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Elements:" << std::endl;
        KratosComponents<Element>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Conditions:" << std::endl;
        KratosComponents<Condition>().PrintData(rOStream);
    }
    ///@}
    ///@name Friends
    ///@{
    ///@}

private:
    ///@name Member Variables
    ///@{

#ifdef SD_APP_FORWARD_COMPATIBILITY

    const CrisfieldTrussElement mCrisfieldTrussElement3D2N;
    const CrisfieldTrussElement mCrisfieldTrussElement3D3N;
    const TrussElement mTrussElement3D2N;
    const TrussElement mTrussElement3D3N;
    const BeamElement mBeamElement3D2N;
    const BeamElement mBeamElement3D3N;
    const TimoshenkoBeamElement mTimoshenkoBeamElement3D2N;
    const TimoshenkoBeamElement mTimoshenkoBeamElement3D3N;
    const TimoshenkoLinearBeamElement mTimoshenkoLinearBeamElement2D2N;
    const TimoshenkoLinearBeamElement mTimoshenkoLinearBeamElement3D2N;

    const TotalLagrangian mTotalLagrangian2D3N;
    const TotalLagrangian mTotalLagrangian2D4N;
    const TotalLagrangian mTotalLagrangian2D6N;
    const TotalLagrangian mTotalLagrangian2D8N;
    const TotalLagrangian mTotalLagrangian2D9N;
    const TotalLagrangian mTotalLagrangian3D4N;
    const TotalLagrangian mTotalLagrangian3D10N;
    const TotalLagrangian mTotalLagrangian3D6N;
    const TotalLagrangian mTotalLagrangian3D15N;
    const TotalLagrangian mTotalLagrangian3D8N;
    const TotalLagrangian mTotalLagrangian3D20N;
    const TotalLagrangian mTotalLagrangian3D27N;

    const KinematicLinear mKinematicLinear2D3N;
    const KinematicLinear mKinematicLinear2D4N;
    const KinematicLinear mKinematicLinear2D6N;
    const KinematicLinear mKinematicLinear2D8N;
    const KinematicLinear mKinematicLinear2D9N;
    const KinematicLinear mKinematicLinear3D4N;
    const KinematicLinear mKinematicLinear3D10N;
    const KinematicLinear mKinematicLinear3D8N;
    const KinematicLinear mKinematicLinear3D20N;
    const KinematicLinear mKinematicLinear3D27N;
    const KinematicLinear mKinematicLinear3D6N;
    const KinematicLinear mKinematicLinear3D15N;

    const KinematicLinearAxisymmetric mKinematicLinearAxisymmetric3N;
    const KinematicLinearAxisymmetric mKinematicLinearAxisymmetric4N;
    const KinematicLinearAxisymmetric mKinematicLinearAxisymmetric6N;
    const KinematicLinearAxisymmetric mKinematicLinearAxisymmetric8N;
    const KinematicLinearAxisymmetric mKinematicLinearAxisymmetric9N;

    const PointForce2D mPointForce2D;
    const PointForce3D mPointForce3D;
    const Face2D  mFace2D;
    const Face3D  mFace3D3N;
    const Face3D  mFace3D6N;
    const Face3D  mFace3D4N;
    const Face3D  mFace3D8N;
    const Face3D  mFace3D9N;
    const FacePressure3D  mFacePressure3D3N;
    const FacePressure3D  mFacePressure3D6N;
    const FacePressure3D  mFacePressure3D4N;
    const FacePressure3D  mFacePressure3D8N;
    const FacePressure3D  mFacePressure3D9N;
    const FacePressure3DTotalLagrangian  mFacePressureTotalLagrangian3D3N;
    const FacePressure3DTotalLagrangian  mFacePressureTotalLagrangian3D6N;
    const FacePressure3DTotalLagrangian  mFacePressureTotalLagrangian3D4N;
    const FacePressure3DTotalLagrangian  mFacePressureTotalLagrangian3D8N;
    const FacePressure3DTotalLagrangian  mFacePressureTotalLagrangian3D9N;
    const FaceTraction3D  mFaceTraction3D3N;
    const FaceTraction3D  mFaceTraction3D6N;
    const FaceTraction3D  mFaceTraction3D4N;
    const FaceTraction3D  mFaceTraction3D8N;
    const FaceTraction3D  mFaceTraction3D9N;
    const LineForce mLineForce2D2N;
    const LineForce mLineForce2D3N;
    const LineForce mLineForce3D2N;
    const LineForce mLineForce3D3N;
    const LinePressure mLinePressure2D2N;
    const LinePressure mLinePressure2D3N;
    const LineTraction mLineTraction2D2N;
    const LineTraction mLineTraction2D3N;
    const FaceForce3D mFaceForce3D3N;
    const FaceForce3D mFaceForce3D6N;
    const FaceForce3D mFaceForce3D4N;
    const FaceForce3D mFaceForce3D8N;
    const FaceForce3D mFaceForce3D9N;

#else // SD_APP_FORWARD_COMPATIBILITY

    const CrisfieldTrussElement mCrisfieldTrussElement3D2N;
    const CrisfieldTrussElement mCrisfieldTrussElement3D3N;
    const TrussElement mTrussElement3D2N;
    const TrussElement mTrussElement3D3N;
    const BeamElement mBeamElement3D2N;
    const BeamElement mBeamElement3D3N;
    const TimoshenkoBeamElement mTimoshenkoBeamElement3D2N;
    const TimoshenkoBeamElement mTimoshenkoBeamElement3D3N;
    const TimoshenkoLinearBeamElement mTimoshenkoLinearBeamElement2D2N;
    const TimoshenkoLinearBeamElement mTimoshenkoLinearBeamElement3D2N;
    const CorotationalLinearBeamElement mCorotationalLinearBeamElement2D2N;
    const CorotationalLinearBeamElement mCorotationalLinearBeamElement3D2N;
    const ShellIsotropic mIsoShellElement;
    const ShellAnisotropic mAnisoShellElement;

    STRUCTURAL_APPLICATION_DEFINE_ELEMENT_ALL_GEOMETRIES(TotalLagrangian)

    STRUCTURAL_APPLICATION_DEFINE_ELEMENT_ALL_2D_GEOMETRIES(TotalLagrangianAxisymmetric)

    STRUCTURAL_APPLICATION_DEFINE_ELEMENT_ALL_GEOMETRIES(UpdatedLagrangian)

    STRUCTURAL_APPLICATION_DEFINE_ELEMENT_ALL_GEOMETRIES(FiniteStrain)

    STRUCTURAL_APPLICATION_DEFINE_ELEMENT_ALL_2D_GEOMETRIES(FiniteStrainAxisymmetric)

    STRUCTURAL_APPLICATION_DEFINE_ELEMENT_ALL_GEOMETRIES(KinematicLinear)
    STRUCTURAL_APPLICATION_DEFINE_ELEMENT_ALL_GEOMETRIES(ComplexKinematicLinear)
    STRUCTURAL_APPLICATION_DEFINE_ELEMENT_ALL_GEOMETRIES(GComplexKinematicLinear)

    STRUCTURAL_APPLICATION_DEFINE_ELEMENT_ALL_2D_GEOMETRIES(KinematicLinearAxisymmetric)

    STRUCTURAL_APPLICATION_DEFINE_ELEMENT_ALL_2D_GEOMETRIES(KinematicLinearAntiPlane)
    STRUCTURAL_APPLICATION_DEFINE_ELEMENT_ALL_2D_GEOMETRIES(ComplexKinematicLinearAntiPlane)
    STRUCTURAL_APPLICATION_DEFINE_ELEMENT_ALL_2D_GEOMETRIES(GComplexKinematicLinearAntiPlane)

    STRUCTURAL_APPLICATION_DEFINE_ELEMENT_ALL_GEOMETRIES(UpdatedKinematicLinear)

    const Ebst mEbst3D3N;
    const EbstVel mEbstVel3D3N;
    const EASElementQ4E4 mEASElementQ4E4;

    const DummyElement mDummySurfaceElement2D3N;
    const DummyElement mDummySurfaceElement2D6N;
    const DummyElement mDummySurfaceElement2D4N;
    const DummyElement mDummySurfaceElement2D8N;
    const DummyElement mDummySurfaceElement2D9N;

    const DummyElement mDummySurfaceElement3D3N;
    const DummyElement mDummySurfaceElement3D6N;
    const DummyElement mDummySurfaceElement3D4N;
    const DummyElement mDummySurfaceElement3D8N;
    const DummyElement mDummySurfaceElement3D9N;

    const DummyElement mDummyVolumeElement3D4N;
    const DummyElement mDummyVolumeElement3D10N;
    const DummyElement mDummyVolumeElement3D8N;
    const DummyElement mDummyVolumeElement3D20N;
    const DummyElement mDummyVolumeElement3D27N;

    STRUCTURAL_APPLICATION_DEFINE_ELEMENT_ALL_GEOMETRIES(DummyElement)

    const Face2D  mFace2D;
    const Face3D  mFace3D3N;
    const Face3D  mFace3D6N;
    const Face3D  mFace3D4N;
    const Face3D  mFace3D8N;
    const Face3D  mFace3D9N;
    const FacePressure3D  mFacePressure3D3N;
    const FacePressure3D  mFacePressure3D6N;
    const FacePressure3D  mFacePressure3D4N;
    const FacePressure3D  mFacePressure3D8N;
    const FacePressure3D  mFacePressure3D9N;
    const FacePressure3DTotalLagrangian  mFacePressureTotalLagrangian3D3N;
    const FacePressure3DTotalLagrangian  mFacePressureTotalLagrangian3D6N;
    const FacePressure3DTotalLagrangian  mFacePressureTotalLagrangian3D4N;
    const FacePressure3DTotalLagrangian  mFacePressureTotalLagrangian3D8N;
    const FacePressure3DTotalLagrangian  mFacePressureTotalLagrangian3D9N;
    const FaceTraction3D  mFaceTraction3D3N;
    const FaceTraction3D  mFaceTraction3D6N;
    const FaceTraction3D  mFaceTraction3D4N;
    const FaceTraction3D  mFaceTraction3D8N;
    const FaceTraction3D  mFaceTraction3D9N;
    const LineForce mLineForce2D2N;
    const LineForce mLineForce2D3N;
    const LineForce mLineForce3D2N;
    const LineForce mLineForce3D3N;
    const LinePressure mLinePressure2D2N;
    const LinePressure mLinePressure2D3N;
    const LinePressureDistributed mLinePressureDistributed2D2N;
    const LinePressureDistributed mLinePressureDistributed2D3N;
    const LineTraction mLineTraction2D2N;
    const LineTraction mLineTraction2D3N;
    const FaceForce3D mFaceForce3D3N;
    const FaceForce3D mFaceForce3D6N;
    const FaceForce3D mFaceForce3D4N;
    const FaceForce3D mFaceForce3D8N;
    const FaceForce3D mFaceForce3D9N;
    const PointForce3D mPointForce3D;
    const ComplexPointForce3D mComplexPointForce3D;
    const PointForce2D mPointForce2D;
    const ComplexPointForce2D mComplexPointForce2D;
    const ElasticConstraint mElasticPointConstraint;
    const ElasticConstraint mElasticLineConstraint2N;
    const ElasticConstraint mElasticLineConstraint3N;
    const ElasticConstraint mElasticFaceConstraint3N;
    const ElasticConstraint mElasticFaceConstraint6N;
    const ElasticConstraint mElasticFaceConstraint4N;
    const ElasticConstraint mElasticFaceConstraint8N;
    const ElasticConstraint mElasticFaceConstraint9N;
    const ElasticLineSprings mElasticLineSprings2N;
    const ElasticLineSprings mElasticLineSprings3N;
    const ComplexElasticLineSprings mComplexElasticLineSprings2N;
    const ComplexElasticLineSprings mComplexElasticLineSprings3N;
    const ElasticFaceSprings mElasticFaceSprings3N;
    const ElasticFaceSprings mElasticFaceSprings6N;
    const ElasticFaceSprings mElasticFaceSprings4N;
    const ElasticFaceSprings mElasticFaceSprings8N;
    const ElasticFaceSprings mElasticFaceSprings9N;
    const ComplexElasticFaceSprings mComplexElasticFaceSprings3N;
    const ComplexElasticFaceSprings mComplexElasticFaceSprings6N;
    const ComplexElasticFaceSprings mComplexElasticFaceSprings4N;
    const ComplexElasticFaceSprings mComplexElasticFaceSprings8N;
    const ComplexElasticFaceSprings mComplexElasticFaceSprings9N;
    const PointPointJointCondition mPointPointJointCondition;
    const PointPointLagrangeCondition mPointPointLagrangeCondition;
    const NitscheIsotropicConstraint mNitscheIsotropicConstraint2D2N;
    const NitscheIsotropicConstraint mNitscheIsotropicConstraint2D3N;
    const NitscheIsotropicConstraint mNitscheIsotropicConstraint3D3N;
    const NitscheIsotropicConstraint mNitscheIsotropicConstraint3D6N;
    const NitscheIsotropicConstraint mNitscheIsotropicConstraint3D4N;
    const NitscheIsotropicConstraint mNitscheIsotropicConstraint3D8N;
    const NitscheIsotropicConstraint mNitscheIsotropicConstraint3D9N;
    const RollerConstraint mRollerConstraint2D2N;
    const RollerConstraint mRollerConstraint2D3N;
    const RollerConstraint mRollerConstraint3D3N;
    const RollerConstraint mRollerConstraint3D6N;
    const RollerConstraint mRollerConstraint3D4N;
    const RollerConstraint mRollerConstraint3D8N;
    const RollerConstraint mRollerConstraint3D9N;

    typedef MeanDisplacementConstraint<0> MeanDisplacementConstraintX;
    typedef MeanDisplacementConstraint<1> MeanDisplacementConstraintY;
    typedef MeanDisplacementConstraint<2> MeanDisplacementConstraintZ;
    STRUCTURAL_APPLICATION_DEFINE_CONDITION_ALL_GEOMETRIES(MeanDisplacementConstraintX)
    STRUCTURAL_APPLICATION_DEFINE_CONDITION_ALL_GEOMETRIES(MeanDisplacementConstraintY)
    STRUCTURAL_APPLICATION_DEFINE_CONDITION_ALL_GEOMETRIES(MeanDisplacementConstraintZ)

    const DummyCondition mDummyLineCondition2D2N;
    const DummyCondition mDummyLineCondition2D3N;
    const DummyCondition mDummySurfaceCondition2D3N;
    const DummyCondition mDummySurfaceCondition2D6N;
    const DummyCondition mDummySurfaceCondition2D4N;
    const DummyCondition mDummySurfaceCondition2D8N;
    const DummyCondition mDummySurfaceCondition2D9N;

    const DummyCondition mDummySurfaceCondition3D3N;
    const DummyCondition mDummySurfaceCondition3D6N;
    const DummyCondition mDummySurfaceCondition3D4N;
    const DummyCondition mDummySurfaceCondition3D8N;
    const DummyCondition mDummySurfaceCondition3D9N;

    const DummyCondition mDummyConditionPoint2D;
    const DummyCondition mDummyConditionPoint3D;
    const DummyCondition mDummyConditionLine2N;
    const DummyCondition mDummyConditionLine3N;
    const DummyCondition mDummyCondition2D3N;
    const DummyCondition mDummyCondition2D4N;
    const DummyCondition mDummyCondition2D6N;
    const DummyCondition mDummyCondition2D8N;
    const DummyCondition mDummyCondition2D9N;
    const DummyCondition mDummyCondition3D4N;
    const DummyCondition mDummyCondition3D10N;
    const DummyCondition mDummyCondition3D8N;
    const DummyCondition mDummyCondition3D20N;
    const DummyCondition mDummyCondition3D27N;
    const DummyCondition mDummyCondition3D6N;
    const DummyCondition mDummyCondition3D15N;

    const Isotropic3D mIsotropic3D;
    const DummyConstitutiveLaw mDummyConstitutiveLaw;
    const CamClay3D mCamClay3D;
#endif

//             const UPCTestElement mUPCTestElement3D20N;
    ///@}
    ///@name Private Operators
    ///@{
    ///@}
    ///@name Private Operations
    ///@{
    ///@}
    ///@name Private  Access
    ///@{
    ///@}
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    KratosStructuralApplication& operator=(KratosStructuralApplication const& rOther);

    /// Copy constructor.
    KratosStructuralApplication(KratosStructuralApplication const& rOther);

    ///@}
}; // Class KratosStructuralApplication
///@}
}  // namespace Kratos.
#endif // KRATOS_STRUCTURAL_APPLICATION_H_INCLUDED  defined

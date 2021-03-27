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
//   Date:                $Date: 2009-03-20 08:55:31 $
//   Revision:            $Revision: 1.23 $
//
//



// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/variables.h"

#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/quadrilateral_2d_9.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/prism_3d_6.h"
#include "geometries/prism_3d_15.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_3d_2.h"
#include "geometries/line_3d_3.h"
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"

#include "structural_application.h"
#include "structural_application_variables.h"


namespace Kratos
{

KratosStructuralApplication::KratosStructuralApplication()
#ifdef SD_APP_FORWARD_COMPATIBILITY
    : KratosApplication("StructuralApplication")
    , mCrisfieldTrussElement3D2N( 0, Element::GeometryType::Pointer( new Line3D2<Node<3> >( Element::GeometryType::PointsArrayType( 2 ) ) ) )
    , mCrisfieldTrussElement3D3N( 0, Element::GeometryType::Pointer( new Line3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) )
    , mTrussElement3D2N( 0, Element::GeometryType::Pointer( new Line3D2<Node<3> >( Element::GeometryType::PointsArrayType( 2 ) ) ) )
    , mTrussElement3D3N( 0, Element::GeometryType::Pointer( new Line3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) )
    , mBeamElement3D2N( 0, Element::GeometryType::Pointer( new Line3D2 <Node<3> >( Element::GeometryType::PointsArrayType( 2 ) ) ) )
    , mBeamElement3D3N( 0, Element::GeometryType::Pointer( new Line3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) )
    , mTimoshenkoBeamElement3D2N( 0, Element::GeometryType::Pointer( new Line3D2 <Node<3> >( Element::GeometryType::PointsArrayType( 2 ) ) ) )
    , mTimoshenkoBeamElement3D3N( 0, Element::GeometryType::Pointer( new Line3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) )
    , mTotalLagrangian2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) )
    , mTotalLagrangian2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) )
    , mTotalLagrangian2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6 ) ) ) )
    , mTotalLagrangian2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) )
    , mTotalLagrangian3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) )
    , mTotalLagrangian3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <Node<3> >( Element::GeometryType::PointsArrayType( 10 ) ) ) )
    , mTotalLagrangian3D6N( 0, Element::GeometryType::Pointer( new Prism3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6 ) ) ) )
    , mTotalLagrangian3D15N( 0, Element::GeometryType::Pointer( new Prism3D15 <Node<3> >( Element::GeometryType::PointsArrayType( 15 ) ) ) )
    , mTotalLagrangian3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) )
    , mTotalLagrangian3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <Node<3> >( Element::GeometryType::PointsArrayType( 20 ) ) ) )
    , mTotalLagrangian3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <Node<3> >( Element::GeometryType::PointsArrayType( 27 ) ) ) )
    , mKinematicLinear2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) )
    , mKinematicLinear2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) )
    , mKinematicLinear2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6 ) ) ) )
    , mKinematicLinear2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) )
    , mKinematicLinear2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9 ) ) ) )
    , mKinematicLinear3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) )
    , mKinematicLinear3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <Node<3> >( Element::GeometryType::PointsArrayType( 10 ) ) ) )
    , mKinematicLinear3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) )
    , mKinematicLinear3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <Node<3> >( Element::GeometryType::PointsArrayType( 20 ) ) ) )
    , mKinematicLinear3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <Node<3> >( Element::GeometryType::PointsArrayType( 27 ) ) ) )
    , mKinematicLinear3D6N( 0, Element::GeometryType::Pointer( new Prism3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6 ) ) ) )
    , mKinematicLinear3D15N( 0, Element::GeometryType::Pointer( new Prism3D15 <Node<3> >( Element::GeometryType::PointsArrayType( 15 ) ) ) )
#else
    : KratosApplication()
    , mCrisfieldTrussElement3D2N( 0, Element::GeometryType::Pointer( new Line3D2<Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) )
    , mCrisfieldTrussElement3D3N( 0, Element::GeometryType::Pointer( new Line3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mTrussElement3D2N( 0, Element::GeometryType::Pointer( new Line3D2<Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) )
    , mTrussElement3D3N( 0, Element::GeometryType::Pointer( new Line3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mBeamElement3D2N( 0, Element::GeometryType::Pointer( new Line3D2 <Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) )
    , mBeamElement3D3N( 0, Element::GeometryType::Pointer( new Line3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mTimoshenkoBeamElement3D2N( 0, Element::GeometryType::Pointer( new Line3D2 <Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) )
    , mTimoshenkoBeamElement3D3N( 0, Element::GeometryType::Pointer( new Line3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mCorotationalLinearBeamElement2D2N( 0, Element::GeometryType::Pointer( new Line2D2 <Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) )
    , mCorotationalLinearBeamElement3D2N( 0, Element::GeometryType::Pointer( new Line3D2 <Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) )
    , mTimoshenkoLinearBeamElement2D2N( 0, Element::GeometryType::Pointer( new Line2D2 <Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) )
    , mTimoshenkoLinearBeamElement3D2N( 0, Element::GeometryType::Pointer( new Line3D2 <Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) )
    , mIsoShellElement( 0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mAnisoShellElement( 0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    // mAnisoLinearShellElement( 0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mMembraneElement( 0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mTotalLagrangian2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mTotalLagrangian2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mTotalLagrangian2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) )
    , mTotalLagrangian2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mTotalLagrangian3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mTotalLagrangian3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <Node<3> >( Element::GeometryType::PointsArrayType( 10, Node<3>() ) ) ) )
    , mTotalLagrangian3D6N( 0, Element::GeometryType::Pointer( new Prism3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) )
    , mTotalLagrangian3D15N( 0, Element::GeometryType::Pointer( new Prism3D15 <Node<3> >( Element::GeometryType::PointsArrayType( 15, Node<3>() ) ) ) )
    , mTotalLagrangian3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mTotalLagrangian3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <Node<3> >( Element::GeometryType::PointsArrayType( 20, Node<3>() ) ) ) )
    , mTotalLagrangian3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <Node<3> >( Element::GeometryType::PointsArrayType( 27, Node<3>() ) ) ) )

    //, mLinearIncompresibleElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    //, mLinearIncompresibleElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4<Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )


    , mMixedLagrangian2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mMixedLagrangian2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mMixedLagrangian2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) )
    , mMixedLagrangian2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mMixedLagrangian3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mMixedLagrangian3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <Node<3> >( Element::GeometryType::PointsArrayType( 10, Node<3>() ) ) ) )
    , mMixedLagrangian3D6N( 0, Element::GeometryType::Pointer( new Prism3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) )
    , mMixedLagrangian3D15N( 0, Element::GeometryType::Pointer( new Prism3D15 <Node<3> >( Element::GeometryType::PointsArrayType( 15, Node<3>() ) ) ) )
    , mMixedLagrangian3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mMixedLagrangian3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <Node<3> >( Element::GeometryType::PointsArrayType( 20, Node<3>() ) ) ) )
    , mMixedLagrangian3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <Node<3> >( Element::GeometryType::PointsArrayType( 27, Node<3>() ) ) ) )

    , mKinematicLinear2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mKinematicLinear2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mKinematicLinear2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) )
    , mKinematicLinear2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mKinematicLinear2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) )
    , mKinematicLinear3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mKinematicLinear3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <Node<3> >( Element::GeometryType::PointsArrayType( 10, Node<3>() ) ) ) )
    , mKinematicLinear3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mKinematicLinear3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <Node<3> >( Element::GeometryType::PointsArrayType( 20, Node<3>() ) ) ) )
    , mKinematicLinear3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <Node<3> >( Element::GeometryType::PointsArrayType( 27, Node<3>() ) ) ) )
    , mKinematicLinear3D6N( 0, Element::GeometryType::Pointer( new Prism3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) )
    , mKinematicLinear3D15N( 0, Element::GeometryType::Pointer( new Prism3D15 <Node<3> >( Element::GeometryType::PointsArrayType( 15, Node<3>() ) ) ) )
    , mUnsaturatedSoilsElement2PhaseSmallStrain3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4<Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mUnsaturatedSoilsElement2PhaseSmallStrain3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mUnsaturatedSoilsElement2PhaseSmallStrain3D6N( 0, Element::GeometryType::Pointer( new Prism3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) )
    , mUnsaturatedSoilsElement2PhaseSmallStrain3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10<Node<3> >( Element::GeometryType::PointsArrayType( 10, Node<3>() ) ) ) )
    , mUnsaturatedSoilsElement2PhaseSmallStrain3D15N( 0, Element::GeometryType::Pointer( new Prism3D15<Node<3> >( Element::GeometryType::PointsArrayType( 15, Node<3>() ) ) ) )
    , mUnsaturatedSoilsElement2PhaseSmallStrain3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20<Node<3> >( Element::GeometryType::PointsArrayType( 20, Node<3>() ) ) ) )
    , mUnsaturatedSoilsElement2PhaseSmallStrain3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27<Node<3> >( Element::GeometryType::PointsArrayType( 27, Node<3>() ) ) ) )
    , mUnsaturatedSoilsElement2PhaseSmallStrainStaggered3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27<Node<3> >( Element::GeometryType::PointsArrayType( 27, Node<3>() ) ) ) )
    , mUnsaturatedSoilsElement3PhaseSmallStrain3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4<Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mUnsaturatedSoilsElement3PhaseSmallStrain3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mUnsaturatedSoilsElement3PhaseSmallStrain3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10<Node<3> >( Element::GeometryType::PointsArrayType( 10, Node<3>() ) ) ) )
    , mUnsaturatedSoilsElement3PhaseSmallStrain3D15N( 0, Element::GeometryType::Pointer( new Prism3D15<Node<3> >( Element::GeometryType::PointsArrayType( 15, Node<3>() ) ) ) )
    , mUnsaturatedSoilsElement3PhaseSmallStrain3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20<Node<3> >( Element::GeometryType::PointsArrayType( 20, Node<3>() ) ) ) )
    , mUnsaturatedSoilsElement3PhaseSmallStrain3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27<Node<3> >( Element::GeometryType::PointsArrayType( 27, Node<3>() ) ) ) )
    , mUnsaturatedSoilsElement3PhaseSmallStrainLiakopolous3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20<Node<3> >( Element::GeometryType::PointsArrayType( 20, Node<3>() ) ) ) )
    , mUnsaturatedSoilsElement3PhaseSmallStrainLiakopolous3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27<Node<3> >( Element::GeometryType::PointsArrayType( 27, Node<3>() ) ) ) )
    , mEbst3D3N( 0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mEbstVel3D3N( 0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mEASElementQ4E4( 0, Element::GeometryType::Pointer( new Quadrilateral2D4<Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )

    , mFace2D( 0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) )
    , mFace3D3N( 0, Element::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mFace3D6N( 0, Element::GeometryType::Pointer( new Triangle3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) )
    , mFace3D4N( 0, Element::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mFace3D8N( 0, Element::GeometryType::Pointer( new Quadrilateral3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mFace3D9N( 0, Element::GeometryType::Pointer( new Quadrilateral3D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) )
    , mFacePressure3D3N( 0, Element::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mFacePressure3D6N( 0, Element::GeometryType::Pointer( new Triangle3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) )
    , mFacePressure3D4N( 0, Element::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mFacePressure3D8N( 0, Element::GeometryType::Pointer( new Quadrilateral3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mFacePressure3D9N( 0, Element::GeometryType::Pointer( new Quadrilateral3D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) )
    , mFacePressureTotalLagrangian3D3N( 0, Element::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mFacePressureTotalLagrangian3D6N( 0, Element::GeometryType::Pointer( new Triangle3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) )
    , mFacePressureTotalLagrangian3D4N( 0, Element::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mFacePressureTotalLagrangian3D8N( 0, Element::GeometryType::Pointer( new Quadrilateral3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mFacePressureTotalLagrangian3D9N( 0, Element::GeometryType::Pointer( new Quadrilateral3D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) )
    , mLineForce2D2N( 0, Element::GeometryType::Pointer( new Line2D2 <Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) )
    , mLineForce2D3N( 0, Element::GeometryType::Pointer( new Line2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mLineForce3D2N( 0, Element::GeometryType::Pointer( new Line3D2 <Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) )
    , mLineForce3D3N( 0, Element::GeometryType::Pointer( new Line3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mLinePressure2D2N( 0, Element::GeometryType::Pointer( new Line2D2 <Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) )
    , mLinePressure2D3N( 0, Element::GeometryType::Pointer( new Line2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mFaceForce3D3N( 0, Element::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mFaceForce3D6N( 0, Element::GeometryType::Pointer( new Triangle3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) )
    , mFaceForce3D4N( 0, Element::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mFaceForce3D8N( 0, Element::GeometryType::Pointer( new Quadrilateral3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mFaceForce3D9N( 0, Element::GeometryType::Pointer( new Quadrilateral3D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) )
    , mMasterContactFace3D( 0, Element::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mMasterContactFace3D3( 0, Element::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mMasterContactFace3D6( 0, Element::GeometryType::Pointer( new Triangle3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) )
    , mMasterContactFace3D8( 0, Element::GeometryType::Pointer( new Quadrilateral3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mMasterContactFace3D9( 0, Element::GeometryType::Pointer( new Quadrilateral3D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) )
    , mSlaveContactFace3D( 0, Element::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mSlaveContactFace3D3( 0, Element::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mSlaveContactFace3D6( 0, Element::GeometryType::Pointer( new Triangle3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) )
    , mSlaveContactFace3D8( 0, Element::GeometryType::Pointer( new Quadrilateral3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mSlaveContactFace3D9( 0, Element::GeometryType::Pointer( new Quadrilateral3D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) )
    , mMasterContactFace3DNewmark( 0, Element::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mMasterContactFace3D3Newmark( 0, Element::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mMasterContactFace3D6Newmark( 0, Element::GeometryType::Pointer( new Triangle3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) )
    , mMasterContactFace3D8Newmark( 0, Element::GeometryType::Pointer( new Quadrilateral3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mMasterContactFace3D9Newmark( 0, Element::GeometryType::Pointer( new Quadrilateral3D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) )
    , mSlaveContactFace3DNewmark( 0, Element::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mSlaveContactFace3D3Newmark( 0, Element::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mSlaveContactFace3D6Newmark( 0, Element::GeometryType::Pointer( new Triangle3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) )
    , mSlaveContactFace3D8Newmark( 0, Element::GeometryType::Pointer( new Quadrilateral3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mSlaveContactFace3D9Newmark( 0, Element::GeometryType::Pointer( new Quadrilateral3D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) )
    , mFaceVel3D3N( 0, Element::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    //, mUPCTestElement3D20N(0,Element::GeometryType::Pointer( new Hexahedra3D20<Node<3> >(Element::GeometryType::PointsArrayType(20, Node<3>()))))
    //, mPointForce3D(0, Element::GeometryType::Pointer(new Geometry <Node<3> >(Element::GeometryType::PointsArrayType(1, Node<3>()))))
    //, mPointForce2D(0, Element::GeometryType::Pointer(new Geometry <Node<3> >(Element::GeometryType::PointsArrayType(1, Node<3>()))))

    , mPointForce3D( 0, Element::GeometryType::Pointer( new Point3D <Node<3> >( Element::GeometryType::PointsArrayType( 1, Node<3>() ) ) ) )
    , mPointForce2D( 0, Element::GeometryType::Pointer( new Point2D <Node<3> >( Element::GeometryType::PointsArrayType( 1, Node<3>() ) ) ) )
    , mPointMoment3D( 0, Element::GeometryType::Pointer( new Point3D <Node<3> >( Element::GeometryType::PointsArrayType( 1, Node<3>() ) ) ) )

    , mElasticPointConstraint( 0, Element::GeometryType::Pointer( new Point3D <Node<3> >( Element::GeometryType::PointsArrayType( 1, Node<3>() ) ) ) )
    , mElasticLineConstraint2N( 0, Element::GeometryType::Pointer( new Line3D2<Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) )
    , mElasticLineConstraint3N( 0, Element::GeometryType::Pointer( new Line3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mElasticFaceConstraint3N( 0, Element::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mElasticFaceConstraint6N( 0, Element::GeometryType::Pointer( new Triangle3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) )
    , mElasticFaceConstraint4N( 0, Element::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mElasticFaceConstraint8N( 0, Element::GeometryType::Pointer( new Quadrilateral3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mElasticFaceConstraint9N( 0, Element::GeometryType::Pointer( new Quadrilateral3D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) )
    , mElasticFaceSprings3N( 0, Element::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mElasticFaceSprings6N( 0, Element::GeometryType::Pointer( new Triangle3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) )
    , mElasticFaceSprings4N( 0, Element::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mElasticFaceSprings8N( 0, Element::GeometryType::Pointer( new Quadrilateral3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mElasticFaceSprings9N( 0, Element::GeometryType::Pointer( new Quadrilateral3D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) )
    , mNitscheIsotropicConstraint2D2N( 0, Element::GeometryType::Pointer( new Line2D2 <Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) )
    , mNitscheIsotropicConstraint2D3N( 0, Element::GeometryType::Pointer( new Line2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mNitscheIsotropicConstraint3D3N( 0, Element::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mNitscheIsotropicConstraint3D6N( 0, Element::GeometryType::Pointer( new Triangle3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) )
    , mNitscheIsotropicConstraint3D4N( 0, Element::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mNitscheIsotropicConstraint3D8N( 0, Element::GeometryType::Pointer( new Quadrilateral3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mNitscheIsotropicConstraint3D9N( 0, Element::GeometryType::Pointer( new Quadrilateral3D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) )

    , mNodeTyingLagrange( 0, Element::GeometryType::Pointer( new Geometry <Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) )
    , mNodeTyingLagrangeZ( 0, Element::GeometryType::Pointer( new Geometry <Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) )
    , mPointPointJointCondition( 0, Element::GeometryType::Pointer( new Line3D2 <Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) )
    , mPointPointLagrangeCondition( 0, Element::GeometryType::Pointer( new Line3D2 <Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) )

    , mSlaveContactPoint2D( 0, Element::GeometryType::Pointer( new Point2D <Node<3> >( Element::GeometryType::PointsArrayType( 1, Node<3>() ) ) ) )
    , mMasterContactFace2D( 0, Element::GeometryType::Pointer( new Line2D2 <Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) )
    , mIsotropic3D()
    , mDummyConstitutiveLaw()
    , mDruckerPrager()
    , mCamClay3D()
#endif
{}

void KratosStructuralApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing KratosStructuralApplication... " << std::endl;

    KRATOS_REGISTER_VARIABLE( DAMAGE_E0 )
    KRATOS_REGISTER_VARIABLE( DAMAGE_EF )

    KRATOS_REGISTER_VARIABLE( MATRIX_A )
    KRATOS_REGISTER_VARIABLE( MATRIX_B )
    KRATOS_REGISTER_VARIABLE( MATRIX_D )
    KRATOS_REGISTER_VARIABLE( COMPOSITE_DIRECTION )
    KRATOS_REGISTER_VARIABLE( JOINT_STIFFNESS )
    KRATOS_REGISTER_VARIABLE( LINING_JOINT_STIFFNESS )
    KRATOS_REGISTER_VARIABLE( ORTHOTROPIC_YOUNG_MODULUS )
    KRATOS_REGISTER_VARIABLE( ORTHOTROPIC_SHEAR_MODULUS )
    KRATOS_REGISTER_VARIABLE( ORTHOTROPIC_POISSON_RATIO )
    KRATOS_REGISTER_VARIABLE( MATERIAL_DIRECTION )
    KRATOS_REGISTER_VARIABLE( GEOMETRIC_STIFFNESS )
    KRATOS_REGISTER_VARIABLE( INSITU_STRESS_SCALE )
    KRATOS_REGISTER_VARIABLE( LAGRANGE_SCALE )
    KRATOS_REGISTER_VARIABLE( REFERENCE_PRESSURE )
    KRATOS_REGISTER_VARIABLE( REFERENCE_WATER_PRESSURE )
    KRATOS_REGISTER_VARIABLE( WATER_PRESSURE_SCALE )
    KRATOS_REGISTER_VARIABLE( OVERCONSOLIDATION_RATIO )
    KRATOS_REGISTER_VARIABLE( EXCESS_PORE_WATER_PRESSURE )
    KRATOS_REGISTER_VARIABLE( PRESSURE_P )
    KRATOS_REGISTER_VARIABLE( PRESSURE_Q )
    KRATOS_REGISTER_VARIABLE( COORDINATES )
    // KRATOS_REGISTER_VARIABLE( STRESSES )
    KRATOS_REGISTER_VARIABLE( STRAIN )
    KRATOS_REGISTER_VARIABLE( FLUID_FLOWS )
    KRATOS_REGISTER_VARIABLE( CONTACT_PENETRATION )
    KRATOS_REGISTER_VARIABLE( DAMPING_RATIO )
    //KRATOS_REGISTER_VARIABLE( KINETIC_ENERGY )
    KRATOS_REGISTER_VARIABLE( POTENCIAL_ENERGY )
    KRATOS_REGISTER_VARIABLE( DEFORMATION_ENERGY )
    KRATOS_REGISTER_VARIABLE( VON_MISES_STRESS )
    KRATOS_REGISTER_VARIABLE( RHS_PRESSURE )


    //  KRATOS_REGISTER_VARIABLE(WRINKLING_APPROACH )
//  KRATOS_REGISTER_VARIABLE(GREEN_LAGRANGE_STRAIN_TENSOR )
//  KRATOS_REGISTER_VARIABLE(PK2_STRESS_TENSOR )
//  KRATOS_REGISTER_VARIABLE(AUXILIARY_MATRIX_1 )
//  KRATOS_REGISTER_VARIABLE(YOUNG_MODULUS )
//  KRATOS_REGISTER_VARIABLE(POISSON_RATIO )
//  KRATOS_REGISTER_VARIABLE(MU )
    KRATOS_REGISTER_VARIABLE( ALPHA )
    KRATOS_REGISTER_VARIABLE( RETRACTION_TIME )
//  KRATOS_REGISTER_VARIABLE(THICKNESS )
//  KRATOS_REGISTER_VARIABLE(NEGATIVE_FACE_PRESSURE )
//  KRATOS_REGISTER_VARIABLE(POSITIVE_FACE_PRESSURE )

    //KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VAUX);
    KRATOS_REGISTER_VARIABLE( CONSTITUTIVE_LAW_NO_INITIALIZE )
//     KRATOS_REGISTER_VARIABLE(DP_EPSILON)
//     KRATOS_REGISTER_VARIABLE(INSITU_STRESS)
//     KRATOS_REGISTER_VARIABLE(DP_ALPHA1)
//     KRATOS_REGISTER_VARIABLE(DP_K)
//     KRATOS_REGISTER_VARIABLE(CALCULATE_INSITU_STRESS)
    //CONTACT_LINK_MASTER is defined in condition.h
    KRATOS_REGISTER_VARIABLE( CONTACT_LINK_MASTER )
    //CONTACT_LINK_SLAVE is defined in condition.h
    // KRATOS_REGISTER_VARIABLE( NEAR_NODE )
    KRATOS_REGISTER_VARIABLE( CONTACT_LINK_SLAVE )
    // KRATOS_REGISTER_VARIABLE( MASTER_CONTACT_LOCAL_POINT )
    // KRATOS_REGISTER_VARIABLE( MASTER_CONTACT_CURRENT_LOCAL_POINT )
    // KRATOS_REGISTER_VARIABLE( MASTER_CONTACT_LAST_CURRENT_LOCAL_POINT )
    // KRATOS_REGISTER_VARIABLE( SLAVE_CONTACT_LOCAL_POINT )
    // KRATOS_REGISTER_VARIABLE( MASTER_CONTACT_GLOBAL_POINT )
    // KRATOS_REGISTER_VARIABLE( MASTER_CONTACT_CURRENT_GLOBAL_POINT )
    // KRATOS_REGISTER_VARIABLE( SLAVE_CONTACT_GLOBAL_POINT )


    KRATOS_REGISTER_VARIABLE( BASE )
    KRATOS_REGISTER_VARIABLE( HEIGHT )
    KRATOS_REGISTER_VARIABLE( CROSS_AREA )
    KRATOS_REGISTER_VARIABLE( AREA )
    KRATOS_REGISTER_VARIABLE( AREA_X )
    KRATOS_REGISTER_VARIABLE( AREA_Y )
    KRATOS_REGISTER_VARIABLE( AREA_Z )
    KRATOS_REGISTER_VARIABLE( INERTIA )
    KRATOS_REGISTER_VARIABLE( LOCAL_INERTIA )
    KRATOS_REGISTER_VARIABLE( INERTIA_X )
    KRATOS_REGISTER_VARIABLE( INERTIA_Y )
    KRATOS_REGISTER_VARIABLE( INERTIA_Z )
    KRATOS_REGISTER_VARIABLE( FC )
    KRATOS_REGISTER_VARIABLE( FT )
    KRATOS_REGISTER_VARIABLE( CONCRETE_YOUNG_MODULUS_C )
    KRATOS_REGISTER_VARIABLE( CONCRETE_YOUNG_MODULUS_T )
    KRATOS_REGISTER_VARIABLE( FRACTURE_ENERGY )
    KRATOS_REGISTER_VARIABLE( CRUSHING_ENERGY )
    KRATOS_REGISTER_VARIABLE( PLASTIC_ENERGY )
    KRATOS_REGISTER_VARIABLE( ELASTIC_ENERGY )
//         KRATOS_REGISTER_VARIABLE( YIELD_STRESS )
    KRATOS_REGISTER_VARIABLE( PLASTIC_MODULUS )
    KRATOS_REGISTER_VARIABLE( PLASTICITY_INDICATOR )
    KRATOS_REGISTER_VARIABLE( LAMNDA ) // Load factor
    KRATOS_REGISTER_VARIABLE( DAMAGE )
    KRATOS_REGISTER_VARIABLE( ORTHOTROPIC_ANGLE )
    KRATOS_REGISTER_VARIABLE( VOLUMEN_FRACTION )
    KRATOS_REGISTER_VARIABLE( MAX_INTERNAL_FRICTION_ANGLE )
    KRATOS_REGISTER_VARIABLE( DILATANCY_ANGLE )
    KRATOS_REGISTER_VARIABLE( MAX_DILATANCY_ANGLE )
    KRATOS_REGISTER_VARIABLE( COHESION )
    KRATOS_REGISTER_VARIABLE( ISOTROPIC_ELASTIC_LIMIT )
    KRATOS_REGISTER_VARIABLE( ORTHOTROPIC_ELASTIC_LIMIT )
    KRATOS_REGISTER_VARIABLE( VECTOR_DAMAGE )
    KRATOS_REGISTER_VARIABLE( ORTHOTROPIC_YOUNG_MODULUS_2D ) // [E1 E2 G12]
    KRATOS_REGISTER_VARIABLE( ORTHOTROPIC_POISSON_RATIO_2D ) // [v12 v21]
    KRATOS_REGISTER_VARIABLE( GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR )
    KRATOS_REGISTER_VARIABLE( ELASTIC_STRAIN_VECTOR )
    KRATOS_REGISTER_VARIABLE( PLASTIC_STRAIN_VECTOR )
    KRATOS_REGISTER_VARIABLE( CURRENT_STRAIN_VECTOR )
    KRATOS_REGISTER_VARIABLE( INTEGRATION_POINT_STRAIN_VECTOR )
    KRATOS_REGISTER_VARIABLE( EQUIVALENT_VOLUMETRIC_ELASTIC_STRAIN )
    KRATOS_REGISTER_VARIABLE( EQUIVALENT_VOLUMETRIC_PLASTIC_STRAIN )
    KRATOS_REGISTER_VARIABLE( EQUIVALENT_DEVIATORIC_ELASTIC_STRAIN )
    KRATOS_REGISTER_VARIABLE( EQUIVALENT_DEVIATORIC_PLASTIC_STRAIN )
    KRATOS_REGISTER_VARIABLE( EQUIVALENT_VOLUMETRIC_STRAIN )
    KRATOS_REGISTER_VARIABLE( EQUIVALENT_DEVIATORIC_STRAIN )
    KRATOS_REGISTER_VARIABLE( EQUIVALENT_STRAIN )
    KRATOS_REGISTER_VARIABLE( ALMANSI_PLASTIC_STRAIN )
    KRATOS_REGISTER_VARIABLE( ALMANSI_ELASTIC_STRAIN )
    KRATOS_REGISTER_VARIABLE( PRESTRESS )
    KRATOS_REGISTER_VARIABLE( PRESTRESS_ZZ )
    KRATOS_REGISTER_VARIABLE( PRESTRESS_FACTOR )
    KRATOS_REGISTER_VARIABLE( INITIAL_STRESS )
    // KRATOS_REGISTER_VARIABLE( MAX_FRECUENCY )

    KRATOS_REGISTER_VARIABLE( DISIPATION )
    KRATOS_REGISTER_VARIABLE( ISOTROPIC_HARDENING_MODULUS )
    KRATOS_REGISTER_VARIABLE( KINEMATIC_HARDENING_MODULUS )
    KRATOS_REGISTER_VARIABLE( ISOTROPIC_HARDENING_TYPE )
    KRATOS_REGISTER_VARIABLE( KINEMATIC_HARDENING_TYPE )
    KRATOS_REGISTER_VARIABLE( DRUCKER_PRAGER_MATCHING_TYPE )
    KRATOS_REGISTER_VARIABLE( HARDENING_POINTS_ON_CURVE )
    KRATOS_REGISTER_VARIABLE( NODAL_STRESS )
    KRATOS_REGISTER_VARIABLE( NODAL_STRESS_VECTOR )
    KRATOS_REGISTER_VARIABLE( NODAL_STRAIN )
    KRATOS_REGISTER_VARIABLE( NODAL_STRAIN_VECTOR )
    KRATOS_REGISTER_VARIABLE( NODAL_VALUES )
    KRATOS_REGISTER_VARIABLE( NODAL_DAMAGE )
    KRATOS_REGISTER_VARIABLE( IS_TARGET )
    KRATOS_REGISTER_VARIABLE( IS_CONTACTOR )
    KRATOS_REGISTER_VARIABLE( COMPUTE_TANGENT_MATRIX )
    KRATOS_REGISTER_VARIABLE( IS_DISCRETE )
    KRATOS_REGISTER_VARIABLE( CONSTRAINT_MATRIX )
    KRATOS_REGISTER_VARIABLE( CONSTRAINT_VECTOR )
    KRATOS_REGISTER_VARIABLE( YIELD_STATE )
    KRATOS_REGISTER_VARIABLE( YIELD_SURFACE )
    KRATOS_REGISTER_VARIABLE( TENSILE_STRENGTH )
    KRATOS_REGISTER_VARIABLE( SHEAR_STRENGTH )
    KRATOS_REGISTER_VARIABLE( VISCOUS_DAMPING )
    KRATOS_REGISTER_VARIABLE( CONSTRAINT_MATRIX )
    KRATOS_REGISTER_VARIABLE( CONSTRAINT_VECTOR )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( JOINT_FORCE_REACTION )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( JOINT_MOMENT_REACTION )
    //KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( INTERNAL_FORCE ) //already put on variables.cpp (warning was appearing on Windows)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( ELASTIC_BEDDING_STIFFNESS )
    KRATOS_REGISTER_VARIABLE( IS_BBAR )
    KRATOS_REGISTER_VARIABLE( NEIGHBOUR_EXPANSION_LEVEL )
    KRATOS_REGISTER_VARIABLE( STRESS_RECOVERY_TYPE )
    KRATOS_REGISTER_VARIABLE( RAYLEIGH_DAMPING_ALPHA )
    KRATOS_REGISTER_VARIABLE( RAYLEIGH_DAMPING_BETA )
    KRATOS_REGISTER_VARIABLE( STABILISATION_FACTOR )
    KRATOS_REGISTER_VARIABLE( SHEAR_MODULUS )
    KRATOS_REGISTER_VARIABLE( SHEAR_MODULUS_EVOLUTION )
    KRATOS_REGISTER_VARIABLE( RECOVERY_STRESSES )
    KRATOS_REGISTER_VARIABLE( THREED_STRESSES )
    KRATOS_REGISTER_VARIABLE( THREED_PRESTRESS )
    KRATOS_REGISTER_VARIABLE( STRESSES_OLD )
    KRATOS_REGISTER_VARIABLE( THREED_STRAIN )
    KRATOS_REGISTER_VARIABLE( PRE_STRAIN_VECTOR )
    KRATOS_REGISTER_VARIABLE( POST_STRAIN_VECTOR )
    KRATOS_REGISTER_VARIABLE( STRESS_TENSOR )
    KRATOS_REGISTER_VARIABLE( STRAIN_TENSOR )
    KRATOS_REGISTER_VARIABLE( ELASTIC_STRAIN_TENSOR )
    KRATOS_REGISTER_VARIABLE( PLASTIC_STRAIN_TENSOR )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( PRESCRIBED_DELTA_DISPLACEMENT )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( LINE_LOAD )
    KRATOS_REGISTER_VARIABLE( IS_CONTACT_NODE )
    KRATOS_REGISTER_VARIABLE( LAGRANGE_MULTIPLIER )
    KRATOS_REGISTER_VARIABLE( HAS_STRAIN_AT_NODE )
    KRATOS_REGISTER_VARIABLE( HAS_STRESSES_AT_NODE )
    KRATOS_REGISTER_VARIABLE( HAS_NODAL_ERROR )
    KRATOS_REGISTER_VARIABLE( FORCE_EQUAL_ORDER_INTERPOLATION )
    KRATOS_REGISTER_VARIABLE( NODAL_ERROR_1 )
    KRATOS_REGISTER_VARIABLE( DUMMY_DOF )

    KRATOS_REGISTER_VARIABLE( PRECONSOLIDATION_PRESSURE )
    KRATOS_REGISTER_VARIABLE( PRECONSOLIDATION_PRESSURE_DEF )
    KRATOS_REGISTER_VARIABLE( PRECONSOLIDATION_PRESSURE_MIN )
    KRATOS_REGISTER_VARIABLE( CSL_SLOPE )
    KRATOS_REGISTER_VARIABLE( VIRGIN_COMPRESSION_INDEX )
    KRATOS_REGISTER_VARIABLE( RECOMPRESSION_INDEX )
    KRATOS_REGISTER_VARIABLE( SWELL_INDEX )
    KRATOS_REGISTER_VARIABLE( VOID_RATIO )
    KRATOS_REGISTER_VARIABLE( SPACING_RATIO )
    KRATOS_REGISTER_VARIABLE( ASSOCIATIVITY )
    KRATOS_REGISTER_VARIABLE( SHAPE_PARAMETER )

    KRATOS_REGISTER_VARIABLE( NEIGHBOUR_WEIGHTS )
    KRATOS_REGISTER_VARIABLE( MATERIAL_DENSITY )
    KRATOS_REGISTER_VARIABLE( MATERIAL_DENSITY_NEW )
    KRATOS_REGISTER_VARIABLE( MATERIAL_DENSITY_FILTERED )
    KRATOS_REGISTER_VARIABLE( ELEMENT_DC )
    KRATOS_REGISTER_VARIABLE( ELEMENT_DC_FILTERED )
    KRATOS_REGISTER_VARIABLE( ELEMENT_DV )
    KRATOS_REGISTER_VARIABLE( ELEMENT_DV_FILTERED )
    KRATOS_REGISTER_VARIABLE( YOUNG_MODULUS_0 )
    KRATOS_REGISTER_VARIABLE( YOUNG_MODULUS_MIN )
    KRATOS_REGISTER_VARIABLE( PENALIZATION_FACTOR )
    KRATOS_REGISTER_VARIABLE( GEOMETRICAL_DOMAIN_SIZE )
    KRATOS_REGISTER_VARIABLE( JACOBIAN_0 )
    // KRATOS_REGISTER_VARIABLE( INTEGRATION_WEIGHT )
    // KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( INTEGRATION_POINT_GLOBAL )
    // KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( INTEGRATION_POINT_LOCAL )
    KRATOS_REGISTER_VARIABLE( ASSOCIATED_NODE )
    KRATOS_REGISTER_VARIABLE( ASSOCIATED_ELEMENT )
    KRATOS_REGISTER_VARIABLE( ASSOCIATED_CONDITION )
    KRATOS_REGISTER_VARIABLE( PRINCIPAL_STRESS )
    KRATOS_REGISTER_VARIABLE( PRINCIPAL_STRAIN )
    KRATOS_REGISTER_VARIABLE( ALGORITHMIC_TANGENT )
    KRATOS_REGISTER_VARIABLE( ELASTIC_TANGENT )
    KRATOS_REGISTER_VARIABLE( IS_MARKED_FOR_REACTION )
    KRATOS_REGISTER_VARIABLE( EXTRAPOLATION_FACTOR )
    KRATOS_REGISTER_VARIABLE( FIX_POROSITY )
    KRATOS_REGISTER_VARIABLE( QUAD_POINT_STATUS )
    KRATOS_REGISTER_VARIABLE( LOCAL_FRAME )
    KRATOS_REGISTER_VARIABLE( STRESS_LIKE_INTERNAL_VARIABLES )
    KRATOS_REGISTER_VARIABLE( STRAIN_LIKE_INTERNAL_VARIABLES )
    KRATOS_REGISTER_VARIABLE ( PRIMARY_HYDRATION_TIME )
    KRATOS_REGISTER_VARIABLE ( PRIMARY_HYDRATION_TIME_GRADIENT )
    KRATOS_REGISTER_VARIABLE ( STIFFNESS_RATIO )
    KRATOS_REGISTER_VARIABLE ( L2_ERROR )
    KRATOS_REGISTER_VARIABLE ( H1_ERROR )
    KRATOS_REGISTER_VARIABLE ( ROTATIONAL_STIFFNESS )


    //KRATOS_REGISTER_VARIABLE(TO_ERASE )
//   KRATOS_REGISTER_VARIABLE(CONTACT_RAMP )
//   KRATOS_REGISTER_VARIABLE(PENALTY )
// //   KRATOS_REGISTER_VARIABLE(INITIAL_PENALTY )
//   KRATOS_REGISTER_VARIABLE(MAXIMUM_PENALTY )
//   KRATOS_REGISTER_VARIABLE(RAMP_CRITERION )
//   KRATOS_REGISTER_VARIABLE(RAMP_FACTOR )
//   KRATOS_REGISTER_VARIABLE(PENALTY_T )
//   KRATOS_REGISTER_VARIABLE(INITIAL_PENALTY_T )
//   KRATOS_REGISTER_VARIABLE(MAXIMUM_PENALTY_T )
//   KRATOS_REGISTER_VARIABLE(RAMP_CRITERION_T )
//   KRATOS_REGISTER_VARIABLE(RAMP_FACTOR_T )
//   KRATOS_REGISTER_VARIABLE(FRICTION_COEFFICIENT )
//   KRATOS_REGISTER_VARIABLE(LAMBDAS )
//   KRATOS_REGISTER_VARIABLE(LAMBDAS_T )
//   KRATOS_REGISTER_VARIABLE(GAPS )
//   KRATOS_REGISTER_VARIABLE(DELTA_LAMBDAS )
//   KRATOS_REGISTER_VARIABLE(DELTA_LAMBDAS_T )
//   KRATOS_REGISTER_VARIABLE(MAX_UZAWA_ITERATIONS)
//   KRATOS_REGISTER_VARIABLE(CONTACT_SLAVE_INTEGRATION_POINT_INDEX )
//   KRATOS_REGISTER_VARIABLE(CONTACT_LINK_M )
//   KRATOS_REGISTER_VARIABLE(CONTACT_DOUBLE_CHECK )
//   KRATOS_REGISTER_VARIABLE(IS_CONTACT_MASTER )
//   KRATOS_REGISTER_VARIABLE(IS_CONTACT_SLAVE )
//   KRATOS_REGISTER_VARIABLE(K_CONTACT )
//   KRATOS_REGISTER_VARIABLE(K_CONTACT_T )
//   KRATOS_REGISTER_VARIABLE(STICK)
//   KRATOS_REGISTER_VARIABLE(FIRST_TIME_STEP)
//   KRATOS_REGISTER_VARIABLE(QUASI_STATIC_ANALYSIS )
//   KRATOS_REGISTER_VARIABLE( NORMAL_STRESS )
//   KRATOS_REGISTER_VARIABLE( TANGENTIAL_STRESS )
//   KRATOS_REGISTER_VARIABLE( NORMAL_CONTACT_STRESS )
//   KRATOS_REGISTER_VARIABLE( TANGENTIAL_CONTACT_STRESS )
//   KRATOS_REGISTER_VARIABLE( CONTACT_STICK )


//   KRATOS_REGISTER_VARIABLE(WATER_PRESSURE)
//   KRATOS_REGISTER_VARIABLE(WATER_PRESSURE_DT)
//   KRATOS_REGISTER_VARIABLE(WATER_PRESSURE_ACCELERATION)
//   KRATOS_REGISTER_VARIABLE(WATER_PRESSURE_NULL)
//   KRATOS_REGISTER_VARIABLE(WATER_PRESSURE_NULL_DT)
//   KRATOS_REGISTER_VARIABLE(WATER_PRESSURE_NULL_ACCELERATION)
//   KRATOS_REGISTER_VARIABLE(WATER_PRESSURE_EINS)
//   KRATOS_REGISTER_VARIABLE(WATER_PRESSURE_EINS_DT)
//   KRATOS_REGISTER_VARIABLE(WATER_PRESSURE_EINS_ACCELERATION)
//
//   KRATOS_REGISTER_VARIABLE(AIR_PRESSURE)
//   KRATOS_REGISTER_VARIABLE(AIR_PRESSURE_DT)
//   KRATOS_REGISTER_VARIABLE(AIR_PRESSURE_ACCELERATION)
//   KRATOS_REGISTER_VARIABLE(AIR_PRESSURE_NULL)
//   KRATOS_REGISTER_VARIABLE(AIR_PRESSURE_NULL_DT)
//   KRATOS_REGISTER_VARIABLE(AIR_PRESSURE_NULL_ACCELERATION)
//   KRATOS_REGISTER_VARIABLE(AIR_PRESSURE_EINS)
//   KRATOS_REGISTER_VARIABLE(AIR_PRESSURE_EINS_DT)
//   KRATOS_REGISTER_VARIABLE(AIR_PRESSURE_EINS_ACCELERATION)
//
//   KRATOS_REGISTER_VARIABLE(DISPLACEMENT_OLD)
//   KRATOS_REGISTER_VARIABLE(DISPLACEMENT_DT)
// //   KRATOS_REGISTER_VARIABLE(ACCELERATION)
//   KRATOS_REGISTER_VARIABLE(DISPLACEMENT_NULL)
//   KRATOS_REGISTER_VARIABLE(DISPLACEMENT_NULL_DT)
//   KRATOS_REGISTER_VARIABLE(ACCELERATION_NULL)
//   KRATOS_REGISTER_VARIABLE(DISPLACEMENT_EINS)
//   KRATOS_REGISTER_VARIABLE(DISPLACEMENT_EINS_DT)
//   KRATOS_REGISTER_VARIABLE(ACCELERATION_EINS)
//   KRATOS_REGISTER_VARIABLE(ELASTIC_LEFT_CAUCHY_GREEN_OLD)
//
//   KRATOS_REGISTER_VARIABLE(ACTIVATION_LEVEL)


    #ifdef SD_APP_FORWARD_COMPATIBILITY

    // register legacy variables
    KRATOS_REGISTER_VARIABLE( LAMBDA_OLD )
    KRATOS_REGISTER_VARIABLE( LAMBDA_NULL )
    KRATOS_REGISTER_VARIABLE( LAMBDA_EINS )
    KRATOS_REGISTER_VARIABLE( LAMBDA_DT )
    KRATOS_REGISTER_VARIABLE( LAMBDA_NULL_DT )
    KRATOS_REGISTER_VARIABLE( LAMBDA_EINS_DT )
    KRATOS_REGISTER_VARIABLE( LAMBDA_DT2 )
    KRATOS_REGISTER_VARIABLE( LAMBDA_NULL_DT2 )
    KRATOS_REGISTER_VARIABLE( LAMBDA_EINS_DT2 )
    KRATOS_REGISTER_VARIABLE( LAMBDAS_T )
    KRATOS_REGISTER_VARIABLE( DELTA_LAMBDAS_T )
    KRATOS_REGISTER_VARIABLE( CONTACT_LINK_M )
    KRATOS_REGISTER_VARIABLE( AUXILIARY_MATRIX_1 )
    KRATOS_REGISTER_VARIABLE( ELASTIC_LEFT_CAUCHY_GREEN_OLD )

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( DISPLACEMENT_OLD )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( DISPLACEMENT_DT )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( DISPLACEMENT_NULL )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( DISPLACEMENT_NULL_DT )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( ACCELERATION_NULL )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( DISPLACEMENT_EINS )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( DISPLACEMENT_EINS_DT )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( ACCELERATION_EINS )

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( ROTATION_OLD )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( ROTATION_DT )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( ROTATION_NULL )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( ROTATION_NULL_DT )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( ANGULAR_ACCELERATION_NULL )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( ROTATION_EINS )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( ROTATION_EINS_DT )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( ANGULAR_ACCELERATION_EINS )

    KRATOS_REGISTER_VARIABLE( PENALTY_T )
    KRATOS_REGISTER_VARIABLE( INITIAL_PENALTY_T )
    KRATOS_REGISTER_VARIABLE( LAMBDAS )
    KRATOS_REGISTER_VARIABLE( GAPS )
    KRATOS_REGISTER_VARIABLE( DELTA_LAMBDAS )
    KRATOS_REGISTER_VARIABLE( STICK )
    KRATOS_REGISTER_VARIABLE( CONTACT_STICK )

    KRATOS_REGISTER_VARIABLE( USE_DISTRIBUTED_PROPERTIES )

    KRATOS_REGISTER_VARIABLE( CONTACT_RAMP )
    KRATOS_REGISTER_VARIABLE( CONTACT_SLAVE_INTEGRATION_POINT_INDEX )
    KRATOS_REGISTER_VARIABLE( CONTACT_DOUBLE_CHECK )
//        KRATOS_REGISTER_VARIABLE( FIRST_TIME_STEP )
//        KRATOS_REGISTER_VARIABLE( QUASI_STATIC_ANALYSIS )

    KRATOS_REGISTER_VARIABLE( INTEGRATION_ORDER )
    KRATOS_REGISTER_VARIABLE( INTEGRATION_WEIGHT )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( INTEGRATION_POINT_GLOBAL )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( INTEGRATION_POINT_GLOBAL_IN_REFERENCE_CONFIGURATION )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( INTEGRATION_POINT_GLOBAL_IN_CURRENT_CONFIGURATION )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( INTEGRATION_POINT_LOCAL )

//        KRATOS_REGISTER_VARIABLE( OSS_SWITCH )
    KRATOS_REGISTER_VARIABLE( WRINKLING_APPROACH )
    KRATOS_REGISTER_VARIABLE( CALCULATE_INSITU_STRESS )
//        KRATOS_REGISTER_VARIABLE( PERIODIC_PAIR_INDEX )
//        KRATOS_REGISTER_VARIABLE( STATIONARY )

//        KRATOS_REGISTER_VARIABLE( PARTITION_MASK )

//        KRATOS_REGISTER_VARIABLE( FACE_HEAT_FLUX )

//        KRATOS_REGISTER_VARIABLE( NODAL_VOLUME )

    KRATOS_REGISTER_VARIABLE( PRESSURE_DT )
    KRATOS_REGISTER_VARIABLE( WATER_PRESSURE_DT )
    KRATOS_REGISTER_VARIABLE( WATER_PRESSURE_ACCELERATION )
    KRATOS_REGISTER_VARIABLE( WATER_PRESSURE_NULL )
    KRATOS_REGISTER_VARIABLE( WATER_PRESSURE_NULL_DT )
    KRATOS_REGISTER_VARIABLE( WATER_PRESSURE_NULL_ACCELERATION )
    KRATOS_REGISTER_VARIABLE( WATER_PRESSURE_EINS )
    KRATOS_REGISTER_VARIABLE( WATER_PRESSURE_EINS_DT )
    KRATOS_REGISTER_VARIABLE( WATER_PRESSURE_EINS_ACCELERATION )

    KRATOS_REGISTER_VARIABLE( AIR_PRESSURE_DT )
    KRATOS_REGISTER_VARIABLE( AIR_PRESSURE_ACCELERATION )
    KRATOS_REGISTER_VARIABLE( AIR_PRESSURE_NULL )
    KRATOS_REGISTER_VARIABLE( AIR_PRESSURE_NULL_DT )
    KRATOS_REGISTER_VARIABLE( AIR_PRESSURE_NULL_ACCELERATION )
    KRATOS_REGISTER_VARIABLE( AIR_PRESSURE_EINS )
    KRATOS_REGISTER_VARIABLE( AIR_PRESSURE_EINS_DT )
    KRATOS_REGISTER_VARIABLE( AIR_PRESSURE_EINS_ACCELERATION )
    KRATOS_REGISTER_VARIABLE( SUCTION )
//        KRATOS_REGISTER_VARIABLE( FLAG_VARIABLE )

//        KRATOS_REGISTER_VARIABLE( DP_EPSILON )
//        KRATOS_REGISTER_VARIABLE( DP_ALPHA1 )
//        KRATOS_REGISTER_VARIABLE( DP_K )

//        KRATOS_REGISTER_VARIABLE( EQ_STRAIN_RATE )
//        KRATOS_REGISTER_VARIABLE( RHS_WATER )
//        KRATOS_REGISTER_VARIABLE( RHS_AIR )

    KRATOS_REGISTER_VARIABLE( PERMEABILITY_28_DAYS )
    KRATOS_REGISTER_VARIABLE( PERMEABILITY_1_DAY )
    KRATOS_REGISTER_VARIABLE( PERMEABILITY_TRANSITION )

//        KRATOS_REGISTER_VARIABLE( TEMPERATURE_OLD_IT )
//        KRATOS_REGISTER_VARIABLE( EFFECTIVE_VISCOSITY )
//        KRATOS_REGISTER_VARIABLE( KINEMATIC_VISCOSITY)
//        KRATOS_REGISTER_VARIABLE( DYNAMIC_VISCOSITY)
//        KRATOS_REGISTER_VARIABLE( WEIGHT_FATHER_NODES )

    KRATOS_REGISTER_VARIABLE( PARENT_ELEMENT_ID )
    KRATOS_REGISTER_VARIABLE( INTEGRATION_POINT_INDEX )

    KRATOS_REGISTER_VARIABLE( IS_CONTACT_MASTER )
    KRATOS_REGISTER_VARIABLE( IS_CONTACT_SLAVE )
//        KRATOS_REGISTER_VARIABLE( IS_BOUNDARY )
//        KRATOS_REGISTER_VARIABLE( IS_VISITED )

    KRATOS_REGISTER_VARIABLE( LAGRANGE_MULTIPLIER_CONSTRAINT )
    KRATOS_REGISTER_VARIABLE( LAGRANGE_MULTIPLIER_CONSTRAINT_REACTION )

    KRATOS_REGISTER_VARIABLE( IS_SHAPE_FUNCTION_REQUIRED )
    KRATOS_REGISTER_VARIABLE( RESET_CONFIGURATION )
    #endif

    /// Register elements and conditions

    #ifdef SD_APP_FORWARD_COMPATIBILITY
    KRATOS_REGISTER_ELEMENT( "CrisfieldTrussElement3D2N", mCrisfieldTrussElement3D2N )
    KRATOS_REGISTER_ELEMENT( "CrisfieldTrussElement3D3N", mCrisfieldTrussElement3D3N )
    KRATOS_REGISTER_ELEMENT( "TrussElement3D2N", mTrussElement3D2N )
    KRATOS_REGISTER_ELEMENT( "TrussElement3D3N", mTrussElement3D3N )
    KRATOS_REGISTER_ELEMENT( "BeamElement3D2N", mBeamElement3D2N )
    KRATOS_REGISTER_ELEMENT( "BeamElement3D3N", mBeamElement3D3N )
    KRATOS_REGISTER_ELEMENT( "TimoshenkoBeamElement3D2N", mTimoshenkoBeamElement3D2N )
    KRATOS_REGISTER_ELEMENT( "TimoshenkoBeamElement3D3N", mTimoshenkoBeamElement3D3N )
    KRATOS_REGISTER_ELEMENT( "TimoshenkoLinearBeamElement2D2N", mTimoshenkoLinearBeamElement2D2N )
    KRATOS_REGISTER_ELEMENT( "TimoshenkoLinearBeamElement3D2N", mTimoshenkoLinearBeamElement3D2N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangian2D3N", mTotalLagrangian2D3N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangian2D4N", mTotalLagrangian2D4N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangian2D6N", mTotalLagrangian2D6N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangian2D8N", mTotalLagrangian2D8N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangian3D4N", mTotalLagrangian3D4N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangian3D10N", mTotalLagrangian3D10N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangian3D6N", mTotalLagrangian3D6N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangian3D15N", mTotalLagrangian3D15N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangian3D8N", mTotalLagrangian3D8N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangian3D20N", mTotalLagrangian3D20N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangian3D27N", mTotalLagrangian3D27N )

    KRATOS_REGISTER_ELEMENT( "KinematicLinear2D3N", mKinematicLinear2D3N )
    KRATOS_REGISTER_ELEMENT( "KinematicLinear2D4N", mKinematicLinear2D4N )
    KRATOS_REGISTER_ELEMENT( "KinematicLinear2D6N", mKinematicLinear2D6N )
    KRATOS_REGISTER_ELEMENT( "KinematicLinear2D8N", mKinematicLinear2D8N )
    KRATOS_REGISTER_ELEMENT( "KinematicLinear2D9N", mKinematicLinear2D9N )
    KRATOS_REGISTER_ELEMENT( "KinematicLinear3D4N", mKinematicLinear3D4N )
    KRATOS_REGISTER_ELEMENT( "KinematicLinear3D10N", mKinematicLinear3D10N )
    KRATOS_REGISTER_ELEMENT( "KinematicLinear3D8N", mKinematicLinear3D8N )
    KRATOS_REGISTER_ELEMENT( "KinematicLinear3D20N", mKinematicLinear3D20N )
    KRATOS_REGISTER_ELEMENT( "KinematicLinear3D27N", mKinematicLinear3D27N )
    KRATOS_REGISTER_ELEMENT( "KinematicLinear3D6N", mKinematicLinear3D6N )
    KRATOS_REGISTER_ELEMENT( "KinematicLinear3D15N", mKinematicLinear3D15N )
    #else
    KRATOS_REGISTER_ELEMENT( "CrisfieldTrussElement3D2N", mCrisfieldTrussElement3D2N )
    KRATOS_REGISTER_ELEMENT( "CrisfieldTrussElement3D3N", mCrisfieldTrussElement3D3N )
    KRATOS_REGISTER_ELEMENT( "TrussElement3D2N", mTrussElement3D2N )
    KRATOS_REGISTER_ELEMENT( "TrussElement3D3N", mTrussElement3D3N )
    //KRATOS_REGISTER_ELEMENT( "LinearIncompresibleElement2D3N", mLinearIncompresibleElement2D3N )
    //KRATOS_REGISTER_ELEMENT( "LinearIncompresibleElement3D4N", mLinearIncompresibleElement3D4N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangian", mTotalLagrangian3D4N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangian2D3N", mTotalLagrangian2D3N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangian2D4N", mTotalLagrangian2D4N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangian2D6N", mTotalLagrangian2D6N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangian2D8N", mTotalLagrangian2D8N )
    KRATOS_REGISTER_ELEMENT( "BeamElement3D2N", mBeamElement3D2N )
    KRATOS_REGISTER_ELEMENT( "BeamElement3D3N", mBeamElement3D3N )
    KRATOS_REGISTER_ELEMENT( "TimoshenkoBeamElement3D2N", mTimoshenkoBeamElement3D2N )
    KRATOS_REGISTER_ELEMENT( "TimoshenkoBeamElement3D3N", mTimoshenkoBeamElement3D3N )
    KRATOS_REGISTER_ELEMENT( "CorotationalLinearBeamElement2D2N", mCorotationalLinearBeamElement2D2N )
    KRATOS_REGISTER_ELEMENT( "CorotationalLinearBeamElement3D2N", mCorotationalLinearBeamElement3D2N )
    KRATOS_REGISTER_ELEMENT( "TimoshenkoLinearBeamElement2D2N", mTimoshenkoLinearBeamElement2D2N )
    KRATOS_REGISTER_ELEMENT( "TimoshenkoLinearBeamElement3D2N", mTimoshenkoLinearBeamElement3D2N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangian3D4N", mTotalLagrangian3D4N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangian3D10N", mTotalLagrangian3D10N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangian3D6N", mTotalLagrangian3D6N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangian3D15N", mTotalLagrangian3D15N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangian3D8N", mTotalLagrangian3D8N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangian3D20N", mTotalLagrangian3D20N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangian3D27N", mTotalLagrangian3D27N )

    KRATOS_REGISTER_ELEMENT( "MixedLagrangian3D4N", mMixedLagrangian3D4N )
    KRATOS_REGISTER_ELEMENT( "MixedLagrangian3D10N", mMixedLagrangian3D10N )
    KRATOS_REGISTER_ELEMENT( "MixedLagrangian3D6N", mMixedLagrangian3D6N )
    KRATOS_REGISTER_ELEMENT( "MixedLagrangian3D15N", mMixedLagrangian3D15N )
    KRATOS_REGISTER_ELEMENT( "MixedLagrangian3D8N", mMixedLagrangian3D8N )
    KRATOS_REGISTER_ELEMENT( "MixedLagrangian3D20N", mMixedLagrangian3D20N )
    KRATOS_REGISTER_ELEMENT( "MixedLagrangian3D27N", mMixedLagrangian3D27N )


    KRATOS_REGISTER_ELEMENT( "KinematicLinear2D3N", mKinematicLinear2D3N )
    KRATOS_REGISTER_ELEMENT( "KinematicLinear2D4N", mKinematicLinear2D4N )
    KRATOS_REGISTER_ELEMENT( "KinematicLinear2D6N", mKinematicLinear2D6N )
    KRATOS_REGISTER_ELEMENT( "KinematicLinear2D8N", mKinematicLinear2D8N )
    KRATOS_REGISTER_ELEMENT( "KinematicLinear2D9N", mKinematicLinear2D9N )
    KRATOS_REGISTER_ELEMENT( "KinematicLinear3D4N", mKinematicLinear3D4N )
    KRATOS_REGISTER_ELEMENT( "KinematicLinear3D10N", mKinematicLinear3D10N )
    KRATOS_REGISTER_ELEMENT( "KinematicLinear3D8N", mKinematicLinear3D8N )
    KRATOS_REGISTER_ELEMENT( "KinematicLinear3D20N", mKinematicLinear3D20N )
    KRATOS_REGISTER_ELEMENT( "KinematicLinear3D27N", mKinematicLinear3D27N )
    KRATOS_REGISTER_ELEMENT( "KinematicLinear3D6N", mKinematicLinear3D6N )
    KRATOS_REGISTER_ELEMENT( "KinematicLinear3D15N", mKinematicLinear3D15N )
    KRATOS_REGISTER_ELEMENT( "MembraneElement", mMembraneElement )
    KRATOS_REGISTER_ELEMENT( "IsoShellElement", mIsoShellElement )
    KRATOS_REGISTER_ELEMENT( "AnisoShellElement", mAnisoShellElement )
    // KRATOS_REGISTER_ELEMENT( "AnisoLinearShellElement", mAnisoLinearShellElement )
    KRATOS_REGISTER_ELEMENT( "UnsaturatedSoilsElement2PhaseSmallStrain3D4N", mUnsaturatedSoilsElement2PhaseSmallStrain3D4N )
    KRATOS_REGISTER_ELEMENT( "UnsaturatedSoilsElement2PhaseSmallStrain3D10N", mUnsaturatedSoilsElement2PhaseSmallStrain3D10N )
    KRATOS_REGISTER_ELEMENT( "UnsaturatedSoilsElement2PhaseSmallStrain3D20N", mUnsaturatedSoilsElement2PhaseSmallStrain3D20N )
    KRATOS_REGISTER_ELEMENT( "UnsaturatedSoilsElement2PhaseSmallStrain3D27N", mUnsaturatedSoilsElement2PhaseSmallStrain3D27N )
    KRATOS_REGISTER_ELEMENT( "UnsaturatedSoilsElement2PhaseSmallStrain3D15N", mUnsaturatedSoilsElement2PhaseSmallStrain3D15N )
    KRATOS_REGISTER_ELEMENT( "UnsaturatedSoilsElement2PhaseSmallStrain3D8N", mUnsaturatedSoilsElement2PhaseSmallStrain3D8N )
    KRATOS_REGISTER_ELEMENT( "UnsaturatedSoilsElement2PhaseSmallStrainStaggered3D27N", mUnsaturatedSoilsElement2PhaseSmallStrainStaggered3D27N )
    KRATOS_REGISTER_ELEMENT( "UnsaturatedSoilsElement3PhaseSmallStrain3D4N", mUnsaturatedSoilsElement3PhaseSmallStrain3D4N )
    KRATOS_REGISTER_ELEMENT( "UnsaturatedSoilsElement3PhaseSmallStrain3D10N", mUnsaturatedSoilsElement3PhaseSmallStrain3D10N )
    KRATOS_REGISTER_ELEMENT( "UnsaturatedSoilsElement3PhaseSmallStrain3D20N", mUnsaturatedSoilsElement3PhaseSmallStrain3D20N )
    KRATOS_REGISTER_ELEMENT( "UnsaturatedSoilsElement3PhaseSmallStrain3D27N", mUnsaturatedSoilsElement3PhaseSmallStrain3D27N )
    KRATOS_REGISTER_ELEMENT( "UnsaturatedSoilsElement3PhaseSmallStrain3D15N", mUnsaturatedSoilsElement3PhaseSmallStrain3D15N )
    KRATOS_REGISTER_ELEMENT( "UnsaturatedSoilsElement3PhaseSmallStrain3D8N", mUnsaturatedSoilsElement3PhaseSmallStrain3D8N )
    KRATOS_REGISTER_ELEMENT( "UnsaturatedSoilsElement3PhaseSmallStrainLiakopolous3D20N", mUnsaturatedSoilsElement3PhaseSmallStrainLiakopolous3D20N )
    KRATOS_REGISTER_ELEMENT( "UnsaturatedSoilsElement3PhaseSmallStrainLiakopolous3D27N", mUnsaturatedSoilsElement3PhaseSmallStrainLiakopolous3D27N )
    KRATOS_REGISTER_ELEMENT( "Ebst3D3N", mEbst3D3N )
    KRATOS_REGISTER_ELEMENT( "EbstVel3D3N", mEbstVel3D3N )
    KRATOS_REGISTER_ELEMENT( "EASElementQ4E4", mEASElementQ4E4 )

    KRATOS_REGISTER_CONDITION( "Face2D", mFace2D )
    KRATOS_REGISTER_CONDITION( "Face3D", mFace3D3N )
    KRATOS_REGISTER_CONDITION( "Face3D3N", mFace3D3N )
    KRATOS_REGISTER_CONDITION( "Face3D6N", mFace3D6N )
    KRATOS_REGISTER_CONDITION( "Face3D4N", mFace3D4N )
    KRATOS_REGISTER_CONDITION( "Face3D8N", mFace3D8N )
    KRATOS_REGISTER_CONDITION( "Face3D9N", mFace3D9N )
    KRATOS_REGISTER_CONDITION( "FacePressure3D3N", mFacePressure3D3N )
    KRATOS_REGISTER_CONDITION( "FacePressure3D6N", mFacePressure3D6N )
    KRATOS_REGISTER_CONDITION( "FacePressure3D4N", mFacePressure3D4N )
    KRATOS_REGISTER_CONDITION( "FacePressure3D8N", mFacePressure3D8N )
    KRATOS_REGISTER_CONDITION( "FacePressure3D9N", mFacePressure3D9N )
    KRATOS_REGISTER_CONDITION( "FacePressureTotalLagrangian3D3N", mFacePressureTotalLagrangian3D3N )
    KRATOS_REGISTER_CONDITION( "FacePressureTotalLagrangian3D6N", mFacePressureTotalLagrangian3D6N )
    KRATOS_REGISTER_CONDITION( "FacePressureTotalLagrangian3D4N", mFacePressureTotalLagrangian3D4N )
    KRATOS_REGISTER_CONDITION( "FacePressureTotalLagrangian3D8N", mFacePressureTotalLagrangian3D8N )
    KRATOS_REGISTER_CONDITION( "FacePressureTotalLagrangian3D9N", mFacePressureTotalLagrangian3D9N )
    KRATOS_REGISTER_CONDITION( "LineForce2D2N", mLineForce2D2N )
    KRATOS_REGISTER_CONDITION( "LineForce2D3N", mLineForce2D3N )
    KRATOS_REGISTER_CONDITION( "LineForce3D2N", mLineForce3D2N )
    KRATOS_REGISTER_CONDITION( "LineForce3D3N", mLineForce3D3N )
    KRATOS_REGISTER_CONDITION( "LinePressure2D2N", mLinePressure2D2N )
    KRATOS_REGISTER_CONDITION( "LinePressure2D3N", mLinePressure2D3N )
    KRATOS_REGISTER_CONDITION( "FaceForce3D3N", mFaceForce3D3N )
    KRATOS_REGISTER_CONDITION( "FaceForce3D6N", mFaceForce3D6N )
    KRATOS_REGISTER_CONDITION( "FaceForce3D4N", mFaceForce3D4N )
    KRATOS_REGISTER_CONDITION( "FaceForce3D8N", mFaceForce3D8N )
    KRATOS_REGISTER_CONDITION( "FaceForce3D9N", mFaceForce3D9N )
    KRATOS_REGISTER_CONDITION( "MasterContactFace3D", mMasterContactFace3D )
    KRATOS_REGISTER_CONDITION( "MasterContactFace3D3", mMasterContactFace3D3 )
    KRATOS_REGISTER_CONDITION( "MasterContactFace3D6", mMasterContactFace3D6 )
    KRATOS_REGISTER_CONDITION( "MasterContactFace3D8", mMasterContactFace3D8 )
    KRATOS_REGISTER_CONDITION( "MasterContactFace3D9", mMasterContactFace3D9 )
    KRATOS_REGISTER_CONDITION( "MasterContactFace3DNewmark", mMasterContactFace3DNewmark )
    KRATOS_REGISTER_CONDITION( "MasterContactFace3D3Newmark", mMasterContactFace3D3Newmark )
    KRATOS_REGISTER_CONDITION( "MasterContactFace3D6Newmark", mMasterContactFace3D6Newmark )
    KRATOS_REGISTER_CONDITION( "MasterContactFace3D8Newmark", mMasterContactFace3D8Newmark )
    KRATOS_REGISTER_CONDITION( "MasterContactFace3D9Newmark", mMasterContactFace3D9Newmark )
    KRATOS_REGISTER_CONDITION( "SlaveContactFace3D", mSlaveContactFace3D )
    KRATOS_REGISTER_CONDITION( "SlaveContactFace3D3", mSlaveContactFace3D3 )
    KRATOS_REGISTER_CONDITION( "SlaveContactFace3D6", mSlaveContactFace3D6 )
    KRATOS_REGISTER_CONDITION( "SlaveContactFace3D8", mSlaveContactFace3D8 )
    KRATOS_REGISTER_CONDITION( "SlaveContactFace3D9", mSlaveContactFace3D9 )
    KRATOS_REGISTER_CONDITION( "SlaveContactFace3DNewmark", mSlaveContactFace3DNewmark )
    KRATOS_REGISTER_CONDITION( "SlaveContactFace3D3Newmark", mSlaveContactFace3D3Newmark )
    KRATOS_REGISTER_CONDITION( "SlaveContactFace3D6Newmark", mSlaveContactFace3D6Newmark )
    KRATOS_REGISTER_CONDITION( "SlaveContactFace3D8Newmark", mSlaveContactFace3D8Newmark )
    KRATOS_REGISTER_CONDITION( "SlaveContactFace3D9Newmark", mSlaveContactFace3D9Newmark )
    KRATOS_REGISTER_CONDITION( "PointForce3D", mPointForce3D )
    KRATOS_REGISTER_CONDITION( "PointMoment3D", mPointMoment3D )

    KRATOS_REGISTER_CONDITION( "ElasticPointConstraint", mElasticPointConstraint )
    KRATOS_REGISTER_CONDITION( "ElasticLineConstraint2N", mElasticLineConstraint2N )
    KRATOS_REGISTER_CONDITION( "ElasticLineConstraint3N", mElasticLineConstraint3N )
    KRATOS_REGISTER_CONDITION( "ElasticFaceConstraint3N", mElasticFaceConstraint3N )
    KRATOS_REGISTER_CONDITION( "ElasticFaceConstraint6N", mElasticFaceConstraint6N )
    KRATOS_REGISTER_CONDITION( "ElasticFaceConstraint4N", mElasticFaceConstraint4N )
    KRATOS_REGISTER_CONDITION( "ElasticFaceConstraint8N", mElasticFaceConstraint8N )
    KRATOS_REGISTER_CONDITION( "ElasticFaceConstraint9N", mElasticFaceConstraint9N )

    KRATOS_REGISTER_CONDITION( "ElasticFaceSprings3N", mElasticFaceSprings3N )
    KRATOS_REGISTER_CONDITION( "ElasticFaceSprings6N", mElasticFaceSprings6N )
    KRATOS_REGISTER_CONDITION( "ElasticFaceSprings4N", mElasticFaceSprings4N )
    KRATOS_REGISTER_CONDITION( "ElasticFaceSprings8N", mElasticFaceSprings8N )
    KRATOS_REGISTER_CONDITION( "ElasticFaceSprings9N", mElasticFaceSprings9N )

    KRATOS_REGISTER_CONDITION( "NitscheIsotropicConstraint2D2N", mNitscheIsotropicConstraint2D2N )
    KRATOS_REGISTER_CONDITION( "NitscheIsotropicConstraint2D3N", mNitscheIsotropicConstraint2D3N )
    KRATOS_REGISTER_CONDITION( "NitscheIsotropicConstraint3D3N", mNitscheIsotropicConstraint3D3N )
    KRATOS_REGISTER_CONDITION( "NitscheIsotropicConstraint3D6N", mNitscheIsotropicConstraint3D6N )
    KRATOS_REGISTER_CONDITION( "NitscheIsotropicConstraint3D4N", mNitscheIsotropicConstraint3D4N )
    KRATOS_REGISTER_CONDITION( "NitscheIsotropicConstraint3D8N", mNitscheIsotropicConstraint3D8N )
    KRATOS_REGISTER_CONDITION( "NitscheIsotropicConstraint3D9N", mNitscheIsotropicConstraint3D9N )

    KRATOS_REGISTER_CONDITION( "NodeTyingLagrange", mNodeTyingLagrange )
    KRATOS_REGISTER_CONDITION( "NodeTyingLagrangeZ", mNodeTyingLagrangeZ )
    KRATOS_REGISTER_CONDITION( "PointPointJointCondition", mPointPointJointCondition )
    KRATOS_REGISTER_CONDITION( "PointPointLagrangeCondition", mPointPointLagrangeCondition )
    KRATOS_REGISTER_CONDITION( "FaceVel3D3N", mFaceVel3D3N )


//         KRATOS_REGISTER_ELEMENT("UPCTestElement3D20N", mUPCTestElement3D20N)
    KRATOS_REGISTER_CONDITION( "PointForce2D", mPointForce2D )


    KRATOS_REGISTER_CONDITION( "SlaveContactPoint2D", mSlaveContactPoint2D )
    KRATOS_REGISTER_CONDITION( "MasterContactFace2D", mMasterContactFace2D )

// KratosComponents<ConstitutiveLaw >::Add("Isotropic3D", mIsotropic3D);
    Serializer::Register( "Isotropic3D", mIsotropic3D );
    Serializer::Register( "DummyConstitutiveLaw", mDummyConstitutiveLaw );
    Serializer::Register( "DruckerPrager", mDruckerPrager );
    Serializer::Register( "CamClay3D", mCamClay3D );
//        std::cout << "registered objects:" << std::endl;
//
//        for(Serializer::RegisteredObjectsContainerType::iterator i = Serializer::GetRegisteredObjects().begin() ; i != Serializer::GetRegisteredObjects().end() ; i++)
//            std::cout << i->first << std::endl;
    #endif
}

}  // namespace Kratos.



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
//   Date:                $Date: 2012-09-17 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_PILE_UTILITY_INCLUDED )
#define  KRATOS_PILE_UTILITY_INCLUDED

// System includes

// External includes
#include "boost/smart_ptr.hpp"
#include "boost/timer.hpp"
// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "custom_conditions/pile_kinematic_linear.h"
#include "utilities/math_utils.h"
#include "geometries/point_3d.h"


namespace Kratos
{
/**
    * Steering Utility
        */

class PileUtility
{
public:
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;
    typedef Geometry<Node<3> >::IntegrationPointsArrayType IntegrationPointsArrayType;
    typedef Geometry<Node<3> > GeometryType;
    typedef Properties PropertiesType;

    /**
     * class pointer definition
     */
    KRATOS_CLASS_POINTER_DEFINITION( PileUtility );

    /**
     * Constructor.
     */
    PileUtility()
    {
    }

    /**
     * Destructor.
     */
    virtual ~PileUtility()
    {
    }

    /**
     * Initializes mesh tying by means of lagrange multipliers
     */
    void InitializePileUtility( ModelPart& model_part, std::vector<unsigned int>& pile_elements, std::vector<unsigned int>& soil_elements )
    {
        ElementsArrayType::Pointer piles( new ElementsArrayType() );
        ElementsArrayType::Pointer soil_elems( new ElementsArrayType() );


        Node<3> point( 0.0, 0.0, 0.0 );
        GeometryType::Pointer tempGeometry = GeometryType::Pointer(new Point3D<Node<3> >(point) );
//        KRATOS_WATCH( tempGeometry );
//        KRATOS_WATCH( *tempGeometry );

        int properties_index = model_part.NumberOfProperties();
        PropertiesType::Pointer tempProperties( new PropertiesType( properties_index + 1 ) );
        model_part.AddProperties( tempProperties );


        std::cout << "Initializing PileUtility..." << std::endl;

        for ( unsigned int it = 0; it != pile_elements.size(); it++ )
        {
            piles->push_back( model_part.GetElement( pile_elements[it] ) );
        }

        for ( unsigned int it = 0; it != soil_elements.size(); it++ )
        {
            soil_elems->push_back( model_part.GetElement( soil_elements[it] ) );
        }

        for ( ElementsArrayType::ptr_iterator it = piles->ptr_begin();
                it != piles->ptr_end(); ++it )
        {
            /******KRATOS_WATCH(it);*/
//            KRATOS_WATCH( *it );

            for ( IndexType i = 0; i < ( *it )->GetGeometry().IntegrationPoints().size(); i++ )
            {
                Point<3> PilePoint;
                Point<3> PileLocalPoint = ( *it )->GetGeometry().IntegrationPoints()[i];
                PilePoint = ( *it )->GetGeometry().GlobalCoordinates( PilePoint, PileLocalPoint );
                Point<3> SoilLocalPoint;
                Element::Pointer TargetElement;

                if ( FindPartnerElement( PilePoint, soil_elems, TargetElement, SoilLocalPoint ) )
                {

                    Point<3> SoilGlobalPoint;
                    SoilGlobalPoint = TargetElement->GetGeometry().GlobalCoordinates( SoilGlobalPoint, SoilLocalPoint );

                    int a = ( model_part.Conditions().end() - 1 )->Id();

                    IndexType newId = a + 1;


                    Condition::Pointer newLink = Condition::Pointer( new Pile_Kinematic_Linear( newId, tempGeometry, tempProperties, TargetElement, *it,
                   // Condition::Pointer newLink = Condition::Pointer( new PileCondition( newId, tempGeometry, tempProperties, TargetElement, *it,
                                                 SoilLocalPoint, PileLocalPoint, i ) );

                    model_part.Conditions().push_back( newLink );

                }

            }
        }

        return;
    }//InitializePileUtility

    /**
     * calculates for a point given with the physical coords newNode
     * the element oldElement where it lays in and the natural coords
     * localPoint within this element
     * @return whether a corresponding element and natural coords could be found
     * @param newNode physical coordinates of given point
     * @param OldMeshElementsArray Array of elements wherein the search should be performed
     * @param oldElement corresponding element for newNode
     * @param rResult corresponding natural coords for newNode
     * TODO: find a faster method for outside search (hextree? etc.), maybe outside this
     * function by restriction of OldMeshElementsArray
     */
    bool FindPartnerElement( Point<3>& sourcePoint,
                             const ElementsArrayType::Pointer& SoilElements,
                             Element::Pointer& TargetElement, Point<3>& rResult )
    {
        bool partner_found = false;
        ElementsArrayType::Pointer SoilElementsCandidates( new ElementsArrayType() );
        std::vector<double > OldMinDist;
        bool newMinDistFound = false;

        std::cout << "line 187" << std::endl;
        int counter = 0;

        do
        {
            double minDist = 1.0e120;
            newMinDistFound = false;
            SoilElementsCandidates->clear();
            // (global search)

            for ( ElementsArrayType::ptr_iterator it = SoilElements->ptr_begin();
                    it != SoilElements->ptr_end(); ++it )
            {
                //loop over all nodes in tested element
                for ( unsigned int n = 0; n < ( *it )->GetGeometry().size(); n++ )
                {
                    double dist = (( *it )->GetGeometry().GetPoint( n ).X0() - sourcePoint[0] )
                                  * (( *it )->GetGeometry().GetPoint( n ).X0() - sourcePoint[0] )
                                  + (( *it )->GetGeometry().GetPoint( n ).Y0() - sourcePoint[1] )
                                  * (( *it )->GetGeometry().GetPoint( n ).Y0() - sourcePoint[1] )
                                  + (( *it )->GetGeometry().GetPoint( n ).Z0() - sourcePoint[2] )
                                  * (( *it )->GetGeometry().GetPoint( n ).Z0() - sourcePoint[2] );

                    if ( fabs( dist - minDist ) < 1e-7 )
                    {
                        SoilElementsCandidates->push_back( *it );
                    }
                    else if ( dist < minDist )
                    {
                        bool alreadyUsed = false;

                        for ( unsigned int old_dist = 0; old_dist < OldMinDist.size(); old_dist++ )
                        {
                            if ( fabs( dist - OldMinDist[old_dist] ) < 1e-7 )
                                alreadyUsed = true;
                        }

                        if ( !alreadyUsed )
                        {
                            SoilElementsCandidates->clear();
                            minDist = dist;
                            SoilElementsCandidates->push_back( *it );
                            newMinDistFound = true;
                        }
                    }
                }
            }

            OldMinDist.push_back( minDist );

//                     KRATOS_WATCH(OldElementsSet->size());

            for ( ElementsArrayType::ptr_iterator it = SoilElementsCandidates->ptr_begin();
                    it != SoilElementsCandidates->ptr_end(); ++it )
            {
//                         std::cout << "checking elements list" << std::endl;
                if (( *it )->GetGeometry().IsInside( sourcePoint, rResult ) )
                {
//                             std::cout << "isInside" << std::endl;
                    TargetElement = *( it );
                    partner_found = true;
                    return partner_found;
                }
            }

//                     std::cout << counter << std::endl;
            counter++;

            if ( counter > 27 )
                break;
        }
        while ( newMinDistFound );

        if ( !partner_found )
            std::cout << " !!!! NO PARTNER FOUND !!!! " << std::endl;

        return partner_found;
    }

    /**
     * Calculates for given Loacal coordinates the global coordinates
     * @param Surface surface
     * @param rResult global coordinates
     * @param LocalCoordinates local coordinates
     * @return global coordinates
     */

    GeometryType::CoordinatesArrayType& GlobalCoordinatesSoil( Element::Pointer SoilElements, GeometryType::CoordinatesArrayType& rResult, GeometryType::CoordinatesArrayType const& LocalCoordinates )
    {
        noalias( rResult ) = ZeroVector( 3 );

        for ( IndexType i = 0 ; i < SoilElements->GetGeometry().size() ; i++ )
        {
            double shape_func = SoilElements->GetGeometry().ShapeFunctionValue( i, LocalCoordinates );

            rResult( 0 ) += shape_func *
                            (( SoilElements->GetGeometry()[i] ).X0()
                             + ( SoilElements->GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_X ) );

            rResult( 1 ) += shape_func *
                            (( SoilElements->GetGeometry()[i] ).Y0()
                             + ( SoilElements->GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_Y ) );

            rResult( 2 ) += shape_func *
                            (( SoilElements->GetGeometry()[i] ).Z0()
                             + ( SoilElements->GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_Z ) );
        }

        return rResult;
    }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    GeometryType::CoordinatesArrayType& GlobalCoordinatesPile( const Element::Pointer PileElements, GeometryType::CoordinatesArrayType& rResult, GeometryType::CoordinatesArrayType const& LocalCoordinates )
    {
        noalias( rResult ) = ZeroVector( 3 );

        for ( IndexType i = 0 ; i < PileElements->GetGeometry().size() ; i++ )
        {
            double shape_func = PileElements->GetGeometry().ShapeFunctionValue( i, LocalCoordinates );

            rResult( 0 ) += shape_func *
                            (( PileElements->GetGeometry()[i] ).X0()
                             + ( PileElements->GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_X ) );

            rResult( 1 ) += shape_func *
                            (( PileElements->GetGeometry()[i] ).Y0()
                             + ( PileElements->GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_Y ) );

            rResult( 2 ) += shape_func *
                            (( PileElements->GetGeometry()[i] ).Z0()
                             + ( PileElements->GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_Z ) );
        }

        return rResult;
    }

};//class PileUtility
}  // namespace Kratos.

#endif // KRATOS_PILE_UTILITY_INCLUDED defined 

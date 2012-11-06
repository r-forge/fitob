/* @HEADER@ */
// ************************************************************************
//
//                              Fitob
//            Copyright (2012) Janos Benk (benkjanos@gmail.com)
//
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Janos Benk (benkjanos@gmail.com),
//
// ************************************************************************
/* @HEADER@ */

/*
 * FideumGaussSeidelWB.hpp
 *
 *  Created on: Jan 3, 2011
 *      Author: benk
 */

#ifndef FIDEUMGAUSSSEIDEL_WB_HPP_
#define FIDEUMGAUSSSEIDEL_WB_HPP_

#include "src/utils/fitobdefs.hpp"
#include "src/pdesolver/FitobSmootherBase.hpp"

namespace fitob {

using namespace std;

/**
 * Smoother which operates on grids without boundary points <br>.
 * In addition this smoother assumes second derivative = 0 boundary condition an all boundaries <br>.
 * This boundary condition is realized with ghost points which satisfy the second derivative = 0 condition */
class GaussSeidel_WB : public SmootherBase{
public:

	GaussSeidel_WB( const XMLConfiguration* config , MultigridFGBase *grid ,
			     const ModelCollection* models ,
			     const FitobCalculator* fitobcalculator);

	virtual ~GaussSeidel_WB(){;}

	/** sets the new time step and signals that diag and rhs should be computed newly*/
	virtual void setTimeStep( double microTimestep ) {
		isDiagonalCalculated_ = false;
		calculateRHS_ = false;
		microTimestep_ = microTimestep;
	}

	/** this method is used in multigrid to signal that it should not calculate the residuum */
	void residumIsSet() { calculateRHS_ = true; }

	/** */
	void smoothGrid( const DVector& globalCoords);

	/** calculates the RHS */
	void calculateRHS( DVector& globCoord);

	/** makes one explicit midpoint step , see super class */
	void explicitStep( DVector& globCoord , DVector& unknowns , double underrelaxFactor );

	/** */
    double calcResiduum( const DVector& globalCoords );

private:

	void calculateDiag( DVector& globCoord);

	void oneSmoothGrid( DVector& globCoord);

	/** function notation from the formula 1.10 */
	inline double V_ij_x_sp(int ind , double hh0 , double hh1 , DVector &gridV) {
		//return ((hh0/hh1)*grid()->u()[linInd+offs] - (hh1/hh0)*grid()->u()[linInd-offs])/(hh0+hh1);
		return ((hh0/hh1)*gridV[3*ind+2] - (hh1/hh0)*gridV[3*ind])/(hh0+hh1);
	}

	/** function notation from the formula 1.10 */
	inline double V_ij_x(int ind , double hh0 , double hh1 , double aa_x_ij , DVector &gridV) {
		//return V_ij_x_sp( linInd , offs , hh0 , hh1 ) + aa_x_ij*grid()->u()[linInd];
		return V_ij_x_sp( ind , hh0 , hh1 , gridV) + aa_x_ij*gridV[3*ind+1];
	}

	/** function notation from the formula 1.13 */
	inline double V_ij_xx_sp( int ind , double hh0 , double hh1 , DVector &gridV) {
		//	return ( hh1*grid()->u()[linInd-offs] + hh0*grid()->u()[linInd+offs] ) / (0.5*hh0*hh1*(hh1+hh0));
		return ( hh1*gridV[3*ind] + hh0*gridV[3*ind+2] ) / (0.5*hh0*hh1*(hh1+hh0));
	}

	/** function notation from the formula 1.13 */
	inline double V_ij_xx(int ind , double hh0 , double hh1 , double aa_xx_ij , DVector &gridV) {
		//	return V_ij_xx_sp( linInd , offs , hh0 , hh1 ) + aa_xx_ij*grid()->u()[linInd];
		return V_ij_xx_sp( ind , hh0 , hh1 , gridV) + aa_xx_ij*gridV[3*ind+1];
	}

	/** Update the status variables at one point*/
	inline void updateStage(int tmp_I , DVector& globCoord , int &linearIndex , IVector& axisIndex , DVector &gridV ,
			                 IVector &tmpMark , DVector &h0 , DVector& h1 , DVector& a_x_ij, DVector& a_xx_ij){

	   //FITOB_OUT_LEVEL3(verb()," GaussSeidl::updateStage  ======================" );
	   linearIndex = 0;
	   //FITOB_OUT_LEVEL3(verb()," GaussSeidl::updateStage , tmp_I:" << tmp_I);
	   for (int i = gdim_-1 ; i >= 0 ; i--){
		   axisIndex[i] = tmp_I / (innerOffset[i]);
		   tmp_I = tmp_I % innerOffset[i];
		   //FITOB_OUT_LEVEL3(verb()," GaussSeidl::updateStage: i:"<<i<<",axisIndex[i]:" << axisIndex[i]
           //                  << " innerOffset[i]:" << innerOffset[i] << ", tmp_I:" << tmp_I);
		   // here we calculate the linear index on the grid
		   linearIndex = linearIndex + (axisIndex[i])*grid()->offset()[i];
	       // set the global coordinates
		   globCoord[ grid()->domain()->localToGlobalIndex(i) ] = grid()->axisScaling(i)[axisIndex[i]+1];
	   }
	   // calculate the "h"
           #pragma ivdep
	   for (int i = 0 ; i < nrFactors_ ; i++){
		    int aI =  factorIndex_to_local[i];
			h0[i] = grid()->axisScaling(aI)[axisIndex[aI]+1] - grid()->axisScaling(aI)[axisIndex[aI]];
			h1[i] = grid()->axisScaling(aI)[axisIndex[aI]+2] - grid()->axisScaling(aI)[axisIndex[aI]+1];
			a_xx_ij[i] = - 1.0 / (0.5*h0[i]*h1[i]);
			// update gridV_
			//FITOB_OUT_LEVEL3(verb()," GaussSeidl::updateStage: i:"<<i<<" linearIndex:"<<linearIndex<< ",innerOffset[i]:" << innerOffset[i]);
			//FITOB_OUT_LEVEL3(verb()," GaussSeidl::updateStage: axisIndex[i]:" << axisIndex[i] );
			tmpMark[3*i] = -1; tmpMark[3*i+1] = 0; tmpMark[3*i+2] = -1;
			gridV[3*i+1] = grid()->u()[linearIndex];
			a_x_ij[i] = ( h1[i]/h0[i] - h0[i]/h1[i])/(h0[i]+h1[i]);
            if ( axisIndex[i] > 0) { gridV[3*i] = grid()->u()[linearIndex-innerOffset[i]]; tmpMark[3*i] = 0;}
            if ( axisIndex[i] < grid()->nrPoints(i) - 1 ) { gridV[3*i+2] = grid()->u()[linearIndex+innerOffset[i]]; tmpMark[3*i+2] = 0;}
            // set the extrapolated values
            gridV[3*i]   = (tmpMark[3*i] == 0) ? gridV[3*i] : gridV[3*i+1] - (h0[i]*(gridV[3*i+2]-gridV[3*i+1])/h1[i]);
            gridV[3*i+2] = (tmpMark[3*i+2] == 0) ? gridV[3*i+2] : gridV[3*i+1] + (h1[i]*(gridV[3*i+1]-gridV[3*i])/h0[i]);
            //FITOB_OUT_LEVEL3(verb()," GaussSeidl::updateStage: i:" << i << " = [" << gridV_[3*i] << ","<<gridV_[3*i+1]<<","<<gridV_[3*i+2]<<"]" );
            //FITOB_OUT_LEVEL3(verb()," GaussSeidl::updateStage: i:" << i << " = [" << tmpMark_[3*i] << ","<<tmpMark_[3*i+1]<<","<<tmpMark_[3*i+2]<<"]" );
			// gridV_xy_ needs to updated separately for each x and y combination
	   }
	   //FITOB_OUT_LEVEL3(verb()," GaussSeidl::updateStage  linearIndex: " << linearIndex );
	}

	// =================== correlation =======================================

		/** function notation from the formula 1.10 */
		inline double V1_ij_x_sp(int ind , int offs , double hh0 , double hh1 , DVector &gridV_xy) {
			//return ((hh0/hh1)*grid()->u()[linInd+offs] - (hh1/hh0)*grid()->u()[linInd-offs])/(hh0+hh1);
			return ((hh0/hh1)*gridV_xy[offs + 2] - (hh1/hh0)*gridV_xy[offs])/(hh0+hh1);
		}

		/** function notation from the formula 1.10 */
		inline double V1_ij_x(int ind , int offs , double hh0 , double hh1 , double aa_x_ij , DVector &gridV_xy) {
			//return V_ij_x_sp( linInd , offs , hh0 , hh1 ) + aa_x_ij*grid()->u()[linInd];
			return V1_ij_x_sp( ind , offs , hh0 , hh1 , gridV_xy) + aa_x_ij*gridV_xy[offs+1];
		}

		/** function notation from the formula 1.16 */
		inline double V_ij_xy_sp(int indx , int indy ,
				 double hh0x , double hh1x , double hh0y , double hh1y ,
				 double aa_x_ij , double aa_y_ij , DVector &gridV_xy) {
	       //( (hh0x/hh1x)*V_ij_x( linInd + offsx, offsy , hh0y , hh1y , aa_y_ij )
	       //- (hh1x/hh0x)*V_ij_x( linInd - offsx, offsy , hh0y , hh1y , aa_y_ij )) / (hh0x + hh1x);
			return ((hh0x/hh1x)*V1_ij_x( indx, 6 , hh0y , hh1y , aa_y_ij , gridV_xy)
				  - (hh1x/hh0x)*V1_ij_x( indx, 0 , hh0y , hh1y , aa_y_ij , gridV_xy)) / (hh0x + hh1x);
		}

		/** function notation from the formula 1.16 */
		inline double V_ij_xy(int indx , int indy ,
				 double hh0x , double hh1x , double hh0y , double hh1y ,
				 double aa_x_ij , double aa_y_ij , DVector &gridV_xy) {
			        //V_ij_xy_sp( linInd , offsx , offsy , hh0x , hh1x , hh0y , hh1y , aa_x_ij , aa_y_ij )
					//+ aa_x_ij * V_ij_x_sp( linInd, offsy , hh0y , hh1y )
					//+ aa_x_ij * aa_y_ij * grid()->u()[linInd];
			return  V_ij_xy_sp( indx , indy , hh0x , hh1x , hh0y , hh1y , aa_x_ij , aa_y_ij , gridV_xy)
					+ aa_x_ij * V1_ij_x_sp( indy , 3 , hh0y , hh1y , gridV_xy)
					+ aa_x_ij * aa_y_ij * gridV_xy[4];
		}

	/** method to update the buffer for the correlation calculation */
	inline void updateStage_xy(int factIndX , int factIndY , DVector &gridV_xy , const DVector &gridV ,
			                   int &linearIndex , IVector &axisIndex , IVector &tmpMark_xy){
		//gridV_xy_ needs to updated separately for each x and y combination
		gridV_xy[4] = gridV[3*factIndX+1];
		gridV_xy[1] = gridV[3*factIndX];
		gridV_xy[7] = gridV[3*factIndX+2];
		gridV_xy[3] = gridV[3*factIndY];
		gridV_xy[5] = gridV[3*factIndY+2];

		tmpMark_xy[4] = 0 ; tmpMark_xy[3] = 0 ; tmpMark_xy[5] = 0 ; tmpMark_xy[1] = 0 ; tmpMark_xy[7] = 0;
		tmpMark_xy[0] = -1; tmpMark_xy[2] = -1; tmpMark_xy[6] = -1; tmpMark_xy[8] = -1;
		if ( (axisIndex[factIndX] > 0) && (axisIndex[factIndY] > 0)){
			gridV_xy[0] =  grid()->u()[linearIndex-innerOffset[factIndX]-innerOffset[factIndY]]; tmpMark_xy[0] = 0;  }
		if ( (axisIndex[factIndX] > 0) && (axisIndex[factIndY] < grid()->nrPoints(factIndY) - 1)){
			gridV_xy[2] =  grid()->u()[linearIndex-innerOffset[factIndX]+innerOffset[factIndY]]; tmpMark_xy[2] = 0;  }
		if ( (axisIndex[factIndX] < grid()->nrPoints(factIndX) - 1) && (axisIndex[factIndY] > 0)){
			gridV_xy[6] =  grid()->u()[linearIndex+innerOffset[factIndX]-innerOffset[factIndY]]; tmpMark_xy[6] = 0;  }
		if ( (axisIndex[factIndX] < grid()->nrPoints(factIndX) - 1) && (axisIndex[factIndY] < grid()->nrPoints(factIndY) - 1)){
			gridV_xy[8] =  grid()->u()[linearIndex+innerOffset[factIndX]+innerOffset[factIndY]]; tmpMark_xy[8] = 0;  }

		// extrapolation , if value is not there
		// make the interpolation in one cell, just as a plane (because those values are there for sure)
		// we do not use the h0 and h1 because we have already the differences in those direction (in X and Y direction)
		gridV_xy[0] = (tmpMark_xy[0] == 0) ? gridV_xy[0] : gridV_xy[1] - (gridV_xy[4]-gridV_xy[3]);//h0[factIndX];
		gridV_xy[2] = (tmpMark_xy[2] == 0) ? gridV_xy[2] : gridV_xy[1] + (gridV_xy[5]-gridV_xy[4]);//h1[factIndX];
		gridV_xy[6] = (tmpMark_xy[6] == 0) ? gridV_xy[6] : gridV_xy[7] - (gridV_xy[4]-gridV_xy[3]);//h0[factIndX];
		gridV_xy[8] = (tmpMark_xy[8] == 0) ? gridV_xy[8] : gridV_xy[7] + (gridV_xy[5]-gridV_xy[4]);//h1[factIndX];

		//FITOB_OUT_LEVEL3(verb(),"====================== updateStage_xy linearIndex:" << linearIndex);
        //FITOB_OUT_LEVEL3(verb(),"updateStage_xy " << " = [" << gridV_xy_[6] << ","<<gridV_xy_[7]<<","<<gridV_xy_[8]<<"]" );
        //FITOB_OUT_LEVEL3(verb(),"updateStage_xy " << " = [" << gridV_xy_[3] << ","<<gridV_xy_[4]<<","<<gridV_xy_[5]<<"]" );
        //FITOB_OUT_LEVEL3(verb(),"updateStage_xy " << " = [" << gridV_xy_[0] << ","<<gridV_xy_[1]<<","<<gridV_xy_[2]<<"]" );
		//FITOB_OUT_LEVEL3(verb(),"-----------------------------------------------------------" );
        //FITOB_OUT_LEVEL3(verb(),"updateStage_xy " << " = [" << tmpMark_xy[6] << ","<<tmpMark_xy[7]<<","<<tmpMark_xy[8]<<"]" );
        //FITOB_OUT_LEVEL3(verb(),"updateStage_xy " << " = [" << tmpMark_xy[3] << ","<<tmpMark_xy[4]<<","<<tmpMark_xy[5]<<"]" );
        //FITOB_OUT_LEVEL3(verb(),"updateStage_xy " << " = [" << tmpMark_xy[0] << ","<<tmpMark_xy[1]<<","<<tmpMark_xy[2]<<"]" );
	}

private:

	/** to save flops we store the diagonal entries of the matrix */
	bool isDiagonalCalculated_;

	/** flag to indicate if the right hand side is calculated */
	bool calculateRHS_;

	/** maps the factor index to the local variables which are the axis of the mesh*/
    IVector factorIndex_to_local;

    /** we store here the offsets, which are the same than those */
    IVector innerOffset;

    /** nr of grids points over which there should be smoothing done*/
	int innerPoints;

};

}

#endif /* FIDEUMGAUSSSEIDEL_WB_HPP_ */

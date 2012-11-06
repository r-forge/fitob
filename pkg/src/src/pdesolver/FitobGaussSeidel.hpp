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
 * FitobGaussSeidel.hpp
 *
 *  Created on: Jul 15, 2010
 *      Author: benk
 */

#ifndef FITOBGAUSSSEIDEL_HPP_
#define FITOBGAUSSSEIDEL_HPP_

#include "src/utils/fitobdefs.hpp"
#include "src/pdesolver/FitobSmootherBase.hpp"

namespace fitob {

  using namespace std;

  /** Classic Gauss Seidl smoother <br>
   * This smoother assumes Dirichlet boundary condition <br> */
class GaussSeidel : public SmootherBase{
public:

	GaussSeidel( const XMLConfiguration* config , MultigridFGBase *grid ,
			     const ModelCollection* models ,
			     const FitobCalculator* fitobcalculator);

	virtual ~GaussSeidel(){;}

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
	inline double V_ij_x_sp(int linInd , int offs , double hh0 , double hh1) {
		return ((hh0/hh1)*grid()->u()[linInd+offs] - (hh1/hh0)*grid()->u()[linInd-offs])/(hh0+hh1);
	}

	/** function notation from the formula 1.10 */
	inline double V_ij_x(int linInd , int offs , double hh0 , double hh1 , double aa_x_ij) {
		return V_ij_x_sp( linInd , offs , hh0 , hh1 ) + aa_x_ij*grid()->u()[linInd];
	}

	/** function notation from the formula 1.13 */
	inline double V_ij_xx_sp(int linInd , int offs , double hh0 , double hh1 ) {
		return ( hh1*grid()->u()[linInd-offs] + hh0*grid()->u()[linInd+offs] ) / (0.5*hh0*hh1*(hh1+hh0));
	}

	/** function notation from the formula 1.13 */
	inline double V_ij_xx(int linInd , int offs , double hh0 , double hh1 , double aa_xx_ij) {
		return V_ij_xx_sp( linInd , offs , hh0 , hh1 ) + aa_xx_ij*grid()->u()[linInd];
	}

	/** function notation from the formula 1.16 */
	inline double V_ij_xy_sp(int linInd , int offsx , int offsy ,
			 double hh0x , double hh1x , double hh0y , double hh1y ,
			 double aa_x_ij , double aa_y_ij ) {
       return ( (hh0x/hh1x)*V_ij_x( linInd + offsx, offsy , hh0y , hh1y , aa_y_ij )
    		     - (hh1x/hh0x)*V_ij_x( linInd - offsx, offsy , hh0y , hh1y , aa_y_ij )) / (hh0x + hh1x);
	}

	/** function notation from the formula 1.16 */
	inline double V_ij_xy(int linInd , int offsx , int offsy ,
			 double hh0x , double hh1x , double hh0y , double hh1y ,
			 double aa_x_ij , double aa_y_ij) {
		return V_ij_xy_sp( linInd , offsx , offsy , hh0x , hh1x , hh0y , hh1y , aa_x_ij , aa_y_ij )
				+ aa_x_ij * V_ij_x_sp( linInd, offsy , hh0y , hh1y )
				+ aa_x_ij * aa_y_ij * grid()->u()[linInd];
	}

	/** Update the status variables at one point*/
	inline void updateStage(int tmp_I , DVector& globCoord , int &linearIndex, IVector &axisIndex,
			                DVector &h0 , DVector &h1 , DVector &a_x_ij , DVector &a_xx_ij){
	   linearIndex = 0;
	   //FITOB_OUT_LEVEL3(verb()," GaussSeidl::updateStage , tmp_I:" << tmp_I);
	   for (int i = gdim_-1 ; i >= 0 ; i--){
		   axisIndex[i] = tmp_I / (innerOffset[i]);
		   tmp_I = tmp_I % innerOffset[i];
		   //FITOB_OUT_LEVEL3(verb()," GaussSeidl::updateStage: i:"<<i<<",axisIndex[i]:" << axisIndex[i]
           //                  << " innerOffset[i]:" << innerOffset[i] << ", tmp_I:" << tmp_I);
		   // here we calculate the linear index on the grid
		   linearIndex = linearIndex + (axisIndex[i]+1)*grid()->offset()[i];
	       // set the global coordinates
		   //FITOB_OUT_LEVEL3(verb()," GaussSeidl::updateStage: i:"<<i<<" grid()->axisScaling(i).size():" << grid()->axisScaling(i).size());
		   globCoord[ grid()->domain()->localToGlobalIndex(i) ] = grid()->axisScaling(i)[axisIndex[i]+1];
	   }
	   // calculate the "h"
	   for (int i = 0 ; i < nrFactors_ ; i++){
		    int aI =  factorIndex_to_local[i];
			h0[i] = grid()->axisScaling(aI)[axisIndex[aI]+1] - grid()->axisScaling(aI)[axisIndex[aI]];
			h1[i] = grid()->axisScaling(aI)[axisIndex[aI]+2] - grid()->axisScaling(aI)[axisIndex[aI]+1];
			a_x_ij[i] = ( h1[i]/h0[i] - h0[i]/h1[i])/(h0[i]+h1[i]);
			a_xx_ij[i] = - 1.0 / (0.5*h0[i]*h1[i]);
	   }
	}

	/** to save flops we store the diagonal entries of the matrix */
	bool isDiagonalCalculated_;

	/** flag to indicate if the right hand side is calculated */
	bool calculateRHS_;

    IVector factorIndex_to_local;
    IVector innerOffset;

    /** the total inner points */
	int innerPoints;

};

}

#endif /* FITOBGAUSSSEIDEL_HPP_ */

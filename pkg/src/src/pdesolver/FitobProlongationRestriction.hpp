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
 * FitobProlongationRestriction.hpp
 *
 *  Created on: Jul 6, 2010
 *      Author: benk
 */

#ifndef FITOBPROLONGATIONRESTRICTION_HPP_
#define FITOBPROLONGATIONRESTRICTION_HPP_

#include "src/utils/fitobdefs.hpp"
#include "src/pdesolver/FitobMultigridFGBase.hpp"

namespace fitob {

  using namespace std;

  /** */
  class ProlongationRestriction : public VerbClass{

  public:

	ProlongationRestriction();

	virtual ~ProlongationRestriction() {;}

	/** Function to make linear prolongation on the stretched grid, from the coarser to the finer grid
	 * @param gFine fine grid
	 * @param vFine the belonging vector to the fine grid
	 * @param gCoarse coarse grid
	 * @param vCoarse the belonging vector to the coarse grid
	 * @param coef the coefficient to be multiplied with */
	static void makeLinearProlongation( MultigridFGBase* gFine , DVector& vFine ,
			                            const MultigridFGBase* gCoarse , const DVector& vCoarse ,
			                            double coef , double coefOLD );

	/** Function to make direct restriction. <br> Direct means that we take 1:1 the
	 * values from the fine grid to the coarser grid.
	 * @param gFine fine grid
	 * @param vFine the belonging vector to the fine grid
	 * @param gCoarse coarse grid
	 * @param vCoarse the belonging vector to the coarse grid
	 * @param coef the coefficient to be multiplied with */
	static void makeDirectRestriction( const MultigridFGBase* gFine , const DVector& vFine ,
			                           MultigridFGBase* gCoarse , DVector& vCoarse ,
			                           double coef , double coefOLD );

  private :

	/** prolongation for grids with boundary points */
	static void makeLinearProlongation_withBoundary(
			                           MultigridFGBase* gFine , DVector& vFine ,
			                           const MultigridFGBase* gCoarse , const DVector& vCoarse ,
			                           double coef , double coefOLD );
	/** restriction for grids with boundary points */
	static void makeDirectRestriction_withBoundary(
			                           const MultigridFGBase* gFine , const DVector& vFine ,
			                           MultigridFGBase* gCoarse , DVector& vCoarse ,
			                           double coef , double coefOLD );
	/** prolongation for grids without boundary points */
	static void makeLinearProlongation_withoutBoundary(
			                           MultigridFGBase* gFine , DVector& vFine ,
			                           const MultigridFGBase* gCoarse , const DVector& vCoarse ,
			                           double coef , double coefOLD );
	/** restriction for grids without boundary points */
	static void makeDirectRestriction_withoutBoundary(
			                           const MultigridFGBase* gFine , const DVector& vFine ,
			                           MultigridFGBase* gCoarse , DVector& vCoarse ,
			                           double coef , double coefOLD );

};

}

#endif /* FITOBPROLONGATIONRESTRICTION_HPP_ */

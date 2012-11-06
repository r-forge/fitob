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
 * FitobSmootherBase.hpp
 *
 *  Created on: Jul 6, 2010
 *      Author: benk
 */

#ifndef FITOBSMOOTHERBASE_HPP_
#define FITOBSMOOTHERBASE_HPP_

#include "src/utils/fitobdefs.hpp"
#include "src/utils/FitobXMLConfiguration.hpp"
#include "src/pdesolver/FitobMultigridFGBase.hpp"
#include "src/diffusionmodel/FitobModelCollection.hpp"

namespace fitob {

   using namespace std;

   /** base class for smoothing algorithms */
   class SmootherBase : public VerbClass{

   public:

	SmootherBase( const XMLConfiguration* config , MultigridFGBase *grid ,
			      const ModelCollection* models ,
			      const FitobCalculator* fitobcalculator);

	virtual ~SmootherBase() {;}

	/** get the configuration object */
	const XMLConfiguration* config() const { return configuration_; }

	/** sets this micro time step , and also signals that diagonal and rhs must be calculated newly */
	virtual void setTimeStep( double microTimestep ) { microTimestep_ = microTimestep;}

	/** this method is used in multigrid to signal that it should not calculate the residuum */
	virtual void residumIsSet() = 0;

	/** this is the method which should be overwritten in case of many smoothers */
	virtual void smoothGrid( const DVector& globalCoords ) = 0;

	/** method to calculates the RHS */
	virtual void calculateRHS( DVector& globCoord) = 0;

	/** method to make an explicit midpoint step , the output should be in the unknowns vector (must be presized)
	 * @param globCoord [in/out] global coords , this will be changed
	 * @param unknowns [out] the updated unknowns , this must be presized!!!*/
	virtual void explicitStep( DVector& globCoord , DVector& unknowns , double underrelaxFactor ) = 0;

	/** method to calculate the residuum in L2 norm*/
	virtual double calcResiduum( const DVector& globalCoords ) = 0;

	/** pointer to the multi grid full grid */
	inline MultigridFGBase* grid() { return grid_; }

	inline const ModelCollection* models() { return models_; }

	/** */
	static boost::shared_ptr<SmootherBase> createSmoother(
			const XMLConfiguration* config , MultigridFGBase *grid ,
			const ModelCollection* models ,
			const FitobCalculator* fitobcalculator);

	/** */
	static SmootherBase* createSmootherPointer(
			const XMLConfiguration* config , MultigridFGBase *grid ,
			const ModelCollection* models ,
			const FitobCalculator* fitobcalculator);

   private:

	/** configuration object */
	const XMLConfiguration* configuration_;

	/** the grid on which we will calculate */
	MultigridFGBase *grid_;

	/** Pointer to the model collection which defines the diffusion */
	const ModelCollection* models_;

	/** the object which should contain all information*/
	const FitobCalculator* fitobcalculator_;

   protected:

	/** the timestep for execution */
	double microTimestep_;

	/** number of diffusion factors*/
	int nrFactors_;

	/** true if there is any correlation factor, FALSE if not*/
	bool hasCorrelations_;

	/** the dimension of the grid (the active dimensions of the domain)*/
	int gdim_;

};

}

#endif /* FITOBSMOOTHERBASE_HPP_ */

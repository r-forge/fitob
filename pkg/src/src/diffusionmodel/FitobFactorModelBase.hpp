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
 * FitobFactorModelBase.hpp
 *
 *  Created on: Apr 13, 2010
 *      Author: benk
 */

#ifndef FITOBFACTORMODELBASE_HPP_
#define FITOBFACTORMODELBASE_HPP_

#include  <boost/utility.hpp>

#include "src/utils/fitobdefs.hpp"
#include "src/variables/FitobVariable.hpp"

namespace fitob{

   class ModelCollection;

using namespace std;

    // forward declaration of the variable class
    //class Variable;

/** Base class for any kind the factor's model <br>.
 *  One such factor might be coupled with other factor (e.g. Heston)
 *  Such a model is responsible for the forward estimation of the <br>
 *
 *  The cross dependencies between variables should be solved on the
 *  subclass level (with the variable "averageVarValues").
 * */
class FactorModelBase : public boost::noncopyable , public VerbClass{

public:

	// so that it can call its private methods
	friend class ModelCollection;

	/** empty Ctor
	 * @param globalIndex the global index of this factor in the global variable table
	 * @param factorIndex the index of the factor in the configuration XML */
	FactorModelBase(const Variable* variable , int factorIndex);

	/** empty Dtor */
	virtual ~FactorModelBase() {;}

	/** returns the value of the convection coefficient (PDE)
	 * @param vars global variables in an array (coordinate)*/
	virtual double convectionCoef(const DVector &vars) const = 0;

	/** returns the value of the drift coefficient (SDE)
	 * @param vars global variables in an array (coordinate)*/
	virtual double driftCoef(const DVector &vars) const = 0;

	/** returns the value of the diffusion coefficient
	 * @param vars global variables in an array (coordinate)*/
	virtual double diffusionCoef(const DVector &vars) const = 0;

	/** returns the value of the diffusion coefficient
	 * @param initSize
	 * @param timeStep
	 * @param enlagementFactor
	 * @param averageVarValues, for coupled variables
	 * @param endSize [out] */
	virtual void forwardEstimation( double initSize,
			                        double timeStep,
			                        double enlargementFactor,
			                        const DVector &averageVarValues,
			                        double &endSize ) const = 0;

	/** formula used for Monte Carlo simulation
	 * @param nrPath number of Monte-Carlo paths
	 * @param actTime the actual time where the step starts
	 * @param timeStep , the time step
	 * @param nrGlobalVar , nr global variables
	 * @param InValues, one vector with global coordinates for each Monte Carlo path
	 * InOutValues[nrPath*nrGlobalVar + globalIndex]
	 * @param OutValues, the output vector the size of nrPath
	 * @param randN random numbers for each scenario size of nrPath */
	virtual void simulateForward( int nrPaths , double actTime, double timeStep , int nrGlobalVar ,
			DVector& InValues ,
			DVector& OutValues ,
			const DVector& randN ) const = 0;

	/** returns the discount factor for a given tenior P(t1,t2)
	 * @param vars global coordinates
	 * @param t1 start of the tenior
	 * @param t2 end of the tenior */
	virtual double discountFactor( const DVector &vars , double t1, double t2) const {
		FITOB_ERROR_EXIT("FactorModelBase::discountFactor , must be overwritten! , current method can not be used as interest rate model !!! ");
		return 1.0;
	}

	/** returns the global index of this factor*/
	const int inline getGlobalIndex() const {return variable_->getGlobalIndex();}

	/** returns the factor index, as it is defined in the configuration XML */
	const int inline getFactorIndex() const {return factorIndex_;}

	/** return the variable of this factor */
	const Variable* getVariable() const { return variable_;}

	/** returns the model collection*/
	inline static const ModelCollection* getModelCollect() { return modelCollection_;}

private:

	/** sets the variable for the interest rate */
	static void setModelCollect(ModelCollection* mc) { modelCollection_ = mc;}

	/** the global index of this factor as a variable*/
	const Variable* variable_;

	/** the index of this factor in the factor list (as it is defined in the configuration)*/
	const int factorIndex_;

	/** the variable (if there is any, if there is , then it is one) which is the risk-less interest rate*/
	static ModelCollection* modelCollection_;
};
}

#endif /* FITOBFACTORMODELBASE_HPP_ */

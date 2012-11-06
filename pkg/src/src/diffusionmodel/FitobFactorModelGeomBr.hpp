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
 * FitobFactorModelGeomBr.hpp
 *
 *  Created on: Apr 13, 2010
 *      Author: benk
 */

#ifndef FITOBFACTORMODELGEOMBR_HPP_
#define FITOBFACTORMODELGEOMBR_HPP_

#include "src/utils/fitobdefs.hpp"
#include "src/diffusionmodel/FitobFactorModelBase.hpp"
#include "src/utils/FitobXMLConfiguration.hpp"
#include "FitobModelCollection.hpp"

namespace fitob{

  using namespace std;

  /** Geometrical Browian model for the factor*/
  class FactorModelGeomBr : public FactorModelBase {
  public:

	FactorModelGeomBr(const boost::shared_ptr<XMLConfiguration> &ptr , const Variable* variable , int factorIndex);

	virtual ~FactorModelGeomBr() {;}

	/** returns the value of the convection coefficient
	 * @param vars global variables in an array (coordinate)*/
	virtual double convectionCoef(const DVector &vars) const{
		return (mu_<-1.0) ? (vars[getGlobalIndex()] * vars[getModelCollect()->getInterestRateGlobalIndex()]) : (vars[getGlobalIndex()] * mu_);
	}

	/** returns the value of the drift coefficient (SDE)
	 * @param vars global variables in an array*/
	virtual double driftCoef(const DVector &vars) const {
		return (drift_<-1.0) ? (vars[getGlobalIndex()] * vars[getModelCollect()->getInterestRateGlobalIndex()]) : (vars[getGlobalIndex()] * drift_);
	}

	/** returns the value of the diffusion coefficient
	 * @param vars global variables in an array*/
	virtual double diffusionCoef(const DVector &vars) const{
		return vars[getGlobalIndex()] * sigma_;
	}

	/** see base class for docu */
	virtual void forwardEstimation( double initSize,
			                        double timeStep,
			                        double enlargementFactor,
			                        const DVector &averageVarValues,
			                        double &endSize ) const {
		// estimate the end size after the given time step
		if (drift_<-1.0){ //vars[this->getModelCollect()->getInterestRateGlobalIndex()]
			   endSize = initSize * exp( averageVarValues[this->getModelCollect()->getInterestRateGlobalIndex()] * timeStep
							          + enlargementFactor * sigma_ * sqrt(timeStep));
		}else {
		   endSize = initSize * exp( (drift_ - r_) * timeStep
						          + enlargementFactor * sigma_ * sqrt(timeStep));
		}
	}

	/** see base class for docu */
	virtual void simulateForward( int nrPaths , double actTime, double timeStep , int nrGlobalVar ,
			DVector& InValues ,
			DVector& OutValues ,
			const DVector& randN ) const ;

	/** see base class for docu */
	virtual double discountFactor( const DVector &vars , double t1, double t2) const {
		// just return the discount factor with constant interest rate
		return ::exp( - vars[getGlobalIndex()] * (t2-t1) );
	}

  private:

	/** */
	double r_;

	/** */
	double mu_;

	/** */
	double drift_;

	/** */
	double sigma_;

  };
}

#endif /* FITOBFACTORMODELGEOMBR_HPP_ */

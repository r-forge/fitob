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
 * FitobFactorModelHeston.hpp
 *
 *  Created on: Nov 23, 2010
 *      Author: benk
 */

#ifndef FITOBFACTORMODELHESTON_HPP_
#define FITOBFACTORMODELHESTON_HPP_

#include "src/utils/fitobdefs.hpp"
#include "src/diffusionmodel/FitobFactorModelBase.hpp"
#include "src/utils/FitobXMLConfiguration.hpp"
#include "FitobModelCollection.hpp"

namespace fitob {

/** The Heston model , but only the part which is :<br>
 * dS(t) = mu*S*dt + sqrt(v)*S*dW <br>
 * where v is the underlying volatility process , this
 * hidden process must be previously declared as a "Heston-CIR" process
 * and this factor must be declared after the "Heston-CIR" process, so
 * automatically it will know that the hidden volatility process is the
 * previous */
class FactorModelHeston  : public FactorModelBase {
public:

	FactorModelHeston(const boost::shared_ptr<XMLConfiguration> &ptr ,
			const Variable* variable , int factorIndex ,
			const Variable* volaVariable , int volaFactorIndex );

	virtual ~FactorModelHeston() {;}

	/** returns the value of the convection coefficient
	 * @param vars global variables in an array*/
	virtual double convectionCoef(const DVector &vars) const{
		return (mu_<-1.0)?(vars[getGlobalIndex()] * vars[getModelCollect()->getInterestRateGlobalIndex()]) : (vars[getGlobalIndex()] * mu_);
	}

	/** returns the value of the drift coefficient (SDE)
	 * @param vars global variables in an array*/
	virtual double driftCoef(const DVector &vars) const {
		return (drift_<-1.0)?(vars[getGlobalIndex()] * vars[getModelCollect()->getInterestRateGlobalIndex()]) : (vars[getGlobalIndex()] * drift_);
	}

	/** returns the value of the diffusion coefficient
	 * @param vars global variables in an array*/
	virtual double diffusionCoef(const DVector &vars) const{
		//std::cout << " diffusionCoef HESTON =" << vars[getGlobalIndex()] * sqrt(vars[volatilityGlobalIndex_]) << std::endl;
		return vars[getGlobalIndex()] * sqrt(vars[volaVariable_->getGlobalIndex()]);
	}

	/** see base class for docu , (estimates the grid scaling on which it will be calculated)*/
	virtual void forwardEstimation( double initSize,
			                        double timeStep,
			                        double enlargementFactor,
			                        const DVector &averageVarValues,
			                               double &endSize ) const {
		// we take the average volatility and take as a constant volatility for the forward estimation
		if (drift_<-1.0){ //vars[this->getModelCollect()->getInterestRateGlobalIndex()]
			   endSize = initSize * exp( averageVarValues[this->getModelCollect()->getInterestRateGlobalIndex()] * timeStep
					               + enlargementFactor * sqrt( timeStep * averageVarValues[volaVariable_->getGlobalIndex()] ));
		} else {
		       endSize = initSize * exp( (drift_ - r_) * timeStep
						          + enlargementFactor * sqrt( timeStep * averageVarValues[volaVariable_->getGlobalIndex()] ));
		}
	}

	/** see base class for docu */
	virtual void simulateForward( int nrPaths ,double actTime, double timeStep , int nrGlobalVar ,
			DVector& InValues ,
			DVector& OutValues ,
			const DVector& randN ) const ;

private:

	/** the drift of the stock price (Heston process)*/
	double mu_;

	/** the average value of the risk free interest rate*/
	double r_;

	/** the coefficient for mesh estimation */
	double drift_;

	/** the global index of the underlying volatility */
	const Variable* volaVariable_;

	/** the factor index of the underlying volatility */
	int volatilityFactorIndex_;
};

}

#endif /* FITOBFACTORMODELHESTON_HPP_ */

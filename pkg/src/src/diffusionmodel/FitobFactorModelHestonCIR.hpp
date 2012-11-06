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
 * FitobFactorModelHestonCIR.hpp
 *
 *  Created on: Nov 23, 2010
 *      Author: benk
 */

#ifndef FITOBFACTORMODELHESTONCIR_HPP_
#define FITOBFACTORMODELHESTONCIR_HPP_

#include "src/utils/fitobdefs.hpp"
#include "src/diffusionmodel/FitobFactorModelBase.hpp"
#include "src/utils/FitobXMLConfiguration.hpp"

namespace fitob {

/** The CIR model which usually is the underlying volatility process of the Heston model <br>
 * dv = k*(theta-v(t))*dt + sigma*sqrt(v(t))*dW <br>
 * */
class FactorModelHestonCIR : public FactorModelBase {
public:

	/** */
	FactorModelHestonCIR(const boost::shared_ptr<XMLConfiguration> &ptr , const Variable* variable , int factorIndex);

	virtual ~FactorModelHestonCIR() {;}

	/** returns the value of the convection coefficient
	 * @param vars global variables in an array*/
	virtual double convectionCoef(const DVector &vars) const {
		return k_*( theta_ - vars[getGlobalIndex()]);
	}

	/** returns the value of the drift coefficient (SDE)
	 * @param vars global variables in an array*/
	virtual double driftCoef(const DVector &vars) const {
		return k_*( theta_ - vars[getGlobalIndex()]);
	}

	/** returns the value of the diffusion coefficient
	 * @param vars global variables in an array*/
	virtual double diffusionCoef(const DVector &vars) const{
		//std::cout << " diffusionCoef Heston-CIR =" << sigma_* sqrt(vars[getGlobalIndex()]) << std::endl;
		return sigma_ * sqrt(vars[getGlobalIndex()]);
	}

	/** see base class for docu */
	virtual void forwardEstimation( double initSize,
			                               double timeStep,
			                               double enlargementFactor,
			                         const DVector &averageVarValues,
			                               double &endSize ) const ;

	/** see base class for docu */
	virtual void simulateForward( int nrPaths , double actTime, double timeStep , int nrGlobalVar ,
			DVector& InValues ,
			DVector& OutValues ,
			const DVector& randN ) const ;

	/** see base class for docu */
	virtual double discountFactor( const DVector &vars , double t1, double t2) const {
		// See Interest Rate Models book, page 66, 3.25
		double h = ::sqrt(k_*k_ + 2*sigma_*sigma_);
		double A_tT = (2*h*exp( (k_+h)*(t2-t1)/2.0)) / (2*h+(k_+h)*(exp((t2-t1)*h)-1));
		A_tT = ::pow( A_tT , 2*k_*theta_/(sigma_*sigma_));
		double B_tT = (2*(exp(t2-t1)*h)-1) / (2*h+(k_+h)*exp((t2-t1)*h)-1);
		return A_tT * ::exp( - B_tT * vars[getGlobalIndex()]);
	}

private:

	/** the constant mean reversion speed */
	double k_;

	/** the constant mean reversion */
	double theta_;

	/** the diffusion coefficient */
	double sigma_;

	/** there are two different ways of make forward estimation: <br>
	 * a) using the mean and variance <br>
	 * b) using only the variance (avoids the drift, better for grid projection)<br>*/
	bool forwardEstimFlag_;
};

}

#endif /* FITOBFACTORMODELHESTONCIR_HPP_ */

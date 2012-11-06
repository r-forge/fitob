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
 * FitobFactorModelHullWhite.hpp
 *
 *  Created on: Nov 23, 2010
 *      Author: benk
 */

#ifndef FITOBFACTORMODELHULLWHITE_HPP_
#define FITOBFACTORMODELHULLWHITE_HPP_

#include "src/utils/fitobdefs.hpp"
#include "src/diffusionmodel/FitobFactorModelBase.hpp"
#include "src/utils/FitobXMLConfiguration.hpp"

namespace fitob {

	using namespace fitob;
	using namespace std;

/** Class which models the Hull White interest rate in the following way
 *  dr(t) = (theta(t) - a*r(t))dt + sigma * dW
 *  theta(t) = theta_mean*a + pow(sigma,2.0)*(1-exp(-2*a*(T-count*stepsize)))/(2*a);*/
class FactorModelHullWhite : public FactorModelBase {
public:

	FactorModelHullWhite(const boost::shared_ptr<XMLConfiguration> &ptr , const Variable* variable , int factorIndex);

	virtual ~FactorModelHullWhite() {;}

	/** returns the value of the convection coefficient
	 * @param vars global variables in an array*/
	virtual double convectionCoef(const DVector &vars) const{
		return theta(vars[0]) - a_*vars[getGlobalIndex()];
	}

	/** returns the value of the drift coefficient (SDE)
	 * @param vars global variables in an array*/
	virtual double driftCoef(const DVector &vars) const {
		return theta(vars[0]) - a_*vars[getGlobalIndex()];
	}

	/** returns the value of the diffusion coefficient
	 * @param vars global variables in an array*/
	virtual double diffusionCoef(const DVector &vars) const{
		return sigma_;
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
		// see the master thesis of Chao, for further informations (2.14) , P(t1,t2)
		// however that formula might need some correction ...
		double T = t2-t1;
		double ratio = (1.0-exp(-a_*T))/a_;
		double tmpt = exp(-a_*t2)-exp(-a_*t1);
		double complexterm = (1.0/(4.0*a_*a_*a_))*(sigma_*sigma_)*(tmpt*tmpt)*(exp(-2.0*a_*t1)-1);
		// just return the discount factor
		return ::exp( -theta_mean_*T + theta_mean_*ratio - complexterm - vars[getGlobalIndex()]*ratio);
	}


	/** return model specific parameters for the direct solver */
	double getA() const { return a_; }
	double getThetaMean() const { return theta_mean_; }
	double getSigma() const { return sigma_; }
	double getT() const { return T_; }

private:

	/** the mean reversion level depending on the tenior time*/
	inline double theta(double t) const {
		// t can not be greater than T_
		t = (t > T_)?T_:t;
		return theta_mean_ * a_  + sigma_*sigma_*(1-exp(-2*a_*(T_-t)))/(2*a_);
	}

	/** model coefficient*/
	double a_;

	/** mean reversion level of the interest rate*/
	double theta_mean_;

	/** volatility of the interest rate*/
	double sigma_;

	/** end time for mean reversion*/
	double T_;

	/** there are two different ways of make forward estimation: <br>
	 * a) using the mean and variance <br>
	 * b) using only the variance (avoids the drift, better for grid projection)<br>*/
	bool forwardEstimFlag_;
};

}

#endif /* FITOBFACTORMODELHULLWHITE_HPP_ */

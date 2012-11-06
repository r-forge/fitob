/*
 * FitobFactorModelHullWhite.cpp
 *
 *  Created on: Nov 23, 2010
 *      Author: benk
 */

#include "FitobFactorModelHullWhite.hpp"

using namespace fitob;

FactorModelHullWhite::FactorModelHullWhite(const boost::shared_ptr<XMLConfiguration> &ptr
	, const Variable* variable , int factorIndex) : FactorModelBase(variable,factorIndex) {

	// read in the model parameters
	//  dr(t) = (theta(t) - a*r(t))dt + sigma * dW
	//  theta(t) = theta_mean*a + pow(sigma,2.0)*(1-exp(-2*a*(T-count*stepsize)))/(2*a)

	a_ = ptr->getDoubleConfiguration("thetaconfigurations.thetaoperator.DIFFUSION_VARIABLES",factorIndex,"<xmlattr>.a" );

	theta_mean_ = ptr->getDoubleConfiguration("thetaconfigurations.thetaoperator.DIFFUSION_VARIABLES",factorIndex,"<xmlattr>.theta_mean" );

	sigma_ = ptr->getDoubleConfiguration("thetaconfigurations.thetaoperator.DIFFUSION_VARIABLES",factorIndex,"<xmlattr>.sigma" );

	T_ = ptr->getDoubleConfiguration("thetaconfigurations.thetaoperator.DIFFUSION_VARIABLES",factorIndex,"<xmlattr>.T" );

	forwardEstimFlag_ = false;
	string tmp = ptr->getStringConfiguration("thetaconfigurations.thetaoperator.DIFFUSION_VARIABLES",factorIndex,"<xmlattr>.MeanFrwrdEstim");
	if (tmp == "true"){
		forwardEstimFlag_ = true;
	}

	FITOB_OUT_LEVEL3(verb() , " FactorModelHullWhite created with T=" << T_ <<" , a="<<a_<<" , theta_mean="<<theta_mean_<<" , sigma="<<sigma_);
	FITOB_OUT_LEVEL3(verb() , " FactorModelHullWhite globalVarIndex=" << variable->getGlobalIndex() <<" , factorIndex="<<factorIndex);

}


void FactorModelHullWhite::forwardEstimation( double initSize,
		                               double timeStep,
		                               double enlargementFactor,
		                         const DVector &averageVarValues,
		                               double &endSize ) const {
	// estimate the end size after the given time step
	// mean_r = r0*exp(-a.*t) + (tet./a).*(1-exp(-a.*t));
	// var_r =((sigma*sigma)/(2.*a))*(1-exp(-2.*a.*t));
	double mean_r = initSize*exp(-a_*timeStep) + ( theta(timeStep) / a_)*( 1.0 - exp(-a_*timeStep));
	double var_r = ((sigma_*sigma_)/(2*a_))*(1-exp(-2*a_*timeStep));
	double std_r = sqrt(var_r);

	//FITOB_OUT_LEVEL3(4, " mean_r = " << mean_r << " , var_r = " << var_r << " , std_r = " << std_r);

	if (forwardEstimFlag_)
	{
		// first version: this allows negative interest rates
		endSize = mean_r + std_r * enlargementFactor;
		// this is nothing else than  = mean(t) + std(t) * enlargementFactor;
		// todo: we must make sure that that the larger domain will always contain the smaller domain
		// practically that the variance is the stronger term
	}
	else // this is the default version
	{
		// second version: this also allows negative interest rates
		endSize = initSize + std_r * enlargementFactor;
		// this is nothing else than  = init_Value + std(t) * enlargementFactor;
	}
}


void FactorModelHullWhite::simulateForward( int nrPaths , double actTime, double timeStep , int nrGlobalVar ,
		DVector& InValues ,
		DVector& OutValues ,
		const DVector& randN ) const {

	//double mean_r = 0.0;
	//double var_r = 0.0;
	//double std_r = 0.0;
	double sqrtStep = sqrt(timeStep);
#if ( defined(FITOB_OPENMP) )
#pragma omp parallel firstprivate(sqrtStep)
{
#pragma omp for schedule(static)
#endif
	for (int p = 0 ; p < nrPaths ; p++ ){
		/**  dr(t) = (theta(t) - a*r(t))dt + sigma * dW
          *  theta(t) = theta_mean*a + pow(sigma,2.0)*(1-exp(-2*a*(T-count*stepsize)))/(2*a); */
		// estimate the end size after the given time step
		// mean_r = r0*exp(-a.*t) + (tet./a).*(1-exp(-a.*t));
		// var_r =((sigma*sigma)/(2.*a))*(1-exp(-2.*a.*t));
		/*
		mean_r = InValues[nrGlobalVar*p + this->getGlobalIndex() ]* exp(-a_*timeStep) + ( theta(actTime-timeStep) / a_)*( 1.0 - exp(-a_*timeStep));
		var_r = ((sigma_*sigma_)/(2*a_))*(1-exp(-2*a_*timeStep));
		std_r = sqrt(var_r);
		OutValues[ p ] = mean_r + std_r * randN[p]; */
		OutValues[ p ] = InValues[nrGlobalVar*p + this->getGlobalIndex()]
				+ (theta(actTime-timeStep) -a_*InValues[nrGlobalVar*p + this->getGlobalIndex()])*timeStep
		        + sigma_ * sqrtStep * randN[p];
	}
#if ( defined(FITOB_OPENMP) )
}
#endif
}

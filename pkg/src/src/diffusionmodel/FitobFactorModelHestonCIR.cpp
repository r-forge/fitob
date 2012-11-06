/*
 * FitobFactorModelHestonCIR.cpp
 *
 *  Created on: Nov 23, 2010
 *      Author: benk
 */

#include "FitobFactorModelHestonCIR.hpp"

using namespace fitob;

FactorModelHestonCIR::FactorModelHestonCIR(const boost::shared_ptr<XMLConfiguration> &ptr ,
		const Variable* variable , int factorIndex) : FactorModelBase(variable,factorIndex)  {

	//dv = k*(theta-v(t))*dt + sigma*sqrt(v(t))*dW
	setVerb(4);

	k_ = ptr->getDoubleConfiguration("thetaconfigurations.thetaoperator.DIFFUSION_VARIABLES",factorIndex,"<xmlattr>.k" );

	theta_ = ptr->getDoubleConfiguration("thetaconfigurations.thetaoperator.DIFFUSION_VARIABLES",factorIndex,"<xmlattr>.theta" );

	sigma_ = ptr->getDoubleConfiguration("thetaconfigurations.thetaoperator.DIFFUSION_VARIABLES",factorIndex,"<xmlattr>.sigma" );

	forwardEstimFlag_ = false;
	string tmp = ptr->getStringConfiguration("thetaconfigurations.thetaoperator.DIFFUSION_VARIABLES",factorIndex,"<xmlattr>.MeanFrwrdEstim");
	if (tmp == "true"){
		forwardEstimFlag_ = true;
	}

	FITOB_OUT_LEVEL3(verb() , " FactorModelHestonCIR created with k=" << k_ <<" , theta="<<theta_<<" , sigma="<<sigma_);
	FITOB_OUT_LEVEL3(verb() , " FactorModelHestonCIR globalVarIndex=" << variable->getGlobalIndex() <<" , factorIndex="<<factorIndex);
}


void FactorModelHestonCIR::forwardEstimation( double initSize,
		                               double timeStep,
		                               double enlargementFactor,
		                         const DVector &averageVarValues,
		                               double &endSize ) const {
	// see docu for more information (e.g. wiki)
	// mean_v = v0*exp(-k.*t) + theta*(1-exp(-k.*t));
	// var_v = v0 * (sigma^2/k)*(exp(-k*t)-exp(-2*k*t))+ (theta*sigma^2/(2*k))*(1-exp(-k*t))^2

	double mean_v = initSize*exp(-k_*timeStep) + theta_*( 1 - exp(-k_*timeStep));
	double var_v = initSize*(sigma_*sigma_/ k_)*( exp(-k_*timeStep) - exp(-2*k_*timeStep))
			       + (theta_*sigma_*sigma_/(2*k_))*pow( 1 - exp(-k_*timeStep),2);
	double std_v = sqrt(var_v);

	//FITOB_OUT_LEVEL3( verb(), " mean_v = " << mean_v << " , var_v = " << var_v << " , std_v = " << std_v);

	if (forwardEstimFlag_)
	{
		// first version:
		if (enlargementFactor >= 0){
			endSize = mean_v + std_v * enlargementFactor;
		} else {
			// here we use the exp formula in order to not allow negative volatilities
			endSize = mean_v * exp( std_v * enlargementFactor );
		}
		// this is nothing else than  = mean(t) + std(t) * enlargementFactor;
		// todo: we must make sure that that the larger domain will always contain the smaller domain
		// practically that the variance is the stronger term, then this is satisfied
	}
	// this is the default procedure
	else
	{
		// second version:
		if (enlargementFactor >= 0){
			endSize = initSize + std_v * enlargementFactor;
		} else {
			// here we use the exp formula in order to not allow negative volatilities
			endSize = initSize * exp( std_v * enlargementFactor );
		}
		// this is nothing else than  = init_Value + std(t) * enlargementFactor;
	}
}

void FactorModelHestonCIR::simulateForward( int nrPaths , double actTime, double timeStep , int nrGlobalVar ,
		DVector& InValues ,
		DVector& OutValues ,
		const DVector& randN ) const {

	//double mean_v = 0.0;
	//double var_v = 0.0;
	//double std_v = 0.0;
	//* dv = k*(theta-v(t))*dt + sigma*sqrt(v(t))*dW <br>
	// see docu for more information (e.g. wiki)
	// mean_v = v0*exp(-k.*t) + theta*(1-exp(-k.*t));
	// var_v = v0 * (sigma^2/k)*(exp(-k*t)-exp(-2*k*t))+ (theta*sigma^2/(2*k))*(1-exp(-k*t))^2
	double sqrtTime = sqrt(timeStep) , sigmaSqrt = 0.0;
#if ( defined(FITOB_OPENMP) )
#pragma omp parallel firstprivate(sqrtTime,sigmaSqrt)
{
#pragma omp for schedule(static)
#endif
	for (int p = 0 ; p < nrPaths ; p++ ){
		/*
		mean_v = InValues[nrGlobalVar*p + this->getGlobalIndex() ] * exp(-k_*timeStep) + theta_*( 1 - exp(-k_*timeStep));
		var_v = InValues[nrGlobalVar*p + this->getGlobalIndex() ] *
				(sigma_*sigma_/ k_)*( exp(-k_*timeStep) - exp(-2*k_*timeStep))
				    + (theta_*sigma_*sigma_/(2*k_))*pow( 1 - exp(-k_*timeStep),2);
		std_v = sqrt(var_v);
		OutValues[ p ] = mean_v + std_v*randN[p];
        */
		// just do the simple dicretization
		sigmaSqrt = sqrt( InValues[nrGlobalVar*p + this->getGlobalIndex()] );
		OutValues[ p ] = InValues[nrGlobalVar*p + this->getGlobalIndex()] + k_*(theta_ - InValues[nrGlobalVar*p + this->getGlobalIndex()] ) * timeStep +
				sigma_ * sigmaSqrt * sqrtTime * randN[p];
		// if the value turns negative , then turn back to positive
		if ( OutValues[ p ] < 0.0 ) { OutValues[ p ] = - OutValues[ p ]; }
		//FITOB_OUT_LEVEL3( 5 , "HEST-CIR old=" << InValues[nrGlobalVar*p + this->getGlobalIndex()] <<" , new="<<OutValues[p]);
	}
#if ( defined(FITOB_OPENMP) )
}
#endif

}

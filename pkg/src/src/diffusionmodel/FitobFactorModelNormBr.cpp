/*
 * FitobFactorModelNormBr.cpp
 *
 *  Created on: Apr 13, 2010
 *      Author: benk
 */

#include "FitobFactorModelNormBr.hpp"

using namespace fitob;
using namespace std;

FactorModelNormBr::FactorModelNormBr(const boost::shared_ptr<XMLConfiguration> &ptr , const Variable* variable , int factorIndex)
: FactorModelBase( variable ,factorIndex){
    // read in from the configuration object the necessary informations
	setVerb(4);

	// even though this might be coupled on one variable we might read in the average value for forward estimation
	r_ = ptr->getDoubleConfiguration("thetaconfigurations.thetaoperator.RISK_FREE_RATE.<xmlattr>.value" );

	mu_ = ptr->getDoubleConfiguration("thetaconfigurations.thetaoperator.DIFFUSION_VARIABLES",factorIndex,"<xmlattr>.convec" );

	drift_ = ptr->getDoubleConfiguration("thetaconfigurations.thetaoperator.DIFFUSION_VARIABLES",factorIndex,"<xmlattr>.drift" );

	sigma_ = ptr->getDoubleConfiguration("thetaconfigurations.thetaoperator.DIFFUSION_VARIABLES",factorIndex,"<xmlattr>.sigma" );

	FITOB_OUT_LEVEL3(verb() , " FactorModelNormBr created with r="<<r_<<" , mu="<<mu_<<" , drift="<<drift_<<" , sigma="<<sigma_);
}

void FactorModelNormBr::simulateForward( int nrPaths , double actTime, double timeStep , int nrGlobalVar ,
		DVector& InValues ,
		DVector& OutValues ,
		const DVector& randN ) const {

	if (drift_ < -1.0) {
		  double sqrtTime = sqrt(timeStep);
#if ( defined(FITOB_OPENMP) )
#pragma omp parallel firstprivate(sqrtTime)
{
#pragma omp for schedule(static)
#endif
		  for (int p = 0 ; p < nrPaths ; p++ ){
			  OutValues[ p ] = InValues[nrGlobalVar*p + this->getGlobalIndex() ]
			               + InValues[nrGlobalVar*p + getModelCollect()->getInterestRateGlobalIndex()] * timeStep + sigma_ * sqrtTime * randN[p] ;
		  }
#if ( defined(FITOB_OPENMP) )
}
#endif
	} else {
		double sqrtTime = sqrt(timeStep);
#if ( defined(FITOB_OPENMP) )
#pragma omp parallel firstprivate(sqrtTime)
{
#pragma omp for schedule(static)
#endif
		for (int p = 0 ; p < nrPaths ; p++ ){
			OutValues[ p ] = InValues[nrGlobalVar*p + this->getGlobalIndex() ]
		               + drift_ * timeStep + sigma_ * sqrtTime * randN[p] ;
		}
#if ( defined(FITOB_OPENMP) )
}
#endif
	}
}

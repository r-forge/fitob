/*
 * FitobFactorModelGeomBr.cpp
 *
 *  Created on: Apr 13, 2010
 *      Author: benk
 */

#include "FitobFactorModelGeomBr.hpp"

using namespace fitob;
using namespace std;

FactorModelGeomBr::FactorModelGeomBr(const boost::shared_ptr<XMLConfiguration> &ptr , const Variable* variable , int factorIndex):
  FactorModelBase(variable,factorIndex){

	setVerb(4);

    // read in from the configuration object the necessary informations

	// even though this might be coupled on one variable we might read in the average value for forward estimation
	r_ = ptr->getDoubleConfiguration("thetaconfigurations.thetaoperator.RISK_FREE_RATE.<xmlattr>.value" );

	mu_ = ptr->getDoubleConfiguration("thetaconfigurations.thetaoperator.DIFFUSION_VARIABLES",factorIndex,"<xmlattr>.convec" );

	drift_ = ptr->getDoubleConfiguration("thetaconfigurations.thetaoperator.DIFFUSION_VARIABLES",factorIndex,"<xmlattr>.drift" );

	sigma_ = ptr->getDoubleConfiguration("thetaconfigurations.thetaoperator.DIFFUSION_VARIABLES",factorIndex,"<xmlattr>.sigma" );

	FITOB_OUT_LEVEL3(verb() , " FactorModelGeomBr created with variable=" << variable->getGlobalIndex()
			<< " r="<<r_<<" mu="<<mu_<<" drift="<<drift_<<" sigma="<<sigma_);
}

void FactorModelGeomBr::simulateForward( int nrPaths , double actTime, double timeStep , int nrGlobalVar ,
		DVector& InValues ,
		DVector& OutValues ,
		const DVector& randN ) const {

	double sigma2 = sigma_* sigma_*0.5;
	double sqrtTime = sqrt(timeStep);

	if (drift_ < -1.0) {
#if ( defined(FITOB_OPENMP)  )
#pragma omp parallel firstprivate(sigma2,sqrtTime)
{
#pragma omp for schedule(static)
#endif
		  for (int p = 0 ; p < nrPaths ; p++ ){
			  OutValues[ p ] = InValues[ nrGlobalVar*p + this->getGlobalIndex() ]
			    * exp( (InValues[nrGlobalVar*p + getModelCollect()->getInterestRateGlobalIndex()] - sigma2) * timeStep + sigma_ * sqrtTime * randN[p] );
		  }
#if ( defined(FITOB_OPENMP)  )
}
#endif
	} else {
#if ( defined(FITOB_OPENMP)  )
#pragma omp parallel firstprivate(sigma2,sqrtTime)
{
#pragma omp for schedule(static)
#endif
	      for (int p = 0 ; p < nrPaths ; p++ ){
		      OutValues[ p ] = InValues[ nrGlobalVar*p + this->getGlobalIndex() ]
		              * exp( (drift_ - sigma2) * timeStep + sigma_ * sqrtTime * randN[p] );
	      }
#if ( defined(FITOB_OPENMP) )
}
#endif
	}
}

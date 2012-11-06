/*
 * FitobFactorModelHeston.cpp
 *
 *  Created on: Nov 23, 2010
 *      Author: benk
 */

#include "FitobFactorModelHeston.hpp"

using namespace fitob;

FactorModelHeston::FactorModelHeston(const boost::shared_ptr<XMLConfiguration> &ptr
 , const Variable* variable , int factorIndex
 , const Variable* volaVariable , int volaFactorIndex ) : FactorModelBase(variable,factorIndex)  {

	// dS(t) = mu*S*dt + sqrt(v)*S*dW
	setVerb(4);
	// ------- !!! the underlying volatility factor should be the one before this factor (in the XML) !!! ------------
	volaVariable_ = volaVariable;
	volatilityFactorIndex_ = volaFactorIndex;

	r_ = ptr->getDoubleConfiguration("thetaconfigurations.thetaoperator.RISK_FREE_RATE.<xmlattr>.value" );

	mu_ = ptr->getDoubleConfiguration("thetaconfigurations.thetaoperator.DIFFUSION_VARIABLES",factorIndex,"<xmlattr>.convec" );

	drift_ = ptr->getDoubleConfiguration("thetaconfigurations.thetaoperator.DIFFUSION_VARIABLES",factorIndex,"<xmlattr>.drift" );

	FITOB_OUT_LEVEL3(verb() , " FactorModelHeston created with , r="<<r_<<" , mu="<<mu_<<" , drift="<<drift_ <<
			" , volaGlobalIndex = " << volaVariable_->getGlobalIndex() << " , volaFactorIndex:" << volatilityFactorIndex_ );
	FITOB_OUT_LEVEL3(verb() , " FactorModelHeston globalVarIndex=" << variable->getGlobalIndex() <<" , factorIndex="<<factorIndex);
}

/** see base class for docu */
void FactorModelHeston::simulateForward( int nrPaths , double actTime,  double timeStep , int nrGlobalVar ,
		DVector& InValues ,
		DVector& OutValues ,
		const DVector& randN ) const {

	double sigmaSqrt = 0.0;
	double sqrtTime = sqrt(timeStep);

	if (drift_ < -1.0){
#if ( defined(FITOB_OPENMP) )
#pragma omp parallel firstprivate(sigmaSqrt,sqrtTime)
{
#pragma omp for schedule(static)
#endif
		  for (int p = 0 ; p < nrPaths ; p++ ){
		      sigmaSqrt = sqrt( InValues[nrGlobalVar*p + volaVariable_->getGlobalIndex()] );
			  OutValues[ p ] = InValues[ nrGlobalVar*p + this->getGlobalIndex() ]
			    * exp( (InValues[nrGlobalVar*p + getModelCollect()->getInterestRateGlobalIndex()] -
			    		InValues[nrGlobalVar*p + volaVariable_->getGlobalIndex()]/2) * timeStep + sqrtTime * sigmaSqrt * randN[p]);
		  }
#if ( defined(FITOB_OPENMP) )
}
#endif
	} else {
#if ( defined(FITOB_OPENMP) )
#pragma omp parallel firstprivate(sigmaSqrt,sqrtTime)
{
#pragma omp for schedule(static)
#endif
	      for (int p = 0 ; p < nrPaths ; p++ ){
		      sigmaSqrt = sqrt( InValues[nrGlobalVar*p + volaVariable_->getGlobalIndex()] );
		      OutValues[ p ] = InValues[nrGlobalVar*p + this->getGlobalIndex() ]
		              * exp( (drift_ - InValues[nrGlobalVar*p + volaVariable_->getGlobalIndex()]/2) * timeStep + sqrtTime * sigmaSqrt * randN[p] );
		      //OutValues[ p ] = InValues[nrGlobalVar*p + this->getGlobalIndex() ] + drift_ * InValues[nrGlobalVar*p + this->getGlobalIndex() ] * timeStep
		      //		          +  sqrtTime * sigmaSqrt * InValues[nrGlobalVar*p + this->getGlobalIndex() ] * randN[p];
		      //FITOB_OUT_LEVEL3( 5 , "HEST old=" << InValues[nrGlobalVar*p + this->getGlobalIndex()] <<" , new="<<OutValues[p]);
	      }
#if ( defined(FITOB_OPENMP) )
}
#endif
	}
}

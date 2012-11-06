/*
 * FitobFactorModelGeneral.cpp
 *
 *  Created on: Aug 7, 2012
 *      Author: benk
 */

#include "FitobFactorModelGeneral.hpp"

using namespace fitob;

FactorModelGeneral::FactorModelGeneral(	const Variable* variable , int factorIndex, ExpressionBasis* drift,
		 ExpressionBasis* convec, ExpressionBasis* sigma, ExpressionBasis* discount) :
		 FactorModelBase(variable , factorIndex) , drift_(drift) , convec_(convec) , sigma_(sigma) , discount_(discount){
}

FactorModelGeneral::FactorModelGeneral(const FactorModelGeneral* copyFact , int factorIndex)
: FactorModelBase( copyFact->getVariable() , factorIndex) ,
		drift_(copyFact->drift_) , convec_(copyFact->convec_) , sigma_(copyFact->sigma_) , discount_(copyFact->discount_){
}

FactorModelGeneral::~FactorModelGeneral() {
	// todo: delete the pointers ... it is probably better to delete them where they were created ?
}

/** see parent class */
double FactorModelGeneral::convectionCoef(const DVector &vars) const {
	return convec_->eval(vars);
}

/** see parent class */
double FactorModelGeneral::driftCoef(const DVector &vars) const {
	return drift_->eval(vars);
}

/** see parent class */
double FactorModelGeneral::diffusionCoef(const DVector &vars) const {
	return sigma_->eval(vars);
}

/** see parent class */
void FactorModelGeneral::forwardEstimation( double initSize,
		                        double timeStep,
		                        double enlargementFactor,
		                        const DVector &averageVarValues,
		                        double &endSize ) const{
	//averageVarValues[this->getGlobalIndex()] = initSize;
	endSize = initSize;
	// todo: make small time steps ?
	const int nrTimeSt = 10; //::ceil(timeStep/0.1);
	for (int st = 0 ; st < nrTimeSt ; st++ ){
		endSize = endSize +
				drift_->eval(averageVarValues)*(timeStep/(double)nrTimeSt) +
				sigma_->eval(averageVarValues) * enlargementFactor * ::sqrt(timeStep/(double)nrTimeSt);
	}
}

/** see parent class */
void FactorModelGeneral::simulateForward( int nrPaths , double actTime, double timeStep , int nrGlobalVar ,
		DVector& InValues ,
		DVector& OutValues ,
		const DVector& randN ) const{

	double sqrtTime = sqrt(timeStep);
	int tmp = 0;
#if ( defined(FITOB_OPENMP)  )
#pragma omp parallel firstprivate(tmp,sqrtTime)
{
#endif
		DVector tmpGlobalValues(nrGlobalVar,0.0);
#if ( defined(FITOB_OPENMP)  )
#pragma omp for schedule(static)
#endif
	      for (int p = 0 ; p < nrPaths ; p++ ){
	    	  // copy the global variables into one vector
	    	  for (tmp = 0 ; tmp < nrGlobalVar ; tmp++) { tmpGlobalValues[tmp] = InValues[ nrGlobalVar*p + tmp];}
	    	  // just evaluate the linearly discretized SDE
		      OutValues[ p ] = InValues[ nrGlobalVar*p + this->getGlobalIndex() ]
		              + drift_->eval(tmpGlobalValues) * timeStep
		              + sigma_->eval(tmpGlobalValues) * sqrtTime * randN[p];
	      }
#if ( defined(FITOB_OPENMP) )
}
#endif

}

/** see parent class */
double FactorModelGeneral::discountFactor( const DVector &vars , double t1, double t2) const{
	// todo: implement this, with the time ... ?
	// todo: this is not yet clear how to make this
	return 1.0;
}

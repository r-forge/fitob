/*
 * FitobMCMachine.cpp
 *
 *  Created on: Mar 14, 2011
 *      Author: benk
 */

#include "FitobMCMachine.hpp"

using namespace fitob;
using namespace std;

MCMachine::MCMachine(const XMLConfiguration* xmlConfig ,
		boost::shared_ptr<ModelCollection>& modelcollection) : xmlConfig_(xmlConfig) ,
		modelcollection_(modelcollection) , randomNrGenerator_(xmlConfig) {

	setVerb(3);

	// initially no MC steps
	nrMCSteps_ = 0;

	// this is read in from the XML file
	nrMCScenarios_ = xmlConfig->getIntConfiguration("thetaconfigurations.montecarlo.PATH_NR.<xmlattr>.value");
	initialDiffusionSTD_ = xmlConfig->getDoubleConfiguration("thetaconfigurations.montecarlo.INITIAL_DIFFUSION.<xmlattr>.sigma");
	if (nrMCScenarios_ <= 0 ) nrMCScenarios_ = 10;

	// read in the maximal time step for the MC simulation
	maxDT_ = xmlConfig->getDoubleConfiguration("thetaconfigurations.montecarlo.MC_TIME_STEP.<xmlattr>.value");
	if (maxDT_ <= 0 ) maxDT_ = 1.0;

	// get the evaluation modes
	xmlConfig->getIntVectorConfiguration("thetaconfigurations.montecarlo.EVALUATION.<xmlattr>.mode" , ',' , evaluationModes_ );

	// the correlation matrix
	const DVector& corrM = modelcollection->getCorrelationMatrix();

	// number of factors
	int nrFact = modelcollection->nrFactors();
	cholevskyMatrix_.resize( nrFact*nrFact , 0.0 );

	// =========== here comes the Cholesky factorization ===============
/** Matlab Code
n = length( M );
L = zeros( n, n );
for i=1:n
    L(i, i) = sqrt( M(i, i) - L(i, :)*L(i, :)' );
    for j=(i + 1):n
        L(j, i) = ( M(j, i) - L(i, :)*L(j, :)' )/L(i, i);
    end
end */

	int i,j,k;
	double sum;
	FITOB_OUT_LEVEL2(verb()," Cholevsky decomposition nrFact:"<<nrFact);
	for (i = 0 ; i < nrFact ; i++ ){
		sum = 0.0;
		for ( k = 0; k < nrFact ; k++) sum += cholevskyMatrix_[ i*nrFact + k]*cholevskyMatrix_[ i*nrFact + k];
		cholevskyMatrix_[ i*nrFact + i ] = sqrt(  corrM[ i*nrFact + i] - sum );
		FITOB_OUT_LEVEL2(verb()," cholevskyMatrix_["<<(i*nrFact + i)<<"]= " << cholevskyMatrix_[ i*nrFact + i ] );
		for (j = i+1 ; j < nrFact ; j++){
			sum = 0.0;
			for ( k = 0; k < nrFact ; k++) sum += cholevskyMatrix_[ i*nrFact + k]*cholevskyMatrix_[ j*nrFact + k];
			cholevskyMatrix_[ j*nrFact + i ] = (corrM[ i*nrFact + j] - sum) / cholevskyMatrix_[ i*nrFact + i ];
			// the other part is always Zero -> cholevskyMatrix_[ i*nrFact + j ] = 0.0;
			FITOB_OUT_LEVEL2(verb()," cholevskyMatrix_["<<(j*nrFact + i)<<"]= " << cholevskyMatrix_[ j*nrFact + i ] );
		}
	}

}

void MCMachine::addMCStep(const Domain& domain , double Acttime , int evalContextIndex ){

	MCSimulations_.push_back(new MCStep( domain ,  nrMCScenarios_ , Acttime , evalContextIndex ));

	nrMCSteps_++;

	if ( (nrMCSteps_ >= 1) && (initialDiffusionSTD_ > 0.0) ){
		// in the case of the first MCStep we might make here an initial distribution of the scenarios
		// so that not all the scenarios will start from the same spot
		// here we add to all risk factors a normally distributed
		const ModelCollection* modelcollection = modelcollection_.get();
		int nrFactor = modelcollection->nrFactors();
		int s , tmp, factorGloablI , variableIndexPrev;
		double tmp_val = 0.0;

		DVector unCorrelatedRandNum( nrFactor*nrMCScenarios_ );
		randomNrGenerator_.getRandomN( nrFactor*nrMCScenarios_ , unCorrelatedRandNum );

		for (s = 0 ; s < nrMCScenarios_ ; s++){
			for (tmp = 0 ; tmp < nrFactor ; tmp++ ) {
				factorGloablI = modelcollection->getModel(tmp).getGlobalIndex();
				variableIndexPrev = MCSimulations_[nrMCSteps_-1].getVariableIndex( factorGloablI );
				tmp_val = MCSimulations_[nrMCSteps_-1].getSimulationValue( s , variableIndexPrev ) +
						 (  (domain.getGradedAxisMax(factorGloablI) - domain.getGradedAxisMin(factorGloablI))
						    *initialDiffusionSTD_*unCorrelatedRandNum[s*nrFactor + tmp]
						 );
				MCSimulations_[nrMCSteps_-1].setSimulationValue( s , variableIndexPrev , tmp_val  );
			}
		}
	}

	// if this is the first step then just return
	if (nrMCSteps_ < 2) return;

	// take the time stamp from the actual and previous, only if there is a difference
	// make a simulation of the factors

	double dt = Acttime - MCSimulations_[nrMCSteps_-2].getAverage()[0];
	FITOB_OUT_LEVEL2(verb(),"MCMachine::addMCStep dt = " << dt << " , Acttime:" << Acttime << " , prev: "
			<< MCSimulations_[nrMCSteps_-2].getAverage()[0]);

	// copy the values
	MCSimulations_[nrMCSteps_-1].setValues( &(MCSimulations_[nrMCSteps_-2]) );

	if (dt > 1e-10)
	{
		// make the forward simulation of the risk factors
		makeForwardMCStep( Acttime , dt ,  MCSimulations_[nrMCSteps_-2] ,  MCSimulations_[nrMCSteps_-1] );
	}

	// print the actual step
	FITOB_OUT_LEVEL4(verb(),"FitobCalculator::addMCStep NewStep:" << MCSimulations_[nrMCSteps_-1].toString() );
}



void  MCMachine::makeForwardMCStep( double actTime, double dt , MCStep& prevMCStep , MCStep& actMCStep ){

	const ModelCollection* modelcollection = modelcollection_.get();
	int buffSize = 1024;
	int tmp , tmp1 , nrGlobalVar = actMCStep.getDom().nrGlobalVariables();
	// nr of time steps and the constant time step
	int nrSDE_Step = ::ceil(dt/maxDT_);
	double actDt = dt / (double) nrSDE_Step;

	int nrFactor = modelcollection->nrFactors();
	DVector correlatedRandNum( nrSDE_Step * nrFactor * nrMCScenarios_ );
	DVector startValue( nrGlobalVar * buffSize );
	std::vector< DVector > outValue( nrFactor );
	std::vector< DVector > randNr( nrFactor );

	FITOB_OUT_LEVEL2(verb(),"MCMachine::makeForwardMCStep dt = " << dt << " , actTime = " << actTime << " , actDt = " << actDt);

	//  HERE WE MAKE A FORWARD STEP OF THE RISK FACTORS

	// generate uncorrelated random numbers uniformly distributed
	randomNrGenerator_.getRandomN( nrSDE_Step * nrFactor * nrMCScenarios_ , correlatedRandNum );


#if ( defined(FITOB_OPENMP) )
#pragma omp parallel private(tmp,tmp1)
{
#endif
	DVector tmpVars(nrFactor);
	// todo: use BLAS routine later for this purpose, to apply the correlation Cholesky
	double sum = 0.0;
#if ( defined(FITOB_OPENMP) )
#pragma omp for schedule(static)
#endif
	for (int s = 0 ; s < nrSDE_Step * nrMCScenarios_ ; s++)
	{
		// copy to a temporary vector
		for (tmp = 0 ; tmp < nrFactor ; tmp++ ) { tmpVars[tmp] = correlatedRandNum[s*nrFactor + tmp]; }

		// matrix vector multiplication
		for (tmp = 0 ; tmp < nrFactor ; tmp++ ){
			sum = 0.0;
			for (tmp1 = 0 ; tmp1 < nrFactor ; tmp1++ )
				{ sum = sum + cholevskyMatrix_[tmp*nrFactor + tmp1 ] * tmpVars[tmp1]; }
			correlatedRandNum[s*nrFactor + tmp] = sum;
			//FITOB_OUT_LEVEL2(verb()," correlatedRandNum["<<(s*nrFactor + tmp)<<"]= " << sum );
		}
	}
#if ( defined(FITOB_OPENMP) )
}
#endif

	int nrs = 0 , factorGloablI = 0 , variableIndexPrev = 0 , variableIndexAct = 0;

	// iterate with a buffer over the nr of scenarios
	while (nrs < nrMCScenarios_)
	{
	    // for each factor get the initial value
		for (int fa = 0 ; fa < nrFactor ; fa++)
		{
			randNr[fa].resize(buffSize);
			outValue[fa].resize(buffSize);
			// the first local variables are always the diffusion variables
			// this is the convention
			factorGloablI = modelcollection->getModel(fa).getGlobalIndex();
			variableIndexPrev = prevMCStep.getVariableIndex( factorGloablI );
			// copy in the buffer the starting values and the random numbers
			for (tmp = 0; (tmp < buffSize) && (nrs < nrMCScenarios_) ; tmp++ , nrs++ )
			{
				startValue[ nrGlobalVar*tmp + factorGloablI] = prevMCStep.getSimulationValue( nrs , variableIndexPrev );
				outValue[fa][tmp] = startValue[ nrGlobalVar*tmp + factorGloablI];
				randNr[fa][tmp] = correlatedRandNum[ nrFactor*nrs + fa ];
			}
			nrs = nrs - tmp;
		}

		// first copy the initial values to the vector
		for (int nrTimeStep = 0 ; nrTimeStep < nrSDE_Step ; nrTimeStep++ )
		{
			// copy the values
			for (int fa = 0 ; fa < nrFactor ; fa++)
			{
				factorGloablI = modelcollection->getModel(fa).getGlobalIndex();
				variableIndexPrev = prevMCStep.getVariableIndex( factorGloablI );
#if ( defined(FITOB_OPENMP) )
#pragma omp parallel firstprivate(tmp,tmp1)
{
#pragma omp for schedule(static)
#endif
				for (tmp1 = 0 ; tmp1 < tmp ; tmp1++ )
				{
					startValue[ nrGlobalVar*tmp1 + factorGloablI] = outValue[fa][tmp1];
					randNr[fa][tmp1] = correlatedRandNum[nrTimeStep*(nrFactor*nrMCScenarios_) + nrFactor*(nrs+tmp1) + fa ];
				}
#if ( defined(FITOB_OPENMP) )
}
#endif
			}
			// use "startValue" value where the values will be changed
			// call the function of "model" which does the MC simulation
			for (int fa = 0 ; fa < nrFactor ; fa++)
			{
				// make the model's forward step
				modelcollection->getModel(fa).simulateForward( tmp , actTime, actDt , nrGlobalVar , startValue , outValue[fa] , randNr[fa] );
			}

			actTime = actTime + actDt;
		// end of time stepping
		}

		// copy the new values bach to the new MC step
		for (int fa = 0 ; fa < nrFactor ; fa++){
			// the first local variables are always the diffusion variables
			// this is the convention
			factorGloablI = modelcollection->getModel(fa).getGlobalIndex();
			variableIndexAct = actMCStep.getVariableIndex( factorGloablI ) ;

			// copy the simulated values back
			for (tmp = 0; (tmp < buffSize) && (nrs < nrMCScenarios_) ; tmp++ , nrs++ )
			{
				actMCStep.setSimulationValue( nrs , variableIndexAct , outValue[fa][tmp] );
			}
			nrs = nrs - tmp;
		}

		nrs = nrs + tmp;
	} // end of while loop
}

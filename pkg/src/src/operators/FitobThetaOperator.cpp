/*
 * FitobThetaOperator.cpp
 *
 *  Created on: Mar 1, 2010
 *      Author: benk
 */

#include "FitobThetaOperator.hpp"

#include "src/diffusionmodel/FitobModelCollection.hpp"
#include "src/diffusionmodel/FitobFactorModelBase.hpp"
#include "src/scripteval/FitobCalculator.hpp"
#include "src/montecarlo/FitobMCMachine.hpp"

#include <vector>
#include <algorithm>

using namespace fitob;

ThetaOperator::ThetaOperator(ExpressionBasis* thetaExpression) :
OperatorBasis("ThetaOperator",  (ThOperator) ) ,
thetaExpression_(thetaExpression) , thetaTime_(0.0){
	setVerb(0);
}

ThetaOperator::~ThetaOperator() {
}

void ThetaOperator::returnScaledNormalDistribution(
		  int level ,
		  double stdFactor ,
		  int scaleIndex ,
		  DVector& output){

/* --- Matlab code ---
	L = 6;
	N = -1:(1/2^L):1;
	Fakt = 5.3;
	V = tan(((pi/2) - 1/Fakt)*N);
	V1 = 3*atan( (N).^5+0.1*(N) );
	V = V ./ max(V);
    V1 = V1 ./ max(V1);
*/

	// first generate linear distribution
	int nrPoints = powerTwo[level] + 1;
	DVector tmpPoints(nrPoints);
	output.resize(nrPoints);

	for (int ii = 0 ; ii < nrPoints; ii++)
		tmpPoints[ii] = (double)(2*ii - powerTwo[level]) / (double)(powerTwo[level]);

	// ---- use different scaling formulas -----
	switch (scaleIndex){
		case 0:{
            double internFactor = 1.0/7.0;
            // calculate grading
            for (int ii = 0 ; ii < nrPoints; ii++)
            	output[ii] = tan( (FITOB_HALF_PI - internFactor)*tmpPoints[ii]);
            // do the scaling
            for (int ii = 0 ; ii < nrPoints; ii++)
            	output[ii] = stdFactor * output[ii] / output[nrPoints-1];
			break;
		}
		case 1:{
			// do some prework in the formula
			for (int ii = 0 ; ii < nrPoints; ii++)
				output[ii] = pow(tmpPoints[ii],5) + 0.1 * tmpPoints[ii];
            // calculate grading
            for (int ii = 0 ; ii < nrPoints; ii++)
            	output[ii] = 3.0 * atan( output[ii] );
            // do the scaling
            for (int ii = 0 ; ii < nrPoints; ii++)
            	output[ii] = stdFactor * output[ii] / output[nrPoints-1];
			break;
		}
	}
}

void ThetaOperator::forwardEstimation(
		boost::ptr_vector<MeshContext>& contextStack ,
		                        FitobCalculator* calc ,
		                        int& stackIndex ,
		                        double& timeStamp ){

	FITOB_OUT_LEVEL3(verb(),"ThetaOperator::forwardEstimation stackIndex:" << stackIndex << "timeStamp"
			<< timeStamp << " contextStack.size():" << contextStack.size());

	MeshContext& actualContext = contextStack[stackIndex];

	if (thetaExpression_->isConstantExpression(actualContext) == false){
		FITOB_ERROR_EXIT(" ThetaOperator, theta expression is not constant ");
	}
    // evaluate the the expression (theta time) and store it
	thetaTime_ = thetaExpression_->eval(actualContext.minGlobCoord());

	// todo: these values now are not taken in consideration,
	// This might be a future feature, to split up one Theta operator in small theta operators
	// in case of too large macro time step
	// 19.09.2010 -> at todays knowledge this is not necessary, this would only make things worse
	maxThetaTime_ = 10.0; //calc->getXMLConfiguration().get()->getDoubleConfiguration("thetaconfigurations.solver.maxThetaStep.<xmlattr>.value" );
	nr_Theta_ = ceil(thetaTime_/maxThetaTime_);

	const boost::shared_ptr<ModelCollection>& factorModels = calc->getModelCollection();
	const boost::shared_ptr<ScriptModel>& scriptModel = calc->getScriptModel();
	const OperatorSequence& scriptBody = scriptModel->bodyOpSeq();
	const Domain& actDom = actualContext.getDom();
	int nrFactors = factorModels->nrFactors();

	Domain newDom(actDom);

	FITOB_OUT_LEVEL3(verb(),"ThetaOperator::forwardEstimation nrFactors:" << factorModels->nrFactors() << " thetaTime:" << thetaTime_ );
	for (int factorIndex = 0 ; factorIndex < nrFactors ; factorIndex++)
	{
		const FactorModelBase& model = factorModels->getModel(factorIndex);
		int globalIndex = model.getGlobalIndex();

		const DVector& startVal = calc->getStartDomain()->getGradedAxis(globalIndex);

		// get the points of the normal distribution
		DVector stdPoints;
		returnScaledNormalDistribution( actDom.getMaximalLevel() ,
				factorModels->stdFactor(factorIndex) , 1 ,  stdPoints );

		FITOB_OUT_LEVEL3(verb(),"ThetaOperator::forwardEstimation factorIndex:" << factorIndex << " P.size():" << stdPoints.size());
		// here we calculate the scaling
		if (startVal.size() > 1){
			for (unsigned int pointI = 0; pointI < stdPoints.size() ; pointI++ ){
				double endVal = 0.0;
				model.forwardEstimation(startVal[pointI] , timeStamp + thetaTime_ ,
						stdPoints[pointI] , actDom.getAverage() , endVal );
				if (verb()>3){
					std::cout << endVal << ",";
				}
				stdPoints[pointI] = endVal;
			}
			if (verb()>3) std::cout << std::endl;
		} else {
			for (unsigned int pointI = 0; pointI < stdPoints.size() ; pointI++ ){
				double endVal = 0.0;
				model.forwardEstimation(startVal[0] , timeStamp + thetaTime_ ,
						stdPoints[pointI] , actDom.getAverage() , endVal );
				if (verb()>3){
					std::cout << endVal << ",";
				}
				stdPoints[pointI] = endVal;
			}
			if (verb()>3) std::cout << std::endl;
		}

		//and call forwardEstimation_DiffusionEnlarement
		int stackIndex_tmp = 0;
		scriptBody.forwardEstimation_DiffusionEnlarement( contextStack ,
                calc ,stdPoints , globalIndex , stackIndex_tmp  , timeStamp );

		// sort stdPoints vector - not necessary (even in mean reversion case)
		// - this is not necessary since we estimate the values from 0 till T and not for dT
		std::sort( stdPoints.begin(),stdPoints.end() );

		// set the axis in the new domain
		newDom.setAxis(stdPoints,globalIndex);
	}

	// add new context, extend contextStack
	contextStack.push_back(new MeshContext(newDom,timeStamp));
	FITOB_OUT_LEVEL3(verb(),"ThetaOperator::forwardEstimation New Domain:" << newDom.toString() );
	timeStamp = timeStamp + thetaTime_;
	stackIndex = stackIndex + 1;

}

void ThetaOperator::forwardEstimation_DiffusionEnlarement(
		                          boost::ptr_vector<MeshContext>& contextStack ,
		                          FitobCalculator* calc ,
		                          DVector& gradValues,
		                          int globalVariableIndex,
			                      int& stackIndex ,
		                          double timeStamp ) const {
	// if the index already below the index then nothing to do
	// this is the exit condition for instructions below "timeStamp"
	if ( stackIndex >= (int)contextStack.size() ) return; // EXIT IF NEEDED

	// here we only increase the stack index
	stackIndex = stackIndex + 1;
}

void ThetaOperator::backwardCalculation(
		boost::ptr_vector<MeshContext>& contextStack ,
		                          FitobCalculator* calc ,
			                      int& stackIndex ,
		                          double& timeStamp ){
	// actions taken in this method
    //       - solve PDE on the old mesh(from FitobCalculator) in the old context, domain (stackIndex)
	//       - create context first with the old mesh (from FitobCalculator) and old Domain (stackIndex)
	//       - create new mesh with the new mesh and with new Domain (stackIndex - 1)
	//       - set the values on the new mesh
	//       - set the new mesh to be the actual mesh in FitobCalculator


	FITOB_OUT_LEVEL1(verb()," ThetaOperator::backwardCalculation  before PDE solve");
	//in case of sparse grids make combination technique, or choose direct solver
	calc->getSolver()->solvePDE((contextStack[stackIndex].getMesh()) ,
			&(contextStack[stackIndex]) ,
			calc->getModelCollection().get() ,
			calc ,
			thetaTime_ );

	FITOB_OUT_LEVEL1(verb()," ThetaOperator::backwardCalculation  after PDE solve");

	// create new Grid
	boost::shared_ptr<MeshBase> newGrid;
	newGrid = calc->getGridFactory()->createMesh( &(contextStack[stackIndex-1].getDom()) , calc );
	// set the grid to be the grid of this grid context
	FITOB_OUT_LEVEL1(verb()," ThetaOperator::backwardCalculation set new values on the mesh");
	contextStack[stackIndex-1].setMesh(newGrid);
	// set the values of the new grid using the values of the old grid
	newGrid.get()->setValues( &(contextStack[stackIndex]) , contextStack[stackIndex].getDom().getAverage() );
	// do eventually plotting
	FITOB_OUT_LEVEL1(verb()," ThetaOperator::backwardCalculation plotting");
	calc->getPlotter()->plot( &(contextStack[stackIndex-1]) , calc->getScriptModel()->getModelName() );
    // decrement the stack index
	stackIndex = stackIndex - 1;
}


void ThetaOperator::forwardMCSimulation( boost::ptr_vector<MeshContext>& contextStack ,
								  FitobCalculator* calc , MCMachine* MC ,
								  int& stackIndex , double& timeStamp ) {
	// create the new step
	FITOB_OUT_LEVEL1(verb()," ThetaOperator::forwardMCSimulation START timeStamp:" << timeStamp << " , thetaTime_:" << thetaTime_
			<< " , stackIndex:" << stackIndex << " contextStack[stackIndex].getTime():" <<contextStack[stackIndex].getTime());
	MC->addMCStep( contextStack[stackIndex].getDom() , timeStamp + thetaTime_ , stackIndex );
	timeStamp = timeStamp + thetaTime_;
	FITOB_OUT_LEVEL1(verb()," ThetaOperator::forwardMCSimulation END timeStamp:" << timeStamp << " , thetaTime_:" << thetaTime_
			<< " , stackIndex:" << stackIndex << " contextStack[stackIndex].getTime():" <<contextStack[stackIndex].getTime());
	stackIndex = stackIndex + 1;
}


void ThetaOperator::forwardMCEvaluation( boost::ptr_vector<MeshContext>& contextStack ,
								  FitobCalculator* calc , MCMachine* MC ,
								  int& stackIndex , double& timeStamp ) {
	// nothing to do only increment the stack index
	MC->getMCStep(stackIndex+1).setExportValues( &(MC->getMCStep(stackIndex)) );
	stackIndex = stackIndex + 1;
	timeStamp = contextStack[stackIndex].getTime() + thetaTime_;
}


void ThetaOperator::backwardMCEvaluation( boost::ptr_vector<MeshContext>& contextStack ,
								  FitobCalculator* calc , MCMachine* MC ,
								  int& stackIndex , double& timeStamp ) {
	// nothing to do only decrement the index
	MC->getMCStep(stackIndex-1).setExportValues( &(MC->getMCStep(stackIndex)) );
	stackIndex = stackIndex - 1;
	timeStamp = contextStack[stackIndex].getTime();
}

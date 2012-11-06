/*
 * FitobLoopOperator.cpp
 *
 *  Created on: Feb 28, 2010
 *      Author: benk
 */

#include "FitobLoopOperator.hpp"

using namespace std;
using namespace fitob;

LoopOperator::LoopOperator(ExpressionBasis* loopCondition,OperatorBasis* loopBody)
: OperatorBasis("LoopOperator" , (LoopOp) ) ,
  loopCondition_(loopCondition) , loopBody_(loopBody){
	setVerb(0);
}

LoopOperator::~LoopOperator() {
}

void LoopOperator::forwardEstimation(
		boost::ptr_vector<MeshContext>& contextStack ,
		                        FitobCalculator* calc ,
		                        int& stackIndex ,
		                        double& timeStamp ){
	// see if the expression in the condition constant is
	// we use the context on the top of the stack
	MeshContext& actualContext = contextStack[stackIndex];
	FITOB_OUT_LEVEL3(verb(),"LoopOperator::forwardEstimation  Test if it is constant");
	bool constCondition_ = loopCondition_->isConstantExpression(actualContext);
	nrLoop_ = 0;
	if (constCondition_)
	{  // expression is constant, we calculate the LOOP count
		nrLoop_ = (int)(loopCondition_->eval(actualContext.minGlobCoord()));
	}
	else{
		FITOB_ERROR_EXIT(" Expression in the LOOP operator must be a constant (invariant) expression ");
	}

	// call the body N times
	for (int ind = 0; ind < nrLoop_ ; ind++){
		FITOB_OUT_LEVEL3(verb(),"LoopOperator::forwardEstimation bef. :" << stackIndex <<
			  " contextStack.size():" << contextStack.size() << " ind:"<<ind << " nrLoop:" << nrLoop_);
		loopBody_->forwardEstimation( contextStack , calc , stackIndex , timeStamp );
		FITOB_OUT_LEVEL3(verb(),"LoopOperator::forwardEstimation aft. :" << stackIndex <<
			  " contextStack.size():" << contextStack.size() << " ind:"<<ind << " nrLoop:" << nrLoop_);
	}
}

/** function for diffusion axis size (enlargement) estimation */
void LoopOperator::forwardEstimation_DiffusionEnlarement(
		                          boost::ptr_vector<MeshContext>& contextStack ,
		                          FitobCalculator* calc ,
		                          DVector& gradValues,
		                          int globalVariableIndex,
			                      int& stackIndex ,
		                          double timeStamp ) const {
	// if the index already below the index then nothing to do
	// this is the exit condition for instructions below "timeStamp"
	if ( stackIndex >= (int)contextStack.size() ) return; // EXIT IF NEEDED

	// call the body N times
	for (int ind = 0; ind < nrLoop_ ; ind++){
		loopBody_->forwardEstimation_DiffusionEnlarement( contextStack , calc , gradValues ,
				globalVariableIndex , stackIndex , timeStamp);
		FITOB_OUT_LEVEL3(verb(),"LoopOperator::forwardEstimation:" << stackIndex <<
			  " contextStack.size():" << contextStack.size() << " ind:"<<ind << " nrLoop:" << nrLoop_);
	}
}

void LoopOperator::backwardCalculation(
		boost::ptr_vector<MeshContext>& contextStack ,
		                          FitobCalculator* calc ,
			                      int& stackIndex ,
		                          double& timeStamp ){
	// call the body N times
	for (int ind = 0; ind < nrLoop_ ; ind++){
		loopBody_->backwardCalculation( contextStack , calc , stackIndex , timeStamp );
	}
}



void LoopOperator::forwardMCSimulation( boost::ptr_vector<MeshContext>& contextStack ,
								  FitobCalculator* calc , MCMachine* MC ,
								  int& stackIndex , double& timeStamp ) {
	// call the body N times
	for (int ind = 0; ind < nrLoop_ ; ind++){
		loopBody_->forwardMCSimulation( contextStack , calc , MC , stackIndex , timeStamp );
	}
}


void LoopOperator::forwardMCEvaluation( boost::ptr_vector<MeshContext>& contextStack ,
								  FitobCalculator* calc , MCMachine* MC ,
								  int& stackIndex , double& timeStamp ) {
	// call the body N times
	for (int ind = 0; ind < nrLoop_ ; ind++){
		loopBody_->forwardMCEvaluation( contextStack , calc , MC , stackIndex , timeStamp );
	}
}


void LoopOperator::backwardMCEvaluation( boost::ptr_vector<MeshContext>& contextStack ,
								  FitobCalculator* calc , MCMachine* MC ,
								  int& stackIndex , double& timeStamp ) {
	// call the body N times
	for (int ind = 0; ind < nrLoop_ ; ind++){
		loopBody_->backwardMCEvaluation( contextStack , calc , MC , stackIndex , timeStamp );
	}
}

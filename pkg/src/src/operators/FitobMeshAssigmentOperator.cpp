/*
 * FitobMeshAssigmentOperator.cpp
 *
 *  Created on: Feb 28, 2010
 *      Author: benk
 */

#include "FitobMeshAssigmentOperator.hpp"
#include "src/scripteval/FitobCalculator.hpp"
#include "src/montecarlo/FitobMCMachine.hpp"

using namespace std;
using namespace fitob;

MeshAssigmentOperator::MeshAssigmentOperator(const Variable* expVariable , ExpressionBasis* assignExpression) :
OperatorBasis("MeshAssigmentOperator",(MeshAssigmentOp) ) ,
 expVariable_(expVariable) , assignExpression_(assignExpression){
  setVerb(0);
}

MeshAssigmentOperator::~MeshAssigmentOperator() {
}

void MeshAssigmentOperator::forwardEstimation(
		boost::ptr_vector<MeshContext>& contextStack ,
		                        FitobCalculator* calc ,
		                        int& stackIndex ,
		                        double& timeStamp ){
	// nothing to do here, no domain enlargement occurs here
 }


void MeshAssigmentOperator::backwardCalculation(
		boost::ptr_vector<MeshContext>& contextStack ,
		                          FitobCalculator* calc ,
			                      int& stackIndex ,
		                          double& timeStamp ){
	FITOB_OUT_LEVEL3(verb(),"MeshAssigmentOperator::backwardCalculation stackIndex:" << stackIndex << " contextStack.size():" << contextStack.size()
			<< "  timeStamp:" << timeStamp );
	FITOB_OUT_LEVEL3(verb(),"MeshAssigmentOperator::backwardCalculation Domain:" << contextStack[stackIndex].getDom().toString() );

    // if we get here then we just have a simple mesh assignment , no if condition will be considered here
	contextStack[stackIndex].getMesh()->setValues(
			assignExpression_ , contextStack[stackIndex].minGlobCoord());

	//apply constraints
	DVector globCoord = contextStack[stackIndex].getDom().getAverage();
	contextStack[stackIndex].getMesh()->applyConstraints( &(calc->getScriptModel().get()->constrainOpSeq()) , globCoord );

	// do eventually plotting
	calc->getPlotter()->plot( &(contextStack[stackIndex]) , calc->getScriptModel()->getModelName() );
}

void MeshAssigmentOperator::applyConstraints(DVector& globalCoords) const {
	// just call the expression and return the result in the array (position 1)back
	// todo: use other function call
	globalCoords[1] = assignExpression_->eval(globalCoords);
}


void MeshAssigmentOperator::forwardMCSimulation( boost::ptr_vector<MeshContext>& contextStack ,
								  FitobCalculator* calc , MCMachine* MC ,
								  int& stackIndex , double& timeStamp ) {
	// --- NO OP ---
}


void MeshAssigmentOperator::forwardMCEvaluation( boost::ptr_vector<MeshContext>& contextStack ,
								  FitobCalculator* calc , MCMachine* MC ,
								  int& stackIndex , double& timeStamp ) {
	// just set the export variable in the current MC step
	(MC->getMCStep(stackIndex)).applyExpression( &(MC->getMCStep(stackIndex)) ,
			expVariable_ , assignExpression_ , calc );
}


void MeshAssigmentOperator::backwardMCEvaluation( boost::ptr_vector<MeshContext>& contextStack ,
								  FitobCalculator* calc , MCMachine* MC ,
								  int& stackIndex , double& timeStamp ) {
	// just set the export variable in the current MC step
	(MC->getMCStep(stackIndex)).applyExpression( &(MC->getMCStep(stackIndex)) ,
			expVariable_ , assignExpression_ , calc );
}

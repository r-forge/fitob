/*
 * FitobRegressionController.cpp
 *
 *  Created on: Jul 7, 2011
 *      Author: benk
 */


#include "FitobRegressionController.hpp"
#include "src/montecarlo/FitobMCStep.hpp"
#include "src/scripteval/FitobCalculator.hpp"

using namespace fitob;

bool RegressionController::canRegressionBeTriggerd_ = false;

bool RegressionController::isRegressionTriggerd_ = false;

Regularization* RegressionController::regClass_ = NULL;

const MCStep* RegressionController::mcStep_ = NULL;

const Domain* RegressionController::dom_ = NULL;

const FitobCalculator* RegressionController::calc_ = NULL;

void RegressionController::startRegression(  const ExpectedExpression* expectExpr , const ExpressionBasis* expEval ){

	if ( !canRegressionBeTriggerd_ || isRegressionTriggerd_ ) return;

	canRegressionBeTriggerd_ = false;
	isRegressionTriggerd_ = true;

	int verb = 6;

	DVector XCoords( mcStep_->getNrScenario() * dom_->nrRealAxis() );
	DVector YCoords( mcStep_->getNrScenario() );
	DVector globalCoords = mcStep_->getAverage();

	// for each scenario copy the coordinates in the array and evaluate the expression
	int d , dim = dom_->nrRealAxis() , varI ,globI;
	FITOB_OUT_LEVEL3( verb ," RegressionController::startRegression , mcStep_->getNrScenario()=" << mcStep_->getNrScenario()
			<< " , dim = " << dim );
	FITOB_OUT_LEVEL3( verb ," RegressionController::startRegression , EXPRESSION=" << expEval->toString() );
	for (int s = 0 ; s < mcStep_->getNrScenario() ; s++ ){
		// for each grid dimension (which are the stochastic variables and additional import variables)
		for (d = 0 ; d < dim ; d++){
			globI = dom_->localToGlobalIndex( d );
			varI = mcStep_->getVariableIndex( globI );
			XCoords[s*dim + d] = mcStep_->getSimulationValue( s , varI );
			globalCoords[globI] = mcStep_->getSimulationValue( s , varI );
		}
		// copy each export variable
		for (d = 0 ; d < dom_->nrExportVariables() ; d++){
			globI = dom_->exportToGlobal(d); // global index of the export variables
			varI = mcStep_->getVariableIndex( globI );
			globalCoords[globI] = mcStep_->getSimulationValue( s , varI );
		}
		// evaluate the expression
		YCoords[s] = expEval->eval( globalCoords );
		//FITOB_OUT_LEVEL3( verb ," RegressionController::startRegression , RES["<<s<<"]=" << YCoords[s] );
	}

	// create the regularization class
	regClass_ = new Regularization( calc_ , XCoords , YCoords , dom_ , expectExpr );
}

void RegressionController::doneRegression(){
	canRegressionBeTriggerd_ = true;
	isRegressionTriggerd_ = false;
	delete regClass_;
}

/*
 * FitobExpExpression.cpp
 *
 *  Created on: Feb 27, 2010
 *      Author: benk
 */

#include "FitobExpectedExpression.hpp"
#include "src/regression/FitobRegressionController.hpp"

using namespace std;
using namespace fitob;

ExpectedExpression::ExpectedExpression(ExpressionBasis* expr ) :
ExpressionBasis("ExpectedExpression") , exp_(expr) , ITriggeredRegression_(false) {
}

ExpectedExpression::~ExpectedExpression() {
}


double ExpectedExpression::eval(const DVector& globalCoordonates) const {

	//FITOB_OUT_LEVEL2( 3 ,"fl1: " << ITriggeredRegression_ << " , fl2:" << RegressionController::canRegressionBeTriggerd() );
	// first check if regression can be triggered, if yes then make the regression
	if ( (!ITriggeredRegression_) && ( RegressionController::canRegressionBeTriggerd() ) ){
		//FITOB_OUT_LEVEL2( 3 ," START REGRESSION " );
		RegressionController::startRegression( this , exp_ );
	}

	// if this expression triggered the regression then evaluate the grid
	if (ITriggeredRegression_){
		double tmp_d = RegressionController::getMesh()->eval(globalCoordonates);
		//FITOB_OUT_LEVEL2( 3 ," EVAL MESH REGRESSION , res = " << tmp_d );
		return tmp_d;
	}
	else {
		return exp_->eval(globalCoordonates);
	}
}

void ExpectedExpression::eval(const std::vector<DVector>& globalCoordonates , DVector& resVector) const {

	// first check if regression can be triggered, if yes then make the regression
	if ( (!ITriggeredRegression_) && ( RegressionController::canRegressionBeTriggerd() ) ){
		RegressionController::startRegression( this , exp_ );
	}

	// if this expression triggered the regression then evaluate the grid
	if (ITriggeredRegression_){
		return RegressionController::getMesh()->eval(globalCoordonates,resVector);
	}
	else {
		exp_->eval(globalCoordonates , resVector );
    }
}

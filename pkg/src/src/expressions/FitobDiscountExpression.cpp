/*
 * FitobExpExpression.cpp
 *
 *  Created on: Feb 27, 2010
 *      Author: benk
 */

#include "FitobDiscountExpression.hpp"
#include "src/scripteval/FitobCalculator.hpp"

using namespace std;
using namespace fitob;

DiscountExpression::DiscountExpression(ExpressionBasis* t1, ExpressionBasis* t2) :
ExpressionBasis("DiscountExpression") , t1_(t1), t2_(t2) {
}

DiscountExpression::~DiscountExpression() {
}


double DiscountExpression::eval(const DVector& globalCoordonates) const {
	// first eval the times, then call the discount factors
	double t1 = t1_->eval(globalCoordonates);
	double t2 = t2_->eval(globalCoordonates);
	return (getFitobCalculator()->getModelCollection())->getDiscountFactor(globalCoordonates,t1,t2);
}


void DiscountExpression::eval(const std::vector<DVector>& globalCoordonates , DVector& resVector) const {
	// first eval the times, then call the discount factors
	DVector t1(resVector.size());
	DVector t2(resVector.size());
	t1_->eval(globalCoordonates,t1);
	t2_->eval(globalCoordonates,t2);
	for (unsigned int ii = 0 ; ii < resVector.size() ; ii++){
		resVector[ii] = (getFitobCalculator()->getModelCollection())->getDiscountFactor( globalCoordonates[ii] , t1[ii] , t2[ii] );
	}
}

bool DiscountExpression::isConstantExpression(const MeshContext& context) const {
	// this expression is constant when both the teniors are constants and when the interest rate is globally constant
	return (t1_->isConstantExpression(context) && t2_->isConstantExpression(context)
			 && ( (getFitobCalculator()->getModelCollection())->getInterestRateGlobalIndex() < 0 ));
}

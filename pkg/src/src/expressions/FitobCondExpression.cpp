/*
 * FitobCondExpression.cpp
 *
 *  Created on: Aug 23, 2011
 *      Author: benk
 */

#include "FitobCondExpression.hpp"

using namespace fitob;

CondExpression::CondExpression(ExpressionBasis* cond , ExpressionBasis* trueExpr ,
ExpressionBasis* falseExpr):
ExpressionBasis("CondExpression") , cond_(cond) , trueExpr_(trueExpr) , falseExpr_(falseExpr){
// nothing to do here more
}

CondExpression::~CondExpression() {
	// nothing to do here
}

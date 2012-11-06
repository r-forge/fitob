/*
 * FitobConstantExpression.cpp
 *
 *  Created on: Feb 26, 2010
 *      Author: benk
 */

#include "FitobConstantExpression.hpp"

using namespace fitob;
using namespace std;

ConstantExpression::ConstantExpression(const double value):
		ExpressionBasis("ConstantExpression "), constValue_(value){
}

ConstantExpression::~ConstantExpression() {
}

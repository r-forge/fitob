/*
 * FitobVariableExpression.cpp
 *
 *  Created on: Feb 27, 2010
 *      Author: benk
 */

#include "FitobVariableExpression.hpp"

using namespace std;
using namespace fitob;

VariableExpression::VariableExpression(const Variable* variable ,  double signVar , bool MClookback)
:ExpressionBasis("VariableExpression") , variable_(variable), sign_(signVar), MClookback_(MClookback){
}

VariableExpression::~VariableExpression() {
}

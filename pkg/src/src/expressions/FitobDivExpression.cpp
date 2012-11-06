/*
 * FitobDivExpression.cpp
 *
 *  Created on: Feb 27, 2010
 *      Author: benk
 */

#include "FitobDivExpression.hpp"

using namespace std;
using namespace fitob;

DivExpression::DivExpression(ExpressionBasis* left , ExpressionBasis* right) :
ExpressionBasis("DivExpression") , left_(left) , right_(right) {
}

DivExpression::~DivExpression() {
}

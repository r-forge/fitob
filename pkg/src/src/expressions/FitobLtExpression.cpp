/*
 * FitobLtExpression.cpp
 *
 *  Created on: Feb 27, 2010
 *      Author: benk
 */

#include "FitobLtExpression.hpp"

using namespace std;
using namespace fitob;

LtExpression::LtExpression(ExpressionBasis* left , ExpressionBasis* right) :
ExpressionBasis("LtExpression") , left_(left) , right_(right) {
}

LtExpression::~LtExpression() {
}

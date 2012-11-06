/*
 * FitobEqExpression.cpp
 *
 *  Created on: Feb 27, 2010
 *      Author: benk
 */

#include "FitobEqExpression.hpp"

using namespace std;
using namespace fitob;

EqExpression::EqExpression(ExpressionBasis* left , ExpressionBasis* right) :
ExpressionBasis("EqExpression") , left_(left) , right_(right) {
}

EqExpression::~EqExpression() {
}

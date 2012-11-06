/*
 * FitobPowExpression.cpp
 *
 *  Created on: Feb 27, 2010
 *      Author: benk
 */

#include "FitobPowExpression.hpp"

using namespace std;
using namespace fitob;

PowExpression::PowExpression(ExpressionBasis* left , ExpressionBasis* right) :
ExpressionBasis("PowExpression") , left_(left) , right_(right) {
}

PowExpression::~PowExpression() {
}

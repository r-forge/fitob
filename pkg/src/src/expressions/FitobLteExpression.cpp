	/*
 * FitobLteExpression.cpp
 *
 *  Created on: Feb 27, 2010
 *      Author: benk
 */

#include "FitobLteExpression.hpp"

using namespace std;
using namespace fitob;

LteExpression::LteExpression(ExpressionBasis* left , ExpressionBasis* right) :
ExpressionBasis("LteExpression") , left_(left) , right_(right) {
}

LteExpression::~LteExpression() {
}

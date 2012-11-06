/*
 * FitobGtExpression.cpp
 *
 *  Created on: Feb 27, 2010
 *      Author: benk
 */

#include "FitobGtExpression.hpp"

using namespace std;
using namespace fitob;

GtExpression::GtExpression(ExpressionBasis* left , ExpressionBasis* right) :
ExpressionBasis("GtExpression") , left_(left) , right_(right) {
}

GtExpression::~GtExpression() {
}

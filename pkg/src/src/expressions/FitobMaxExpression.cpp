/*
 * FitobMaxExpression.cpp
 *
 *  Created on: Feb 27, 2010
 *      Author: benk
 */

#include "FitobMaxExpression.hpp"

using namespace std;
using namespace fitob;

MaxExpression::MaxExpression(ExpressionBasis* left , ExpressionBasis* right) :
ExpressionBasis("MaxExpression") , left_(left) , right_(right) {
}

MaxExpression::~MaxExpression() {
}

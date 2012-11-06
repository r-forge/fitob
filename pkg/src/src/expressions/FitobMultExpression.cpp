/*
 * FitobMultExpression.cpp
 *
 *  Created on: Feb 27, 2010
 *      Author: benk
 */

#include "FitobMultExpression.hpp"

using namespace std;
using namespace fitob;

MultExpression::MultExpression(ExpressionBasis* left , ExpressionBasis* right) :
ExpressionBasis("MultExpression") , left_(left) , right_(right) {
}

MultExpression::~MultExpression() {
}

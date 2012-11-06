/*
 * FitobGteExpression.cpp
 *
 *  Created on: Feb 27, 2010
 *      Author: benk
 */

#include "FitobGteExpression.hpp"

using namespace std;
using namespace fitob;

GteExpression::GteExpression(ExpressionBasis* left , ExpressionBasis* right) :
ExpressionBasis("GteExpression") , left_(left) , right_(right) {
}

GteExpression::~GteExpression() {
}

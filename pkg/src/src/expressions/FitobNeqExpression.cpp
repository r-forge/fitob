/*
 * FitobNeqExpression.cpp
 *
 *  Created on: Feb 27, 2010
 *      Author: benk
 */

#include "FitobNeqExpression.hpp"

using namespace std;
using namespace fitob;

NeqExpression::NeqExpression(ExpressionBasis* left , ExpressionBasis* right) :
ExpressionBasis("NeqExpression") , left_(left) , right_(right) {
}

NeqExpression::~NeqExpression() {
}

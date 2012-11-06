/*
 * FitobMinExpression.cpp
 *
 *  Created on: Feb 27, 2010
 *      Author: benk
 */

#include "FitobMinExpression.hpp"

using namespace std;
using namespace fitob;

MinExpression::MinExpression(ExpressionBasis* left , ExpressionBasis* right) :
ExpressionBasis("MinExpression") , left_(left) , right_(right) {
}

MinExpression::~MinExpression() {
}

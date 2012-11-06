/*
 * FitobAndExpression.cpp
 *
 *  Created on: Feb 27, 2010
 *      Author: benk
 */

#include "FitobAndExpression.hpp"

using namespace std;
using namespace fitob;

AndExpression::AndExpression(ExpressionBasis* left , ExpressionBasis* right) :
ExpressionBasis("AndExpression") , left_(left) , right_(right) {
}

AndExpression::~AndExpression() {
}

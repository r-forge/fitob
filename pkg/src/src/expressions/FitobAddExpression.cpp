/*
 * FitobAddExpression.cpp
 *
 *  Created on: Feb 27, 2010
 *      Author: benk
 */

#include "FitobAddExpression.hpp"

using namespace std;
using namespace fitob;

AddExpression::AddExpression(ExpressionBasis* left , ExpressionBasis* right) :
ExpressionBasis("AddExpression") , left_(left) , right_(right) {
}

AddExpression::~AddExpression() {
}

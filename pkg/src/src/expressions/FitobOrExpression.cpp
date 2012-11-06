/*
 * FitobOrExpression.cpp
 *
 *  Created on: Feb 27, 2010
 *      Author: benk
 */

#include "FitobOrExpression.hpp"

using namespace std;
using namespace fitob;

OrExpression::OrExpression(ExpressionBasis* left , ExpressionBasis* right) :
ExpressionBasis("OrExpression") , left_(left) , right_(right) {
}

OrExpression::~OrExpression() {
}

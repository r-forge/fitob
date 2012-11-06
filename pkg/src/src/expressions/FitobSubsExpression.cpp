/*
 * FitobSuEbsxpression.cpp
 *
 *  Created on: Feb 27, 2010
 *      Author: benk
 */

#include "FitobSubsExpression.hpp"

using namespace std;
using namespace fitob;

SubsExpression::SubsExpression(ExpressionBasis* left , ExpressionBasis* right) :
ExpressionBasis("SubsExpression") , left_(left) , right_(right) {
}

SubsExpression::~SubsExpression() {
}

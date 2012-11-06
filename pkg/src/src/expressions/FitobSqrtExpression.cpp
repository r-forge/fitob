/*
 * FitobSqrtExpression.cpp
 *
 *  Created on: Feb 27, 2010
 *      Author: benk
 */

#include "FitobSqrtExpression.hpp"

using namespace std;
using namespace fitob;

SqrtExpression::SqrtExpression(ExpressionBasis* expr ) :
ExpressionBasis("SqrtExpression") , exp_(expr) {
}

SqrtExpression::~SqrtExpression() {
}

/*
 * FitobExpExpression.cpp
 *
 *  Created on: Feb 27, 2010
 *      Author: benk
 */

#include "FitobExpExpression.hpp"

using namespace std;
using namespace fitob;

ExpExpression::ExpExpression(ExpressionBasis* expr ) :
ExpressionBasis("ExpExpression") , exp_(expr) {
}

ExpExpression::~ExpExpression() {
}

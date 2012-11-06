/*
 * FitobLogExpression.cpp
 *
 *  Created on: Feb 27, 2010
 *      Author: benk
 */

#include "FitobLogExpression.hpp"

using namespace std;
using namespace fitob;

LogExpression::LogExpression(ExpressionBasis* expr ) :
ExpressionBasis("LogExpression") , exp_(expr) {
}

LogExpression::~LogExpression() {
}

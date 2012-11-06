/*
 * FitobExpressionBasis.cpp
 *
 *  Created on: Feb 25, 2010
 *      Author: benk
 */

#include "FitobExpressionBasis.hpp"
#include "src/scripteval/FitobCalculator.hpp"

using namespace fitob;

// init with null
FitobCalculator* ExpressionBasis::fitobCalc_ = 0;

ExpressionBasis::ExpressionBasis():name_("NoNameExpression") {
}

ExpressionBasis::ExpressionBasis(const string& name): name_(name) {
}

ExpressionBasis::~ExpressionBasis() {
}

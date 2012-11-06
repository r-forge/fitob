/*
 * FitobIfElseAssigmentExpression.cpp
 *
 *  Created on: Jul 20, 2010
 *      Author: benk
 */

#include "FitobIfElseAssigmentExpression.hpp"

using namespace fitob;
using namespace std;

IfElseAssigmentExpression::IfElseAssigmentExpression(
		const ExpressionBasis* cond ,
        const ExpressionBasis* trueB ,
        const ExpressionBasis* falseB) :  ExpressionBasis("FitobIfElseAssigmentExpression") ,
        cond_(cond) , true_(trueB) , false_(falseB) , hasElseBanch_(true) {

   if (falseB == NULL) hasElseBanch_ = false;

   setVerb(6);

}

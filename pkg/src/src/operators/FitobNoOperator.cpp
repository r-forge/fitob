/*
 * FitobNoOperator.cpp
 *
 *  Created on: Mar 1, 2010
 *      Author: benk
 */

#include "FitobNoOperator.hpp"

using namespace fitob;
using namespace std;

NoOperator::NoOperator() :
 OperatorBasis("NoOperator",  (NoOp) ) {
}

NoOperator::~NoOperator() {
}

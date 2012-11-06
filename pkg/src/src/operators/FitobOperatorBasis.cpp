/*
 * FitobOperatorBasis.cpp
 *
 *  Created on: Feb 25, 2010
 *      Author: benk
 */

#include "FitobOperatorBasis.hpp"

using namespace fitob;
using namespace std;

OperatorBasis::OperatorBasis() : operationName_("NoNameOperator") ,
		operatorType_( (OperatorType)6) , operatorTime_(0.0)
{
	setVerb(6);
}

OperatorBasis::OperatorBasis(const string& name, const OperatorType& opType)
: operationName_(name) , operatorType_(opType) , operatorTime_(0.0)
{
	setVerb(6);
}

OperatorBasis::~OperatorBasis()
{
}

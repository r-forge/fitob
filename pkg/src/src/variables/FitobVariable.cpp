/*
 * FitobVariable.cpp
 *
 *  Created on: Feb 26, 2010
 *      Author: benk
 */

#include "FitobVariable.hpp"

using namespace fitob;
using namespace std;

Variable::Variable()
: name_("NOVAR"), type_(Export) , globalIndex_(0){
}

Variable::Variable(const string& name, const VariableType type, const int globalIndex)
: name_(name), type_(type) , globalIndex_(globalIndex){
	setVerb(0);
	FITOB_OUT_LEVEL3(verb(),"Variable::Variable name:" << name << " glbIndex:" << globalIndex );
}

Variable::~Variable() {
}

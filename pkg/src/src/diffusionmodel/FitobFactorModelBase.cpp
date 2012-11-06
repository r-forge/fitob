/*
 * FitobFactorModelBase.cpp
 *
 *  Created on: Apr 13, 2010
 *      Author: benk
 */

#include "FitobFactorModelBase.hpp"
#include "FitobModelCollection.hpp"

using namespace fitob;
using namespace std;

ModelCollection* FactorModelBase::modelCollection_ = NULL;

FactorModelBase::FactorModelBase(const Variable* variable , int factorIndex) :
    variable_(variable) , factorIndex_(factorIndex){

	setVerb(4);
	//FITOB_OUT_LEVEL3(5,"globalIndex:" << variable_->getGlobalIndex() << " factorIndex:" << factorIndex);
}

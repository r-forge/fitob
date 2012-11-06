/*
 * FitobModelCollection.cpp
 *
 *  Created on: Apr 13, 2010
 *      Author: benk
 */

#include "FitobModelCollection.hpp"

#include "src/diffusionmodel/FitobFactorModelGeomBr.hpp"
#include "src/diffusionmodel/FitobFactorModelNormBr.hpp"
#include "src/diffusionmodel/FitobFactorModelHullWhite.hpp"
#include "src/diffusionmodel/FitobFactorModelHestonCIR.hpp"
#include "src/diffusionmodel/FitobFactorModelHeston.hpp"
#include "src/diffusionmodel/FitobFactorModelGeneral.hpp"

using namespace fitob;
using namespace std;

int ModelCollection::lastGlobalIndex_ = -1;
int ModelCollection::lastFactorIndex_ = -1;

ModelCollection::ModelCollection(
		const boost::shared_ptr<XMLConfiguration> &ptr ,
		boost::shared_ptr<ScriptModel> &scriptModel )
{

	setVerb(3);
	// here we have to build the array of Factor models

    // make the model collection for all factor visible
	FactorModelBase::setModelCollect( this );

	// get the number of underlyings
	nrFactors_ = ptr->nrXMLNodes("thetaconfigurations.thetaoperator.DIFFUSION_VARIABLES");

	/* Not all variables will appear in the script (e.g. underlying volatility factor)
	   those have to be added to the global variable list */

	factorNames_.resize(nrFactors_);
	factorTypes_.resize(nrFactors_);
    for (int i=0 ; i < nrFactors_; i++)
    {
    	const string tmpStr = ptr->getStringConfiguration("thetaconfigurations.thetaoperator.DIFFUSION_VARIABLES",i,"<xmlattr>.type");
    	// get the global index of the factor as a variable (if is not yet in the list then will be added)
    	// in case of underlying volatility process the factor might not be in the global list
    	factorNames_[i] = ptr->getStringConfiguration("thetaconfigurations.thetaoperator.DIFFUSION_VARIABLES",i,"<xmlattr>.name");
    	factorTypes_[i] = tmpStr;
    	int scriptIndex = -1;
    	// check the factor models in "scriptModel"
    	for (int tmpF = 0 ; tmpF < scriptModel->getNrDefinedModel() ; tmpF++){
    		// we only test the name of the variables with the names which were defined in the script
    		if ( scriptModel->getDefinedModel(tmpF)->getVariable()->getVariableName() == factorNames_[i] ){
    			scriptIndex = tmpF;
    		}
    	}
    	const Variable* var = scriptModel->addVariableVar(factorNames_[i]) ;
        // add the generated factor model to the end of the vector
    	if (scriptIndex > -1) {
    		// use the model defined in the script, and insert it with new index into the factor collection
    		modelContainer_.insert( modelContainer_.end() ,
    			(FactorModelBase*) (new FactorModelGeneral((FactorModelGeneral*)scriptModel->getDefinedModel(scriptIndex),i)) );
    	}
    	else {
    		// create the factor model as it is defined in the XML file
    		modelContainer_.insert( modelContainer_.end() ,
    			makeFactorModel( tmpStr , var , i , scriptModel , ptr) );
    	}
    	ModelCollection::lastGlobalIndex_ = var->getGlobalIndex();
    	ModelCollection::lastFactorIndex_ = i;
    }

	// and read in the correlation matrix
	correlations_.resize(nrFactors_ * nrFactors_);
	standardDeviationFactors_.resize(nrFactors_);
	for (int i = 0 ; i < nrFactors_ ; i++)
	{
		// read in the standard deviation factor for forward estimation
		DVector one_correlation_line;
		standardDeviationFactors_[i] =
				ptr->getDoubleConfiguration("thetaconfigurations.thetaoperator.STANDARD_DEVIATION_FACTORS", i ,"<xmlattr>.value");

		FITOB_OUT_LEVEL3( verb() , " ModelCollection factor=" << i << " , standardDeviationFactors_[i] = " << standardDeviationFactors_[i] );

		// read in one row in the correlation matrix
		ptr->getDoubleVectorConfiguration("thetaconfigurations.thetaoperator.CORRELATIONS", i ,
	                                   "<xmlattr>.value", ',' , one_correlation_line);
		// copy the row from the correlation matrix
		for (int j = 0 ; j < nrFactors_ ; j++){
			correlations_[i*nrFactors_ + j ] = one_correlation_line[j];
			FITOB_OUT_LEVEL3( verb() , " correlation [" << i << "," << j << "]=" << one_correlation_line[j] );
		}
	}

	// read in the risk free interest rate (average value in the case of coupled interest rate)
	r_ = ptr->getDoubleConfiguration("thetaconfigurations.thetaoperator.RISK_FREE_RATE.<xmlattr>.value" );

	// read if the interest rate is coupled to a variable or not
	string variable_coupled_to_interestrate =
		ptr->getStringConfiguration("thetaconfigurations.thetaoperator.RISK_FREE_RATE.<xmlattr>.variable_coupled");

	// initialized to -1, this is default constant interest rate
	r_index_ = -1;
	// check it there is any interest rate model defined in the script
	if ( scriptModel->hasInterestRateModel() ){
		for (int i=0 ; i < nrFactors_; i++)
		{
			// get the interest rate model from the script definition
			if ( factorNames_[i] == scriptModel->getInterestRateModel()->getVariable()->getVariableName() )
			{
				r_index_ = i; // store the factor index of the interest rate
				intRateVariable_ = modelContainer_[i].getVariable();
				FITOB_OUT_LEVEL3( verb() , " found interest rate coupled to variable = " << factorNames_[i] << " global Index=" << r_index_);
				break;
			}
		}
	}
	else
	{
		// if the script does not contain any interest rate models and the XML has one
		if ( variable_coupled_to_interestrate != ""){
			// if this string is not empty then look for the string which should represent one axis
			for (int i=0 ; i < nrFactors_; i++)
			{
				if ( factorNames_[i] == variable_coupled_to_interestrate)
				{
					r_index_ = i; // store the factor index of the interest rate
					intRateVariable_ = modelContainer_[i].getVariable();
					FITOB_OUT_LEVEL3( verb() , " found interest rate coupled to variable = " << factorNames_[i] << " global Index=" << r_index_);
					break;
				}
			}
			// if r_index stays -1 then just use the value
		}
	}
}

void ModelCollection::completeScriptInformation(
		const boost::shared_ptr<XMLConfiguration> &ptr,
		boost::shared_ptr<ScriptModel> &scriptModel){

	   // set the models which are defined in the script
	   // and the interest rate model

	   FITOB_OUT_LEVEL3( verb() , " ModelCollection::completeScriptInformation ");

	   for (int i=0 ; i < nrFactors_; i++)
	    {
	    	int scriptIndex = -1;
	    	// check the factor models in "scriptModel"
	    	FITOB_OUT_LEVEL3( verb() , " ModelCollection::completeScriptInformation  test factor = " << factorNames_[i]);
	    	FITOB_OUT_LEVEL3( verb() , " ModelCollection::completeScriptInformation  scriptModel->getNrDefinedModel() = " << scriptModel->getNrDefinedModel());
	    	for (int tmpF = 0 ; tmpF < scriptModel->getNrDefinedModel() ; tmpF++){
	    		FITOB_OUT_LEVEL3( verb() , " ModelCollection::completeScriptInformation  tmpF = " << tmpF);
	    		// we only test the name of the variables with the names which were defined in the script
	    		FITOB_OUT_LEVEL3( verb() , " ModelCollection::completeScriptInformation  var = " << scriptModel->getDefinedModel(tmpF));
	    		if ( scriptModel->getDefinedModel(tmpF)->getVariable()->getVariableName() == factorNames_[i] ){
	    			scriptIndex = tmpF;
	    		}
	    	}
	    	const Variable* var = scriptModel->addVariableVar(factorNames_[i]) ;
	        // add the generated factor model to the end of the vector
	    	if (scriptIndex > -1) {
	    		// use the model defined in the script, and insert it with new index into the factor collection
	    		FITOB_OUT_LEVEL3( verb() , " ModelCollection::completeScriptInformation set " << factorNames_[i] << " to a user defined model ");
	    		// todo: memory leak ?
	    		modelContainer_.replace(i ,
	                  (FactorModelBase*) (new FactorModelGeneral((FactorModelGeneral*)scriptModel->getDefinedModel(scriptIndex),i)) );
	    	}
	    	ModelCollection::lastGlobalIndex_ = var->getGlobalIndex();
	    	ModelCollection::lastFactorIndex_ = i;
	    }

		// check it there is any interest rate model defined in the script
		if ( scriptModel->hasInterestRateModel() ){
			for (int i=0 ; i < nrFactors_; i++)
			{
				// get the interest rate model from the script definition
				if ( factorNames_[i] == scriptModel->getInterestRateModel()->getVariable()->getVariableName() )
				{
					r_index_ = i; // store the factor index of the interest rate
					intRateVariable_ = modelContainer_[i].getVariable();
					FITOB_OUT_LEVEL3( verb() , " found interest rate coupled to variable = " << factorNames_[i] << " global Index=" << r_index_);
					break;
				}
			}
		}
}

FactorModelBase* ModelCollection::makeFactorModel(
		const string& typeName ,
		const Variable* var ,
		int factorIndex,
		boost::shared_ptr<ScriptModel> &scriptModel ,
		const boost::shared_ptr<XMLConfiguration> &ptr ) {

	if (typeName == "GB"){
		// create geometrical Browian motion model
	   return (FactorModelBase*)(new FactorModelGeomBr(ptr, var ,factorIndex));
	}
	if (typeName == "NB"){
		// create lognormal Browian motion model
	   return (FactorModelBase*)(new FactorModelNormBr(ptr, var ,factorIndex));
	}
	if (typeName == "Hull-White"){
		// create Hull-White model
	   return (FactorModelBase*)(new FactorModelHullWhite(ptr, var ,factorIndex));
	}
	if (typeName == "Heston"){
		// create Hull-White model
	   return (FactorModelBase*)(new FactorModelHeston( ptr ,  var , factorIndex ,
			   scriptModel->getVariable(ModelCollection::lastGlobalIndex_) , ModelCollection::lastFactorIndex_ ));
	}
	if (typeName == "Heston-CIR"){
	   // create Hull-White volatility model
	   return (FactorModelBase*)(new FactorModelHestonCIR(ptr, var ,factorIndex));
	}
	return 0; //Wall
}

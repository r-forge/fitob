/*
 * FitobMCStep.cpp
 *
 *  Created on: Mar 14, 2011
 *      Author: benk
 */

#include "FitobMCStep.hpp"
#include "src/regression/FitobRegressionController.hpp"
#include <boost/lexical_cast.hpp>

using namespace fitob;
using namespace std;

MCStep::MCStep( const Domain& domain , int nrScearions, double time, int evalContextIndex ):
		nrScearions_(nrScearions) ,
		nrVariables_(domain.nrRealAxis() + domain.nrExportVariables()),
		domain_(domain) , evalContextIndex_(evalContextIndex) {

	setVerb(0);
	FITOB_OUT_LEVEL2(verb()," MCStep::MCStep nrVariables_:" << nrVariables_);
	FITOB_OUT_LEVEL2(verb()," Domain : " << domain.toString());

	constGlobalCoordinates_ = domain_.getAverage();
	constGlobalCoordinates_[0] = time;

	// create the simulation vector
	simulationVector_.resize( nrScearions_* nrVariables_ , 0.0 );

	// step has no valid mash
	hasValidMesh_ = false;

	varIndex_to_globalIndex_.resize(domain_.nrRealAxis() + domain_.nrExportVariables(), -1 );
	globalIndex_to_varIndex_.resize(domain_.nrGlobalVariables(), -1 );

	// set the mapping between global index an variable index
	for (int f = 0 ; f < domain_.nrExportVariables() ; f++ ){
		varIndex_to_globalIndex_[f] = f+1;
		globalIndex_to_varIndex_[f+1] = f;
		FITOB_OUT_LEVEL2(verb(), f+1 << " <--> " << f);
	}
	for (int f = domain_.nrExportVariables() ; f < domain_.nrExportVariables()+domain_.nrRealAxis() ; f++ ){
		varIndex_to_globalIndex_[f] = domain_.localToGlobalIndex( f - domain_.nrExportVariables());
		globalIndex_to_varIndex_[domain_.localToGlobalIndex( f - domain_.nrExportVariables())] = f;
		FITOB_OUT_LEVEL2(verb(), domain_.localToGlobalIndex( f - domain_.nrExportVariables()) << " <--> " << f);
	}

	// set all the local variables with their average value
	for (int v = 0 ; v <  nrVariables_ ; v++ ){
		int indexGl = getGlobalIndex(v);
		double avrgVal = constGlobalCoordinates_[indexGl];
		FITOB_OUT_LEVEL2(verb(),"v: " << v << " , indexGl:" << indexGl << " , avrgVal:" << avrgVal);
		for (int i = 0 ; i < nrScearions_ ; i++){
			setSimulationValue(i , v , avrgVal );
		}
	}
}

void MCStep::calculateAverage(){

// todo: use openMP

	DVector average( nrVariables_ , 0.0);

	// calculate the average value
	for (int i = 0 ; i < nrScearions_ ; i++){
		for (int v = 0 ; v <  nrVariables_ ; v++ ){
			average[v] = average[v] + getSimulationValue(i , v );
		}
	}

	// the first are the export variables
	int e = 0;
	for ( ; e < domain_.nrExportVariables() ; e++){
		constGlobalCoordinates_[1+e] = average[e]/(double)nrScearions_;
		//FITOB_OUT_LEVEL4( verb() ,"MCStep::calculateAverage constGlobalCoordinates_[1+e] = " << constGlobalCoordinates_[1+e] );
	}

	// here we set the local variables (axis in the PDE context)
	for ( int l = 0 ; l < domain_.nrRealAxis() ; l++){
		constGlobalCoordinates_[ domain_.localToGlobalIndex(l) ] = average[e+l]/(double)nrScearions_;
	}
}


void MCStep::applyExpression(const MCStep* sourceMC ,
		const Variable* targetVar , const ExpressionBasis* expr ,
		const FitobCalculator* calc ){

	FITOB_ERROR_TEST( this->getNrScenario() == sourceMC->getNrScenario(), " MCStep::applyExpression , nr. of scenarios must agree IN:"
			 << sourceMC->getNrScenario() << " , this:" << this->getNrScenario() );

	DVector glob_tmp = constGlobalCoordinates_;

	int s , var , targetVarGlIndex = targetVar->getGlobalIndex() , varIndexMC;
	double res;

	// trigger for the regression the signal that it should recompute
	RegressionController::triggerableRegression( this , &this->getDom() , calc );

	//todo: if this will go parallel than first make a serial call with 0.0 coordinates, to that regression
	// will be serially trigerd and so can be solved parallel
	res = expr->eval(glob_tmp);

	// get the variable index of the target variable
	varIndexMC = getVariableIndex( targetVarGlIndex );
	for (s = 0 ; s < getNrScenario() ; s++ )
	{
		// for each scenario
		for (var = 0 ; var < sourceMC->nrVariables() ; var++){
			// set the global variable vector
			glob_tmp[ sourceMC->getGlobalIndex( var ) ] = sourceMC->getSimulationValue( s , var );
			//FITOB_OUT_LEVEL4(verb()," glob_tmp[" <<sourceMC->getGlobalIndex( var )<< "]="<< glob_tmp[ sourceMC->getGlobalIndex( var ) ]);
		}
		// calculate and set the value
		res = expr->eval(glob_tmp);
		//FITOB_OUT_LEVEL4( verb() ," res= " << res << " expression:" << expr->toString());
		setSimulationValue( s , varIndexMC , res );
	}

	// todo: before we delete the regression results we might plot the grid

	// tell the regression, that it can drop the actual results
	RegressionController::doneRegression();

	// calculate average values
	this->calculateAverage();

	FITOB_OUT_LEVEL4(verb()," MCStep::applyExpression AFTER :" << this->toString() );
}


void MCStep::setValues(const MCStep* sourceMC ){

	const Domain& sourceDom = sourceMC->getDom();

	// copy the variables, which are different for scenarios
	int v , sourceVarIndex , s;
	for (v = 0 ; v < this->nrVariables() ; v++){
		int globalIndex = this->getGlobalIndex(v);

		sourceVarIndex = sourceMC->getVariableIndex(globalIndex);

		// see if the variable is a constant
		if (sourceVarIndex < 1){
			// this is constant variable
			// get average value
			double avrg = sourceDom.getAverage()[globalIndex];
			for (s = 0 ; s < getNrScenario() ; s++ ){
				setSimulationValue( s , v , avrg );
			}
		}else{
			// it is also variable in the source MCSTep
			for (s = 0 ; s < getNrScenario() ; s++ ){
				setSimulationValue( s , v , sourceMC->getSimulationValue( s , sourceVarIndex ) );
			}
		}
	}

	// calculate average values
	this->calculateAverage();
}


void MCStep::setExportValues(const MCStep* sourceMC ){

	// copy the variables, which are different for scenarios
	int v , sourceVarIndex , s;
	for (v = 0 ; v < this->getDom().nrExportVariables() ; v++){
		int globalIndex = this->getGlobalIndex(v);

		sourceVarIndex = sourceMC->getVariableIndex(globalIndex);

		// it is also variable in the source MCSTep
		for (s = 0 ; s < getNrScenario() ; s++ ){
			setSimulationValue( s , v , sourceMC->getSimulationValue( s , sourceVarIndex ) );
			//FITOB_OUT_LEVEL4(verb()," s= " << s << " v=" << v << " res:" << sourceMC->getSimulationValue( s , sourceVarIndex ) );
		}
	}

	// calculate average values
	this->calculateAverage();
}


const string MCStep::toString() const {
	string ret = "\n  MCStep nrPath:" + boost::lexical_cast<std::string>(getNrScenario()) +
			" , nrVar: " + boost::lexical_cast<std::string>(this->nrVariables()) + " \n";
	int var,s;
	for (s = 0 ; s < getNrScenario() ; s++ ){
		for (var = 0 ; var < this->nrVariables() ; var++){
			ret = ret + boost::lexical_cast<std::string>(this->getSimulationValue( s , var )) + ",";
		}
		ret = ret + "\n";
	}
	return ret;
}

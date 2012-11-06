/*
 * FitobCalculator.cpp
 *
 *  Created on: Apr 15, 2010
 *      Author: benk
 */

#include "FitobCalculator.hpp"
#include "src/expressions/FitobExpressionBasis.hpp"

#include <ctime>
#include <stdio.h>
#include <time.h>

using namespace fitob;
using namespace std;

int FitobCalculator::norm_level_measure = 3;

FitobCalculator::FitobCalculator(
		const string& configurationXMLFile ,
        const string& scriptFile ,
        const int inputLevel ,
        const int convergenceRunNr) {

	setVerb(6);

	// - create XML configuration -
	xmlconfiguration_ = boost::shared_ptr<XMLConfiguration>(new XMLConfiguration(configurationXMLFile));

	// - create script model , but do not parse the script (since it is important that
	// the variables (the factors) will be in the same order as the global variables)
	scriptparser_ =  boost::shared_ptr<ScriptParser>(new ScriptParser(scriptFile));
	scriptmodel_ = boost::shared_ptr<ScriptModel>(new ScriptModel(xmlconfiguration_));

	// - create factor model
   	modelcollection_ = boost::shared_ptr<ModelCollection>(new ModelCollection( xmlconfiguration_ , scriptmodel_ ));

   	// - only now parse the script
   	// so that the factor will be in the same order as the global variables !!!
   	scriptmodel_->parseScript(scriptparser_->getParseInfo());

   	// this method needs to be called, since the script is now complete and it might contain model definitions
   	modelcollection_->completeScriptInformation( xmlconfiguration_ , scriptmodel_ );

   	// - create grid Factory
   	gridFactory_ = boost::shared_ptr<GridFactory>( new GridFactory(xmlconfiguration_) );

   	// - generate the PDE solver
   	pdesolver_ = SolverBase::generateSolver(xmlconfiguration_.get());

   	// - generate the plotter object
   	gridPlotter_ = boost::shared_ptr<GridPlotter>( new GridPlotter(xmlconfiguration_));

   	// ----  set the Calculator for the expression evaluation  ---- , IMPORTANT for discount factor, for expected value
   	ExpressionBasis::setFitobCalculator(this);

   	// read in the required level

   	// ==== DETERMINING THE LEVEL OF THE COMPUTATIONS =============
   	// IDEA: in the forward estimation use the largest level (caz it is cheap), but in the backward calculation
   	// overwrite the level and scaling vector for axis where we have a lower level
   	dimensionAdaptive_ = false;
   	dimensionAdaptiveLevels_.resize(0);

   	int globalLevel =
   			xmlconfiguration_->getIntConfiguration("thetaconfigurations.gridproperties.MAX_LEVEL.<xmlattr>.value");
   	FITOB_OUT_LEVEL3(verb(),"FitobCalculator input globalLevel:" << globalLevel << " inputLevel:" << inputLevel);
   	// if we have an input level then set the global level to the input level (for convergence analysis purpose)
   	if (inputLevel > 0) globalLevel = inputLevel;

	int maxV = 0;

	// if we have a manual dimension adaptive case then read in the level vector
    if ( (convergenceRunNr < 0) &&
    	 ("true" == xmlconfiguration_->getStringConfiguration("thetaconfigurations.gridproperties.DIMENSION_ADPTIVITY.<xmlattr>.value") ) )
    {
       dimensionAdaptive_ = true;
	   xmlconfiguration_->getIntVectorConfiguration("thetaconfigurations.gridproperties.DIMENSION_ADPTIVITY.<xmlattr>.level-vect",
			                     ',' , dimensionAdaptiveLevels_);

	   // set for the factory the adaptive truncation level, in case of T-CT
	   DVector adaptTruncation;
	   xmlconfiguration_->getDoubleVectorConfiguration("thetaconfigurations.gridproperties.combigrid.COMBI_DIAGONAL_CUT_LEVEL.<xmlattr>.vector",
			              ',' , adaptTruncation);
	   if (adaptTruncation.size() > 0)
		   gridFactory_->setDimensionAdaptiveTuncationLevel(adaptTruncation);

	   FITOB_OUT_LEVEL3(verb(), " dimension adaptive = TRUE ");
	   FITOB_ERROR_TEST( dimensionAdaptiveLevels_.size() > 0 , "FitobCalculator dimensionAdaptiveLevels_.size():" << dimensionAdaptiveLevels_.size());
	   for (unsigned int tt = 0 ; tt < dimensionAdaptiveLevels_.size() ; tt++)
		   maxV = (dimensionAdaptiveLevels_[tt] > maxV) ? dimensionAdaptiveLevels_[tt] : maxV ;
	   globalLevel = maxV;
    }

   	// this means that we have a CONVERGENCE ANALYSIS with dimension adaptivity
   	if ( (convergenceRunNr >= 0) &&
		("true" == xmlconfiguration_->getStringConfiguration("thetaconfigurations.gridproperties.convergence-tool.<xmlattr>.dimAdaptive")) )
    {
            dimensionAdaptive_ = true;
        	xmlconfiguration_->getIntVectorConfiguration("thetaconfigurations.gridproperties.convergence-tool.CONV_DIM_ADAPT",
        	        convergenceRunNr , "<xmlattr>.level_vect" , ',' , dimensionAdaptiveLevels_);
        	//FITOB_OUT_LEVEL3(verb(), " dimensionAdaptiveLevels_.size() = " << dimensionAdaptiveLevels_.size() );
        	// set for the factory the adaptive truncation level, in case of T-CT
        	DVector adaptTruncation;
        	xmlconfiguration_->getDoubleVectorConfiguration("thetaconfigurations.gridproperties.convergence-tool.CONV_DIM_ADAPT",
        			convergenceRunNr , "<xmlattr>.trunc_level" , ',' , adaptTruncation);
        	//FITOB_OUT_LEVEL3(verb(), " adaptTruncation.size() = " << adaptTruncation.size() );
        	if (adaptTruncation.size() > 0)
        		gridFactory_->setDimensionAdaptiveTuncationLevel(adaptTruncation);

       		maxV = 0;
     	    for (unsigned int tt = 0 ; tt < dimensionAdaptiveLevels_.size() ; tt++)
     		   maxV = (dimensionAdaptiveLevels_[tt] > maxV) ? dimensionAdaptiveLevels_[tt] : maxV ;
     	    globalLevel = maxV;
     	   FITOB_OUT_LEVEL3(verb(), " dimension adaptive = TRUE , CONV ANALYSIS , globalLevel=" << globalLevel );
    }

   	IVector levels(scriptmodel_->getNrImportVariables(),globalLevel);
   	FITOB_OUT_LEVEL3(verb(),"FitobCalculator scriptmodel_->getNrImportVariables():" << scriptmodel_->getNrImportVariables());
   	FITOB_OUT_LEVEL3(verb(),"FitobCalculator scriptmodel_->getNrExportVariables():" << scriptmodel_->getNrExportVariables());
   	FITOB_OUT_LEVEL3(verb(),"FitobCalculator global Level:" << globalLevel);
   	// ==== END DETERMINING THE LEVEL OF THE COMPUTATIONS =============


   	DVector minValues(scriptmodel_->getNrImportVariables(), 0.0);
   	DVector maxValues(scriptmodel_->getNrImportVariables(), 0.0);
   	DVector minNormMeasure(scriptmodel_->getNrImportVariables(), 0.0);
   	DVector maxNormMeasure(scriptmodel_->getNrImportVariables(), 0.0);
   	DVector evalValues(scriptmodel_->getNrImportVariables(), 0.0);
	// - create initial domain
   	const int nrDomainAxisSpec =
   			xmlconfiguration_->nrXMLNodes("thetaconfigurations.domain");
   	for (int kk = 0 ; kk < nrDomainAxisSpec ; kk++ )
   	{
   		string varName = xmlconfiguration_->getStringConfiguration("thetaconfigurations.domain",kk,"<xmlattr>.name");
   		// get the import variable index
   		int axisIndex = scriptmodel_->getImportVariableIndex(varName);
   		// set the values in the arrays
   		minValues[axisIndex] = xmlconfiguration_->getDoubleConfiguration("thetaconfigurations.domain",kk,"<xmlattr>.min");
   		maxValues[axisIndex] = xmlconfiguration_->getDoubleConfiguration("thetaconfigurations.domain",kk,"<xmlattr>.max");
   		minNormMeasure[axisIndex] = xmlconfiguration_->getDoubleConfiguration("thetaconfigurations.domain",kk,"<xmlattr>.minNormDomain");
   		maxNormMeasure[axisIndex] = xmlconfiguration_->getDoubleConfiguration("thetaconfigurations.domain",kk,"<xmlattr>.maxNormDomain");
   		evalValues[axisIndex] = xmlconfiguration_->getDoubleConfiguration("thetaconfigurations.domain",kk,"<xmlattr>.evaluation");

   		// the norm measuring domain must be within the start domain, if not then treat this case
   		if ( (minValues[axisIndex] > minNormMeasure[axisIndex]) || (maxValues[axisIndex] < maxNormMeasure[axisIndex])){
   			if (convergenceRunNr >= 0) {
   				// limit the norm eval domain to the start domain
   				FITOB_WARNING(" Norm domain must be within start domain !!! ... automatically enforcing that condition minD:"
   					<< minNormMeasure[axisIndex] << " , maxD:" << maxNormMeasure[axisIndex]);
   			}
   			minNormMeasure[axisIndex] = minValues[axisIndex];
   			maxNormMeasure[axisIndex] = maxValues[axisIndex];
   		}
   	}

   	// create the evaluation domain and the start domain for the forward estimation
   	evalDomain_ = boost::shared_ptr<Domain>(new Domain( levels , evalValues , evalValues , scriptmodel_->getNrExportVariables() ));

   	// for higher dimensions we approximate with 3 points per dimension
   	if (evalDomain_->nrRealAxis() > 7) { FitobCalculator::norm_level_measure = FITOB_IMAX( 1 , 10-evalDomain_->nrRealAxis() );}
   	// create the domain on which the norm will be measured
   	IVector norm_levels(levels.size(),FitobCalculator::norm_level_measure);
   	normMeasureDomain_ = boost::shared_ptr<Domain>(new Domain( norm_levels , minNormMeasure , maxNormMeasure , scriptmodel_->getNrExportVariables() ));

   	// calculate the evaluation point
   	DVector evalPoint(0);
   	for (int i = 0 ; i < (int)minValues.size() ; i++ ){
   		if ( fabs(minValues[i]-maxValues[i]) > 1e-10 ){
   			evalPoint.push_back(evalValues[i]);
   		}
   	}
   	startDomain_ = boost::shared_ptr<Domain>(new Domain( levels , minValues , maxValues , scriptmodel_->getNrExportVariables() , &(evalPoint) ));

   	//FITOB_OUT_LEVEL3(verb(), "FitobCalculator evalDomain:" << evalDomain_->toString());
   	//FITOB_OUT_LEVEL3(verb(), "FitobCalculator startDomain:" << startDomain_->toString());
}

void FitobCalculator::forwardEstimation(){
	// here we do the forward estimation
	double timeStamp = 0.0;
	// set the stack index to zero
	stackIndex_ = 0;
	// put the first context into the stack
	Domain startDom(startDomain_.get());
	contextStack_.push_back(new MeshContext(startDom,timeStamp));
	// do the forward estimation
	scriptmodel_->bodyOpSeq().forwardEstimation( contextStack_ , this , stackIndex_ , timeStamp );
	// store end time
	endTime_ = timeStamp;
}

void FitobCalculator::backwardCalculation(){
	// do here the backward calculation
	double timeStamp = endTime_;
	// get the correct stack size
	stackIndex_ = contextStack_.size()-1;
	// generate mesh at the first level,
	FITOB_OUT_LEVEL2(verb(),"FitobCalculator::backwardCalculation create initial Grid" );
	//FITOB_OUT_LEVEL4(verb(),"FitobCalculator::backwardCalculation Domain:" << contextStack_[stackIndex_].getDom().toString() );
	boost::shared_ptr<fitob::MeshBase> firstGrid = gridFactory_->createMesh( &(contextStack_[stackIndex_].getDom()) , this );
	FITOB_OUT_LEVEL2(verb(),"FitobCalculator::backwardCalculation set Grid " );
	contextStack_[stackIndex_].setMesh( firstGrid );
	// do the backward calculation , call operation sequence
	FITOB_OUT_LEVEL2(verb(),"FitobCalculator::backwardCalculation call op sequence backward calculation " );
	scriptmodel_->bodyOpSeq().backwardCalculation( contextStack_ , this , stackIndex_ , timeStamp );
}

double FitobCalculator::evaluateScript(){

	std::clock_t _start_time = std::clock();

	time_t rawtime;
	struct tm * timeinfo;
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	string start_time_string = asctime (timeinfo);

	// do the forward estimation
	forwardEstimation();
	// backward calculation
	backwardCalculation();
#if defined(FITOB_MPI)
	// evaluate the grid on the evaluation point
	double res = 0.0;
	if ( FITOB_MPI_Comm_rank() == 0 ) {
		res = contextStack_[0].eval(evalDomain_.get()->getAverage());
	}
	FITOB_MPI_Bcast_Double( &res , 1 , 0 );
#else
	double res = contextStack_[0].eval(evalDomain_.get()->getAverage());
#endif

	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	string end_time_string = asctime (timeinfo);

	FITOB_OUT_LEVEL1(verb(),"FitobCalculator::evaluateScript  eval = " << res <<
			" , Runtime(sec):" << (double(std::clock() - _start_time) / CLOCKS_PER_SEC) );
	FITOB_OUT_LEVEL1(verb()," Start date: " << start_time_string );
	FITOB_OUT_LEVEL1(verb()," End date: " << end_time_string);
	return res;
}

boost::shared_ptr<FullGrid> FitobCalculator::evaluateScript_FG(double& retVal, int level_in)
{
	std::clock_t _start_time= std::clock();

	time_t rawtime;
	struct tm * timeinfo;
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	string start_time_string = asctime (timeinfo);

	// do the forward estimation
	forwardEstimation();
	FITOB_OUT_LEVEL2(verb(),"FitobCalculator::evaluateScript_FG  start backward calc " );
	// backward calculation
	backwardCalculation();
	FITOB_OUT_LEVEL2(verb(),"FitobCalculator::evaluateScript_FG  end backward calc " );

#if defined(FITOB_MPI)
	// evaluate the grid on the evaluation point
	if ( FITOB_MPI_Comm_rank() == 0 ) {
		retVal = contextStack_[0].eval(evalDomain_.get()->getAverage());
	}
	FITOB_MPI_Bcast_Double( &retVal , 1 , 0 );
#else
	// evaluate the grid on the evaluation point
	retVal = contextStack_[0].eval(evalDomain_.get()->getAverage());
#endif

	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	string end_time_string = asctime (timeinfo);

	FITOB_OUT_LEVEL1(verb(),"FitobCalculator::evaluateScript_FG  eval = " << retVal <<
			" , Runtime(sec):" << (double(std::clock() - _start_time) / CLOCKS_PER_SEC) );
	FITOB_OUT_LEVEL1(verb()," Start date: " << start_time_string );
	FITOB_OUT_LEVEL1(verb()," End date: " << end_time_string);

	// import variables
	int nrImp = normMeasureDomain_->nrImportVariables();
	// todo: configurable resolution for the full grid, which measures the convergence of the method
	// or make it depending from the dimensions
	int level_ref = level_in; //norm_level_measure;
	IVector lev( nrImp , level_ref);
	DVector minVal( nrImp );
	DVector maxVal( nrImp );

	for (unsigned int i = 0 ; i < (unsigned int)nrImp ; i++){
		minVal[i] = normMeasureDomain_->getGradedAxisMin(i+scriptmodel_->getNrExportVariables()+1); // because of global variable indexing
		maxVal[i] = normMeasureDomain_->getGradedAxisMax(i+scriptmodel_->getNrExportVariables()+1);
		FITOB_OUT_LEVEL3(verb() , " Axis i:" << i << " , min:" << minVal[i]<< " , max:" << maxVal[i] );
	}

	//
	Domain dom_tmp = Domain( lev , minVal , maxVal );

	boost::shared_ptr<FullGrid> ret = boost::shared_ptr<FullGrid>(new FullGrid(&dom_tmp));

	// in case of MPI only at rank one
#if defined(FITOB_MPI)
	// evaluate the grid on the evaluation point
	if ( FITOB_MPI_Comm_rank() == 0 ) {
		ret->setValues( &(contextStack_[0]) , evalDomain_.get()->getAverage() );
	}
	// eventually later this grid could be redistributed
	//FITOB_MPI_Bcast_Double( &retVal , 1 , 0 );
#else
	ret->setValues( &(contextStack_[0]) , evalDomain_.get()->getAverage() );
#endif

	return ret;
}

DVector FitobCalculator::evaluateScript_MonteCarlo(){

	// the result vector
	DVector res( 0 , -1.0 );

	FITOB_OUT_LEVEL2(verb(),"FitobCalculator::evaluateScript_MonteCarlo BEFORE forwardEstimation()");
	// do the forward estimation
	forwardEstimation();
	FITOB_OUT_LEVEL2(verb(),"FitobCalculator::evaluateScript_MonteCarlo AFTER forwardEstimation()");

	// do the forward simulation of the MonteCarlo paths
	int stackIndex = 0 , i , evaluationMode = -1;
	double timeStamp = 0.0;
	mCMachine_ = boost::shared_ptr<MCMachine>( new MCMachine( xmlconfiguration_.get() , modelcollection_ ) );
	Domain domStart = Domain(evalDomain_.get());
	// since the eval domain is just a point we extend the domain, so that the diffusion axis(risk factor) will be slightly recognized
	// it is important that the Monte-Carlo will have the same starting point as it is specified in the XML file
	DVector avrgV = domStart.getAverage();
	for (int a = 0 ; a < startDomain_->nrRealAxis() ; a++){
		// create a new domain, which is still variable but only 1e-10
		DVector vect = startDomain_->getGradedAxis( startDomain_->localToGlobalIndex(a) );
		for (int i = 0 ; i < (int)vect.size() ; i++){
			vect[i] = avrgV[ startDomain_->localToGlobalIndex(a) ] + 1e-10*( vect[i] - avrgV[ startDomain_->localToGlobalIndex(a) ] );
		}
		// set the axis which is just slightly changes
		domStart.setAxis( vect , startDomain_->localToGlobalIndex(a) );
	}

	mCMachine_->addMCStep( domStart , 0.0 , 0 );
	stackIndex = 1; timeStamp = contextStack_[stackIndex].getTime();
	FITOB_OUT_LEVEL2(verb(),"FitobCalculator::evaluateScript_MonteCarlo BEFORE forwardMCSimulation()");
	scriptmodel_->bodyOpSeq().forwardMCSimulation( contextStack_ , this , mCMachine_.get() , stackIndex , timeStamp );
	FITOB_OUT_LEVEL2(verb(),"FitobCalculator::evaluateScript_MonteCarlo AFTER forwardMCSimulation()");

	// ---
	for ( i = 0 ; i < mCMachine_->getNrEvaluation() ; i++ )
	{
		// get the actual evaluation mode
		evaluationMode = mCMachine_->getEvaluationMode(i);

		if ( evaluationMode > 0 ){
			// forward evaluation
			stackIndex = 0; timeStamp = 0.0;
			FITOB_OUT_LEVEL2(verb(),"FitobCalculator::evaluateScript_MonteCarlo BEFORE forwardMCSimulation()");
			scriptmodel_->bodyOpSeq().forwardMCEvaluation( contextStack_ , this , mCMachine_.get() , stackIndex , timeStamp );
			FITOB_OUT_LEVEL2(verb(),"FitobCalculator::evaluateScript_MonteCarlo AFTER forwardMCSimulation()");
		} else {
			// backward evaluation
			stackIndex = mCMachine_->nrMCSteps()-1; timeStamp = endTime_;
			FITOB_OUT_LEVEL2(verb(),"FitobCalculator::evaluateScript_MonteCarlo BEFORE backwardMCEvaluation()");
			scriptmodel_->bodyOpSeq().backwardMCEvaluation( contextStack_ , this , mCMachine_.get() , stackIndex , timeStamp );
			FITOB_OUT_LEVEL2(verb(),"FitobCalculator::evaluateScript_MonteCarlo AFTER backwardMCEvaluation()");
		}
	}

	//depending on the last evaluation mode calculate
	if ( evaluationMode < 0 ){
		// backward evaluation, take the first MCStep
		const Domain &dom = mCMachine_->getMCStep(0).getDom();
		res.resize(dom.nrExportVariables());
		for (i = 0 ; i < (int)res.size() ; i++){
			res[i] = mCMachine_->getMCStep(0).getAverage()[i+1];
			FITOB_OUT_LEVEL2(verb()," Export["<<i<<"] " << " Name:" << scriptmodel_->getExportVariableName(i) << " , res: "<< res[i]);
		}
	} else {
		// forward evaluation, take the last MCStep
		const Domain &dom = mCMachine_->getMCStep(mCMachine_->nrMCSteps()-1).getDom();
		res.resize(dom.nrExportVariables());
		for (i = 0 ; i < (int)res.size() ; i++){
			res[i] = mCMachine_->getMCStep(mCMachine_->nrMCSteps()-1).getAverage()[i+1];
			FITOB_OUT_LEVEL2(verb()," Export["<<i<<"] " << " Name:" << scriptmodel_->getExportVariableName(i) << " , res: "<< res[i]);
		}
	}
	FITOB_OUT_LEVEL2(verb()," FitobCalculator::evaluateScript_MonteCarlo END ");

	return res;
}

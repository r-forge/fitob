/*
 * FitobRegularization.cpp
 *
 *  Created on: Jul 7, 2011
 *      Author: benk
 */

#include "FitobRegularization.hpp"
#include "src/scripteval/FitobCalculator.hpp"
#include "src/mesh/FitobMeshBase.hpp"
#include "src/mesh/FitobFullGrid.hpp"
#include "src/mesh/FitobExtrapGrid.hpp"
#include "src/mesh/FitobCombiGrid.hpp"
#include "src/mesh/FitobSGppRegMesh.hpp"

// define just for the developement
//#define COMBI_REGRESSION_SOLVER


// include the combi grid library
#include "combigrid.hpp"
#ifdef COMBI_REGRESSION_SOLVER
#include "combigrid_solvers.hpp"
#include "sgpp_parallel.hpp"
#include "sgpp_solver.hpp"
#include "datadriven/application/LearnerBaseSP.hpp"
#endif

using namespace fitob;

Regularization::Regularization( const FitobCalculator* calc , DVector& XCoords , DVector& YCoords ,
		const Domain& dom , const ExpectedExpression* expectExpr) :
    actualMesh_() , fitobDomain_() , expectedExpression_(expectExpr) , domain_(0)
#ifdef COMBI_REGRESSION_SOLVER
    , myLearner_(0) , myLearnerSP_(0)
#endif
    {

	int verb = 6;

	// signal to the expression that this triggered one expression
	expectedExpression_->setITriggeredRegression(true);

	string solverType = calc->getXMLConfiguration()->getStringConfiguration("thetaconfigurations.montecarlo.REGRESSION.<xmlattr>.solver");
	// read in which solver to choose , if is SGpp solver than enter in the special function for this
    if (solverType == "SGppTikhonov"){
    	configureSGppSolver( calc , XCoords , YCoords , dom );
    }
    else
    {
	// otherwise just use the Tikhonov solver from the Combigrid.
    	// create a combigrid Domain
    	FITOB_OUT_LEVEL3( verb ," Regularization::Regularization , create Domain  dim=" << dom.nrRealAxis()  );
    	// - create a domain
    	fitobDomain_ = boost::shared_ptr<Domain>( new Domain(dom));
		double lambda = calc->getXMLConfiguration()->getDoubleConfiguration("thetaconfigurations.montecarlo.REGRESSION.<xmlattr>.lambda");
		string forcelambda_str = calc->getXMLConfiguration()->getStringConfiguration("thetaconfigurations.montecarlo.REGRESSION.<xmlattr>.forcelamdba");
		bool forcelamdba = (forcelambda_str == "true");

		//create the
		std::vector< std::vector<double> > stretching(dom.nrRealAxis());
		for (int d = 0 ; d < dom.nrRealAxis() ; d++){
			int globI = dom.localToGlobalIndex(d);
			stretching[d] = dom.getGradedAxis( globI );
		}
		// create
		domain_ = new combigrid::GridDomain( dom.nrRealAxis() , stretching );
		// cerate a mesh (potentially combi grid )
		actualMesh_ = calc->getGridFactory()->createMesh( fitobDomain_.get() , calc );
		// - test if it is a full grid , then create
		MeshBase* actMesh = actualMesh_.get();

		FullGrid* fg = dynamic_cast<FullGrid*>(actMesh);
		ExtrapGrid* eg = dynamic_cast<ExtrapGrid*>(actMesh);
		CombiGrid* cg = dynamic_cast<CombiGrid*>(actMesh);

		if (fg != 0){
			FITOB_OUT_LEVEL3( verb ," Regularization::Regularization , Tikhonov on FullGrid " );
			combigrid::FullGridD* combiFg = new combigrid::FullGridD( dom.nrRealAxis() , dom.getLevelVector() );
			combiFg->createFullGrid();
			combiFg->setDomain( domain_ );
			// call the Tikhonov regulariziation solver, which also finds the optimal lambda parameter as well
#ifdef COMBI_REGRESSION_SOLVER
			if (forcelamdba){
				combigrid::RunTikhonov::computeFGTikhonov_FG( combiFg , combiFg->getElementVector() , lambda, XCoords, YCoords );
			} else {
				combigrid::RunTikhonov::computeTikhonov_FG_crossvalidation( combiFg , lambda, XCoords, YCoords);
			}
#endif
			// write the solution back
			fg->unknVect() = combiFg->getElementVector();
			combiFg->setDomain( 0 );
			delete combiFg;
		}else{
			if (eg != 0){
				FITOB_OUT_LEVEL3( verb ," Regularization::Regularization , Tikhonov on ExtrapGrid (CombiGrid) " );
				combigrid::AbstractCombiGrid* combigrid = eg->getAbstractCombiGrid();
				combigrid->setDomainAllFG( domain_ );
				// call the Tikhonov regulariziation solver, which also finds the optimal lambda parameter as well
#ifdef COMBI_REGRESSION_SOLVER
				if (forcelamdba){
					combigrid::RunTikhonov::computeTikhonov_CG( combigrid, lambda , XCoords , YCoords);
				} else {
					combigrid::RunTikhonov::computeTikhonov_CG_crossvalidation( combigrid , lambda , XCoords , YCoords);
				}
#endif
				combigrid->setDomainAllFG( 0 );
				FITOB_OUT_LEVEL3( verb ," Regularization::Regularization , END Tikhonov on ExtrapGrid (CombiGrid) " );
			}
			else
				if (cg != 0){
					FITOB_OUT_LEVEL3( verb ," Regularization::Regularization , Tikhonov on CombiGrid " );
					combigrid::AbstractCombiGrid* combigrid = cg->getAbstractCombiGrid();
					combigrid->setDomainAllFG( domain_ );
					// call the Tikhonov regulariziation solver, which also finds the optimal lambda parameter as well
#ifdef COMBI_REGRESSION_SOLVER
					if (forcelamdba){
						combigrid::RunTikhonov::computeTikhonov_CG( combigrid, lambda , XCoords , YCoords);
					} else {
						combigrid::RunTikhonov::computeTikhonov_CG_crossvalidation( combigrid , lambda , XCoords , YCoords);
					}
#endif
					combigrid->setDomainAllFG( 0 );
					FITOB_OUT_LEVEL3( verb ," Regularization::Regularization , END Tikhonov on CombiGrid" );
				}
				else{
					// throw error
					FITOB_ERROR_EXIT( " Regularization::Regularization , grid must be eighter FullGrid or CombiGrid(ExtrapGrid) with boundary points" );
				}
		}
    }
}

Regularization::~Regularization(){
	// tell the expression that the results are not used any more
	expectedExpression_->setITriggeredRegression(false);
	if (domain_ != 0) delete domain_;
#ifdef COMBI_REGRESSION_SOLVER
	if (myLearner_ != 0) { delete myLearner_; }
	if (myLearnerSP_ != 0) { delete myLearnerSP_; }
#endif
}


void Regularization::configureSGppSolver( const FitobCalculator* calc ,
		DVector& XCoords , DVector& YCoords , const Domain& dom ){
#ifdef COMBI_REGRESSION_SOLVER
	int verb = 6;
	// variables for the configuration
	int tmpINput = -1;
	double tmpINputDouble = -1.0;
	string tmpInStr = "";

	// create a combigrid Domain
	FITOB_OUT_LEVEL3( verb ," Regularization::Regularization , create Domain  dim=" << dom.nrRealAxis()  );

	// ------------- read in the configurations --------------------
	double lambda = calc->getXMLConfiguration()->getDoubleConfiguration("thetaconfigurations.montecarlo.REGRESSION.<xmlattr>.lambda");
	double domainFactor = calc->getXMLConfiguration()->getDoubleConfiguration("thetaconfigurations.montecarlo.REGRESSION.SGppLearner.domainFactor.<xmlattr>.value");

	sg::base::AdpativityConfiguration adaptConfig;
	sg::solver::SLESolverConfiguration SLESolverConfigRefine , SLESolverConfigFinal;
	sg::solver::SLESolverSPConfiguration SLESolverSPConfigRefine , SLESolverSPConfigFinal;
	sg::parallel::VectorizationType vecType;
    sg::base::RegularGridConfiguration GridConfig;

    // choose the vectorization type
	string vectorization = calc->getXMLConfiguration()->getStringConfiguration(
			 "thetaconfigurations.montecarlo.REGRESSION.SGppLearner.vectType.<xmlattr>.value");
	if (vectorization == "X86SIMD") {
		vecType = sg::parallel::X86SIMD;
	} else if  (vectorization == "OCL") {
		vecType = sg::parallel::OpenCL;
	} else if  (vectorization == "HYBRID_X86SIMD_OCL") {
		vecType = sg::parallel::Hybrid_X86SIMD_OpenCL;
	} else if  (vectorization == "ArBB") {
		vecType = sg::parallel::ArBB;
	} else {
		vecType = sg::parallel::X86SIMD;
	}

	// configuration for the adaptive
	tmpInStr = calc->getXMLConfiguration()->getStringConfiguration("thetaconfigurations.montecarlo.REGRESSION.SGppLearner.adaptConfigmaxLevelType.<xmlattr>.value");
	if (tmpInStr == "true")
		 { adaptConfig.maxLevelType_ = true; }
	else { adaptConfig.maxLevelType_ = false; }
	tmpINput = calc->getXMLConfiguration()->getIntConfiguration("thetaconfigurations.montecarlo.REGRESSION.SGppLearner.adaptConfignoPoints.<xmlattr>.value");
	if (tmpINput != -1 )
		 { adaptConfig.noPoints_ = static_cast<size_t>(tmpINput);}
	else { adaptConfig.noPoints_ = 0; }
	tmpINput = calc->getXMLConfiguration()->getIntConfiguration("thetaconfigurations.montecarlo.REGRESSION.SGppLearner.adaptConfignumRefinements.<xmlattr>.value");
	if (tmpINput != -1 )
		 { adaptConfig.numRefinements_ = static_cast<size_t>(tmpINput);}
	else { adaptConfig.numRefinements_ = 0; }
	tmpINputDouble = calc->getXMLConfiguration()->getDoubleConfiguration("thetaconfigurations.montecarlo.REGRESSION.SGppLearner.adaptConfigpercent.<xmlattr>.value");
	//if (tmpINputDouble != -1 ) { adaptConfig.percent_ = tmpINputDouble;} else { adaptConfig.percent_ = 100.0; }
	adaptConfig.percent_ = tmpINputDouble;
	tmpINputDouble = calc->getXMLConfiguration()->getDoubleConfiguration("thetaconfigurations.montecarlo.REGRESSION.SGppLearner.adaptConfigthreshold.<xmlattr>.value");
	//if (tmpINputDouble != -1 ) { adaptConfig.threshold_ = tmpINputDouble;} else { adaptConfig.threshold_ = 0.0; }
	adaptConfig.threshold_ = tmpINputDouble;

	// Set solver for refinement
	tmpINputDouble = calc->getXMLConfiguration()->getDoubleConfiguration(
			"thetaconfigurations.montecarlo.REGRESSION.SGppLearner.solverConfigRefine_eps.<xmlattr>.value");
	if (tmpINputDouble != -1 )
		 { SLESolverConfigRefine.eps_ = tmpINputDouble;
		   SLESolverSPConfigRefine.eps_ = static_cast<float>(tmpINputDouble); }
	else { SLESolverConfigRefine.eps_ = 1e-5; SLESolverSPConfigRefine.eps_ = 1e-5;}
	tmpINput = calc->getXMLConfiguration()->getIntConfiguration(
			"thetaconfigurations.montecarlo.REGRESSION.SGppLearner.solverConfigRefine_maxIterations.<xmlattr>.value");
	if (tmpINput != -1 )
	     { SLESolverSPConfigRefine.maxIterations_ = SLESolverConfigRefine.maxIterations_ = static_cast<size_t>(tmpINput);}
	else { SLESolverSPConfigRefine.maxIterations_ = SLESolverConfigRefine.maxIterations_ = 10; }
	tmpINputDouble = calc->getXMLConfiguration()->getDoubleConfiguration(
			"thetaconfigurations.montecarlo.REGRESSION.SGppLearner.solverConfigRefine_threshold.<xmlattr>.value");
	//if (tmpINputDouble != -1 ) { SLESolverConfigRefine.threshold_ = tmpINputDouble;} else { SLESolverConfigRefine.threshold_ = 1.0; }
	SLESolverConfigRefine.threshold_ = tmpINputDouble;
	SLESolverSPConfigRefine.threshold_ = static_cast<float>(tmpINputDouble);
	string solverTypeRefine = calc->getXMLConfiguration()->getStringConfiguration(
			 "thetaconfigurations.montecarlo.REGRESSION.SGppLearner.solverConfigRefine_type.<xmlattr>.value");
	if (solverTypeRefine == "BiCGSTAB")
		 { SLESolverSPConfigRefine.type_ = SLESolverConfigRefine.type_ = sg::solver::BiCGSTAB; }
	else { SLESolverSPConfigRefine.type_ = SLESolverConfigRefine.type_ = sg::solver::CG; }


	// Set solver for final step
	tmpINputDouble = calc->getXMLConfiguration()->getDoubleConfiguration(
			"thetaconfigurations.montecarlo.REGRESSION.SGppLearner.solverConfigFinal_eps.<xmlattr>.value");
	if (tmpINputDouble != -1 )
		 { SLESolverConfigFinal.eps_ = tmpINputDouble;
		   SLESolverSPConfigFinal.eps_ = static_cast<float>(tmpINputDouble);}
	else { SLESolverConfigFinal.eps_ = 1e-5; SLESolverConfigFinal.eps_ = 1e-5;}
	tmpINput = calc->getXMLConfiguration()->getIntConfiguration(
			"thetaconfigurations.montecarlo.REGRESSION.SGppLearner.solverConfigFinal_maxIterations.<xmlattr>.value");
	if (tmpINput != -1 )
		 { SLESolverSPConfigFinal.maxIterations_ = SLESolverConfigFinal.maxIterations_ = static_cast<size_t>(tmpINput);}
	else { SLESolverSPConfigFinal.maxIterations_ = SLESolverConfigFinal.maxIterations_ = 10; }
	tmpINputDouble = calc->getXMLConfiguration()->getDoubleConfiguration(
			"thetaconfigurations.montecarlo.REGRESSION.SGppLearner.solverConfigFinal_threshold.<xmlattr>.value");
	//if (tmpINputDouble != -1 ) { SLESolverConfigFinal.threshold_ = tmpINputDouble;} else { SLESolverConfigFinal.threshold_ = 1.0; }
	SLESolverConfigFinal.threshold_ = tmpINputDouble;
	SLESolverSPConfigFinal.threshold_ = static_cast<float>(tmpINputDouble);
	string solverTypeFinal = calc->getXMLConfiguration()->getStringConfiguration(
			 "thetaconfigurations.montecarlo.REGRESSION.SGppLearner.solverConfigFinal_type.<xmlattr>.value");
	if (solverTypeFinal == "BiCGSTAB")
	 	 { SLESolverSPConfigFinal.type_ = SLESolverConfigFinal.type_ = sg::solver::BiCGSTAB; }
	else { SLESolverSPConfigFinal.type_ = SLESolverConfigFinal.type_ = sg::solver::CG; }

	// ---- grid configuration ------
	GridConfig.dim_ = dom.nrRealAxis();
	string gridtype = calc->getXMLConfiguration()->getStringConfiguration(
			 "thetaconfigurations.montecarlo.REGRESSION.SGppLearner.gridType.<xmlattr>.value");
	if (gridtype == "linearboundary") {
		GridConfig.type_ = sg::base::LinearTrapezoidBoundary;
	} else if (gridtype == "modlinear") {
		GridConfig.type_ = sg::base::ModLinear;
	} else if (gridtype == "linear") {
		GridConfig.type_ = sg::base::Linear;
	} else {
		GridConfig.type_ = sg::base::LinearTrapezoidBoundary;
	}
	tmpINput = calc->getXMLConfiguration()->getIntConfiguration("thetaconfigurations.montecarlo.REGRESSION.SGppLearner.gridConfiglevel.<xmlattr>.value");
	if (tmpINput != -1 ) { GridConfig.level_ = static_cast<size_t>(tmpINput);} else { GridConfig.level_ = 3;}

	string precisionStr = calc->getXMLConfiguration()->getStringConfiguration(
			 "thetaconfigurations.montecarlo.REGRESSION.SGppLearner.precision.<xmlattr>.value");

    int nrMCPoints = YCoords.size();
    int dim = dom.nrRealAxis();

    size_t i;
    DVector minValues(dim) , maxValues(dim);
    for (i = 0; i < (size_t)dim ; i++) { minValues[i] = +1e+100; maxValues[i] = -1e+100; }
    // get the minimal and maximal value per axis
    for (size_t mcInt = 0 ; mcInt < (size_t)nrMCPoints ; mcInt++){
        for (i = 0; i < (size_t)dim ; i++){
        	minValues[i] = (minValues[i] > XCoords[ mcInt*dim + i ]) ? XCoords[ mcInt*dim + i ] : minValues[i];
        	maxValues[i] = (maxValues[i] < XCoords[ mcInt*dim + i ]) ? XCoords[ mcInt*dim + i ] : maxValues[i];
        }
    }

    sg::base::DataMatrix data((size_t)1 , (size_t)dim);
    sg::base::DataVector classes((size_t)1);
    sg::base::DataMatrixSP dataSP((size_t)1 , (size_t)dim);
    sg::base::DataVectorSP classesSP((size_t)1);


    if ( precisionStr == "SP" ){
    	// create the learner
        myLearnerSP_ = new sg::parallel::LearnerVectorizedIdentitySP(vecType, true , true);
    	classesSP.resize((size_t)nrMCPoints);
    	dataSP.resize((size_t)nrMCPoints);
        // copy the data to the DataVectors and DataMatrix
        for (size_t mcInt = 0 ; mcInt < (size_t)nrMCPoints ; mcInt++){
        	classesSP[mcInt] = YCoords[mcInt];
        	for (i = 0; i < (size_t)dim ; i++){
        	    // scaling and copy that data to "data" and "classes"
        		// todo: evtl. use "domainFactor" configuration
        		dataSP.set(mcInt , i , (XCoords[ mcInt*dim + i ]- minValues[i])/(maxValues[i] - minValues[i])  );
        	}
        }
    }
    else
    { // by default we have "DP"
    	// create the learner
        myLearner_ = new sg::parallel::LearnerVectorizedIdentity(vecType, true , true);
    	classes.resize((size_t)nrMCPoints);
    	data.resize((size_t)nrMCPoints);
        // copy the data to the DataVectors and DataMatrix
        for (size_t mcInt = 0 ; mcInt < (size_t)nrMCPoints ; mcInt++){
        	classes[mcInt] = YCoords[mcInt];
        	for (i = 0; i < (size_t)dim ; i++){
        	    // scaling and copy that data to "data" and "classes"
        		// todo: evtl. use "domainFactor" configuration
        		data.set(mcInt , i , (XCoords[ mcInt*dim + i ]- minValues[i])/(maxValues[i] - minValues[i])  );
        	}
        }
    }

    // define the learner class
    FITOB_OUT_LEVEL3( verb ," Regularization::Regularization , START SGppLearner .... " );
    FITOB_OUT_LEVEL3( verb ,"lambda = "  << lambda );
    FITOB_OUT_LEVEL3( verb ,"vectType value = "  << vecType << " , str = " << vectorization);
    FITOB_OUT_LEVEL3( verb ,"gridConfiglevel value = " << GridConfig.level_);
    FITOB_OUT_LEVEL3( verb ,"gridType value = " << GridConfig.type_ << " , str = " << gridtype );
    FITOB_OUT_LEVEL3( verb ,"domainFactor = " << domainFactor );
    FITOB_OUT_LEVEL3( verb ,"adaptConfigmaxLevelType   value = " << adaptConfig.maxLevelType_ );
    FITOB_OUT_LEVEL3( verb ,"adaptConfignoPoints   value = " << adaptConfig.noPoints_ );
    FITOB_OUT_LEVEL3( verb ,"adaptConfignumRefinements   value = " << adaptConfig.numRefinements_ );
    FITOB_OUT_LEVEL3( verb ,"adaptConfigpercent   value = " << adaptConfig.percent_ );
    FITOB_OUT_LEVEL3( verb ,"adaptConfigthreshold   value = " << adaptConfig.threshold_);
    FITOB_OUT_LEVEL3( verb ,"solverConfigRefine_eps   value = " << SLESolverConfigRefine.eps_);
    FITOB_OUT_LEVEL3( verb ,"solverConfigRefine_maxIterations   value = " << SLESolverConfigRefine.maxIterations_);
    FITOB_OUT_LEVEL3( verb ,"solverConfigRefine_threshold   value = " << SLESolverConfigRefine.threshold_);
    FITOB_OUT_LEVEL3( verb ,"solverConfigRefine_type   value = " << SLESolverConfigRefine.type_ );
    FITOB_OUT_LEVEL3( verb ,"solverConfigFinal_eps   value = " << SLESolverConfigFinal.eps_ );
    FITOB_OUT_LEVEL3( verb ,"solverConfigFinal_maxIterations   value = " << SLESolverConfigFinal.maxIterations_ );
    FITOB_OUT_LEVEL3( verb ,"solverConfigFinal_threshold   value = " << SLESolverConfigFinal.threshold_ );
    FITOB_OUT_LEVEL3( verb ,"solverConfigFinal_type   value = " << SLESolverConfigFinal.type_ );
    FITOB_OUT_LEVEL3( verb ,"precisionStr   value = " << precisionStr );
    sg::datadriven::LearnerTiming gtimings;

    // make the regularization with SGpp
    if ( precisionStr == "SP" ){
    	gtimings = myLearnerSP_->train( dataSP , classesSP , GridConfig,
    			SLESolverSPConfigRefine, SLESolverSPConfigFinal, adaptConfig, false, static_cast<float>(lambda));
    	/*gtimings = myLearner->train(    dataSP , classesSP , GridConfig, SolverConfigRefine   ,	SolverConfigFinal   , AdaptConfig, false, lambda);
    	train(sg::base::DataMatrixSP& testDataset,
    		  sg::base::DataVectorSP& classes,
    		  const sg::base::RegularGridConfiguration& GridConfig,
    		  const sg::solver::SLESolverSPConfiguration& SolverConfigRefine,
    		  const sg::solver::SLESolverSPConfiguration& SolverConfigFinal,
    		  const sg::base::AdpativityConfiguration& AdaptConfig,
    		  bool testAccDuringAdapt, const float lambda);*/
    } else {
    	gtimings = myLearner_->train(data, classes, GridConfig, SLESolverConfigRefine, SLESolverConfigFinal, adaptConfig, false, lambda);
    }

    FITOB_OUT_LEVEL3( verb ," Regularization::Regularization , END SGppLearner .... " );
    //FITOB_OUT_LEVEL3( verb ," DATA = " << data.toString() );
    //FITOB_OUT_LEVEL3( verb ," VECT = " << classes.toString() );

	// set the domain accordingly, level is a dummy 3, attention to the
    IVector axisLevel( dom.nrImportVariables() , 3);
    DVector minValuesD( dom.nrImportVariables() , 0.0) , maxValuesD( dom.nrImportVariables() , 0.0 );
    for (i = 0; i < (size_t)dom.nrImportVariables() ; i++){
    	minValuesD[i] = maxValuesD[i] = dom.getAverage()[dom.importToGlobal(i)];
    }
    for (i = 0; i < (size_t)dim ; i++){
    	// compute the index in the import variables
    	int Ind = dom.globalToImport( dom.localToGlobalIndex(i) );
    	minValuesD[Ind] = minValues[i];
    	maxValuesD[Ind] = maxValues[i];
    	//FITOB_OUT_LEVEL3( verb ," minValuesD[" << Ind << "]=" << minValuesD[Ind]);
    	//FITOB_OUT_LEVEL3( verb ," maxValuesD[" << Ind << "]=" << maxValuesD[Ind]);
    }

    // define the domain of the mesh
    //FITOB_OUT_LEVEL3( verb ," Regularization::Regularization , dom.nrGlobalVariables()=" << dom.nrGlobalVariables()
    //		<< " , dom.nrImportVariables()=" << dom.nrImportVariables());
    //FITOB_OUT_LEVEL3( verb ," Regularization::Regularization , DOM_ORIG = " << dom.toString());
	fitobDomain_ = boost::shared_ptr<Domain>( new Domain(axisLevel , minValuesD , maxValuesD , dom.nrExportVariables() , NULL) );
	//FITOB_OUT_LEVEL3( verb ," Regularization::Regularization , DOM_NEW = " << fitobDomain_->toString());
    //FITOB_OUT_LEVEL3( verb ," Regularization::Regularization , fitobDomain_->nrGlobalVariables()=" << fitobDomain_->nrGlobalVariables()
    //		<< " , fitobDomain_->nrImportVariables()=" << fitobDomain_->nrImportVariables());

	// set the mesh so that
	actualMesh_ = boost::shared_ptr<MeshBase>(new SGppRegMesh(fitobDomain_.get() , myLearner_ , myLearnerSP_) );
#endif
}

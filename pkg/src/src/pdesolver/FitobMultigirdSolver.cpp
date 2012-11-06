/*
 * FitobMultigirdSolver.cpp
 *
 *  Created on: Jul 6, 2010
 *      Author: benk
 */

#include "FitobMultigirdSolver.hpp"
#include "src/scripteval/FitobCalculator.hpp"
#include "src/pdesolver/FitobMultigridFGBase.hpp"
#include "src/pdesolver/FitobMultigridFG.hpp"
#include "src/pdesolver/FitobMultigridFG_WB.hpp"
#include "src/pdesolver/FitobSmootherBase.hpp"
#include "src/pdesolver/FitobGaussSeidel.hpp"
#include "src/pdesolver/FitobGaussSeidel_WB.hpp"
#include "src/mesh/FitobSparseGrid.hpp"
#include "src/mesh/FitobFullGrid.hpp"
#include "src/mesh/FitobFullGrid_WB.hpp"
#include "src/mesh/FitobFullGridBase.hpp"

#if defined(FITOB_MPI)
#include <pthread.h>
#endif

using namespace fitob;
using namespace std;

MultigirdSolver::MultigirdSolver(const XMLConfiguration* config) : SolverBase(config) {

	setVerb(4);

	minMicroDT_ = config->getDoubleConfiguration("thetaconfigurations.solver.multigrid-solver.mindt.<xmlattr>.value" );
    maxMicroDT_ = config->getDoubleConfiguration("thetaconfigurations.solver.multigrid-solver.maxdt.<xmlattr>.value" );
	solverEpsilon_ = config->getDoubleConfiguration("thetaconfigurations.solver.multigrid-solver.solverEps.<xmlattr>.value" );

    if (config->getStringConfiguration("thetaconfigurations.solver.multigrid-solver.usepredictor.<xmlattr>.value") == "true")
   	  usePredictor_ = true;
    else usePredictor_ = false;

    if (config->getStringConfiguration("thetaconfigurations.solver.multigrid-solver.useTimeStepControl.<xmlattr>.value") == "true")
   	  makeTimeStepControll_ = true;
    else makeTimeStepControll_ = false;

    if (config->getStringConfiguration("thetaconfigurations.solver.multigrid-solver.timeStepControlInfNorm.<xmlattr>.value") == "true")
   	  useInfNormForTimeStepControll_ = true;
    else useInfNormForTimeStepControll_ = false;

	timeStepControll_Epsilon_ = config->getDoubleConfiguration("thetaconfigurations.solver.multigrid-solver.timeStepControl.<xmlattr>.value" );
	predictorUnderrelaxCoef_ = config->getDoubleConfiguration("thetaconfigurations.solver.multigrid-solver.predictorUnderRelaxCoef.<xmlattr>.value" );

	// read solver specific configuration parameters
	maxCombiDT_ = config->getDoubleConfiguration("thetaconfigurations.solver.multigrid-solver.maxCombiStep.<xmlattr>.value" );
	minCombiDT_ = config->getDoubleConfiguration("thetaconfigurations.solver.multigrid-solver.minCombiStep.<xmlattr>.value" );
	combiStepEps_ = config->getDoubleConfiguration("thetaconfigurations.solver.multigrid-solver.CombiStepControl.<xmlattr>.value" );

	FITOB_OUT_LEVEL3(verb(), " MultigirdSolver::solvePDE , configurations : minMicroDT" << minMicroDT_ << ", maxMicroDT:" << maxMicroDT_ <<
			", solverEpsilon:" << solverEpsilon_ << " usePredictor:" << usePredictor_ <<
			", makeTimeStepControll:" << makeTimeStepControll_ << " , useInfNormForTimeStepControll:" << useInfNormForTimeStepControll_ <<
			", timeStepControll_Epsilon:" << timeStepControll_Epsilon_ << ", predictorUnderrelaxCoef:" << predictorUnderrelaxCoef_ <<
			", maxCombiDT:" << maxCombiDT_)

}

#if defined(FITOB_MPI)
void* MultigirdSolver::MPI_distribute_Load(void *obj){
	// test if we got the right parameter
	SparseGrid* sgrid = (SparseGrid*)obj;
	FITOB_ERROR_TEST( sgrid != NULL , "MultigirdSolver::MPI_distribute_Load cast to SparseGrid failed");

	int nrGridAssigned = 0 , nrGridFinished = 0 , verbMy = 2;
	int nrProc = FITOB_MPI_Comm_size() , nrFullGrids = sgrid->nrFullGrid();
	FITOB_OUT_LEVEL3( verbMy , "  MultigirdSolver::solvePDE  Managing thread started ....");
	IVector sortedFullGrids_stat(nrFullGrids);

	// measure time
	time_t rawtime;
	struct tm * timeinfo;
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	string start_time_string = asctime (timeinfo);

	// recieve local array TODO
	int rankNr = 0 , fgIndex , tmp = -1;
	FITOB_MPI_Recv_Int( &(sortedFullGrids_stat[0]) , nrFullGrids , rankNr );

	// assign the first N grids
	for (int r = 0; r < nrProc ; r++ ){
		// if we exceed the number of full grids then just exit

		if ( r >= nrFullGrids ) { break; }
		// send the index of the full grid
		FITOB_MPI_Send_Int( &sortedFullGrids_stat[r] , 1 , r );
		// if this is not rank 0 then send the vector
		if ( r > 0) {
			FITOB_OUT_LEVEL3( verbMy , "  send FG" << sortedFullGrids_stat[r] << " index:" << r << " to P" << r );
			FITOB_MPI_Send_Double( &(sgrid->getFG(sortedFullGrids_stat[r]).unknVect()[0]) , sgrid->getFG(sortedFullGrids_stat[r]).totalPoints() , r );
		}
		nrGridAssigned++;
	}
	rankNr = 0; fgIndex = 0;  tmp = -1;
	// wait for the FGs to finish
    for (int g = 0; g < nrFullGrids; g++){
    	// Receive the index of the full grid that is ready
    	FITOB_OUT_LEVEL3( verbMy , "  Waiting for finished FGs " );
    	FITOB_MPI_Recv_Int_AnySource(&fgIndex , 1 , rankNr );
    	FITOB_OUT_LEVEL3( verbMy , "  P:" << rankNr << "   FG index:" << fgIndex);
    	if ( rankNr > 0 ){
    		FITOB_MPI_Recv_Double( &(sgrid->getFG(fgIndex).unknVect()[0]) , sgrid->getFG(fgIndex).totalPoints() , rankNr );
		}
    	// if we have unfinished full grids, then assigned the FG with index "nrGridAssigned"
    	if (nrGridAssigned < nrFullGrids){
    		FITOB_MPI_Send_Int( &sortedFullGrids_stat[nrGridAssigned] , 1 , rankNr );
    		FITOB_OUT_LEVEL3( verbMy , "  send FG:" << sortedFullGrids_stat[nrGridAssigned] << " with index:" << nrGridAssigned
    				<< " out of " << nrFullGrids  << " to P" << rankNr );
    		if ( rankNr > 0) {
    			FITOB_MPI_Send_Double( &(sgrid->getFG(sortedFullGrids_stat[nrGridAssigned]).unknVect()[0]) ,
    					sgrid->getFG(sortedFullGrids_stat[nrGridAssigned]).totalPoints() , rankNr );
    		}
        	nrGridAssigned++;
    	}
    	// nothing else to do then send -1 to signal that task done
    	else{
    		FITOB_OUT_LEVEL3( verbMy , "  send P:" << rankNr << "  no further FG to be done" );
    		FITOB_MPI_Send_Int( &tmp , 1 , rankNr );
    	}
    	nrGridFinished++;
    }

	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	string end_time_string = asctime (timeinfo);

	FITOB_OUT_LEVEL1(verbMy,"MPI PDE SOLVE  Start date: " << start_time_string );
	FITOB_OUT_LEVEL1(verbMy,"MPI PDE SOLVE  End date: " << end_time_string);
    // return
    return 0;
}
#endif

void MultigirdSolver::solvePDE(MeshBase *grid , const MeshContext* context ,
              const ModelCollection* models ,
              const FitobCalculator* calc ,
              double timeStep){

	// store the global coordinates
	globalCoordinates_backup = context->minGlobCoord();

	// --- old solver (Gauss-Seidl) call  -----
	//solvePDE_with_GS( dynamic_cast<FullGrid*>(grid) , context , models , calc , timeStep);

    // try to convert to sparse grid
	SparseGrid* sgrid = dynamic_cast<SparseGrid*>(grid);

	if ( sgrid != NULL)
	{
		// apply the maxCombi time step , make here a for loop
		double oneCombiTimeStep = 0.0 , nextStep, actualTime = 0.0;
		double combiErrorIndicator = 0.0;
		// make a simple for loop to make the calculated time step
		int co = 0;
		// take the minimum out of these two
		oneCombiTimeStep = (timeStep < minCombiDT_) ? timeStep : minCombiDT_;
		// if we have a simply extrapolation technique then we can make just one cycle
		if (sgrid->canCombineFullGrids() == false ) oneCombiTimeStep = timeStep;

		sortedFullGrids_.resize(sgrid->nrFullGrid());
		IVector fgSizes(sgrid->nrFullGrid());
		//  ----------- sort the array ------------------------
		for (int tt = 0 ; tt < sgrid->nrFullGrid() ; tt++)
		{
			fgSizes[tt] = sgrid->getFG(tt).totalPoints();  sortedFullGrids_[tt] = tt;
			//FITOB_OUT_LEVEL1(verb(), " subspace:" << tt << " with size:" << fgSizes[tt] );
		}
		// sort only if the nr of ful grids is not too high
		if ( sgrid->nrFullGrid() <= 10000 ){
			bool hasDiff = true; int tmpSwap = 0;
			// simple bubblesort
			while (hasDiff){
				hasDiff = false;
				for (int tt = 1 ; tt < sgrid->nrFullGrid() ; tt++){
					if ( fgSizes[tt-1] < fgSizes[tt]){
						hasDiff = true;
						tmpSwap = fgSizes[tt-1];  fgSizes[tt-1] = fgSizes[tt];  fgSizes[tt] = tmpSwap;
						tmpSwap = sortedFullGrids_[tt-1];  sortedFullGrids_[tt-1] = sortedFullGrids_[tt];  sortedFullGrids_[tt] = tmpSwap;
					}
				}
			}
		}
		for (int tt = 0 ; tt < sgrid->nrFullGrid() ; tt++)
		{
			//FITOB_OUT_LEVEL1(verb(), " sorted subspace:" << tt << " original pos:" << sortedFullGrids_[tt] << " with size:" << fgSizes[tt] );
		}

		// loop as long we have time left
		while ( (timeStep - actualTime) > 1e-15 )
		{
			// in case of sparse grid, decompose the sparse grid and do the solving for each subspace
			FITOB_OUT_LEVEL1(verb(), " MultigirdSolver::solvePDE SOLVING SparseGrid with Combi Technique cycle:" << co
					<< " , actualTime:" << actualTime << " , timestep:" << oneCombiTimeStep);
			// decompose the sparse grid
			FITOB_OUT_LEVEL1(verb(), " MultigirdSolver::solvePDE , deCompose() ");
			sgrid->deCompose();

			// ------------------- MPI parallel section ----------------
#if defined(FITOB_MPI)
 if ( FITOB_MPI_Comm_size() > 1 ){
			// Synchronize all processes
			FITOB_MPI_Barrier();

			// variables needed for the thread
			pthread_t thread_main;
			int  iret;

			if ( FITOB_MPI_Comm_rank() == 0 ) {
				// start the thread that distributes the load
				iret = pthread_create( &thread_main, 0 , &MultigirdSolver::MPI_distribute_Load , sgrid );
				// sent the sorted array to the main thread
				FITOB_MPI_Send_Int( &(sortedFullGrids_[0]) , sgrid->nrFullGrid() , 0 );
			}

			int fullGridIndex = 0 ;
			// iterate as long we do not get an index smaller than zero
			while (fullGridIndex >= 0)
			{
				// Receive the fullgridIndex
				FITOB_OUT_LEVEL3(verb(), "  Wait for root's command ");
				FITOB_MPI_Recv_Int( &fullGridIndex , 1 , 0 );

				// if this is smaller than zero then we stop
				if ( fullGridIndex < 0 ) { break; }

				DVector *vect = 0;
				// if the vector is not on the processor then we have to receive it
				if (FITOB_MPI_Comm_rank() > 0) {
					vect = new DVector(sgrid->getFG(fullGridIndex).totalPoints());
					vect->resize(sgrid->getFG(fullGridIndex).totalPoints());
					FITOB_MPI_Recv_Double( &((*vect)[0]), sgrid->getFG(fullGridIndex).totalPoints() , 0 );
					sgrid->getFG(fullGridIndex).setUnknVect(vect);
				}

				FITOB_OUT_LEVEL1(verb(), "  MultigirdSolver::solvePDE , subspace:" << fullGridIndex << " out of:" << sgrid->nrFullGrid());
				FITOB_OUT_LEVEL1(verb(), "  MultigirdSolver::solvePDE , process nr:" << FITOB_MPI_Comm_rank() << " out of:" << FITOB_MPI_Comm_size() );
				FITOB_OUT_LEVEL1(verb(),  sgrid->getFG(fullGridIndex).toStringFG() );
				// solve for this subspace (full grid)
				solvePDE_with_MG( &(sgrid->getFG(fullGridIndex)) , context , models , calc , oneCombiTimeStep , fullGridIndex);

				// send one INT to signal (index of the full grid) that I am ready
				FITOB_MPI_Send_Int( &fullGridIndex , 1 , 0 );

				// send the double vector to the rank 0
				FITOB_OUT_LEVEL4(verb(), "  MultigirdSolver::solvePDE , send FG:" << fullGridIndex );
				if (fitob::FITOB_MPI_Comm_rank() > 0) {
					FITOB_MPI_Send_Double( &((*vect)[0]) , sgrid->getFG(fullGridIndex).totalPoints() , 0 );
				}
				// delete the vector after we finished computing
				if (vect != 0 ) delete vect;
			}
			FITOB_OUT_LEVEL4(verb(), "  Finished iterating ");

			// rank 0 waits for the managing thread to finish
			if ( FITOB_MPI_Comm_rank() == 0 ) {
				FITOB_OUT_LEVEL4(verb(), "  Waiting for the managing thread to finish");
				pthread_join( thread_main , NULL);
			}
			// Synchronize all processes
			FITOB_MPI_Barrier();
 }
 else
		// ------------------- END of MPI parallel section ----------------
#endif
 {
			// ------------------- NO MPI parallel section ----------------
			// for cycle which can run parallel with OpenMP
            #if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_COMBI) )
            #pragma omp parallel
			{
            #pragma omp for schedule(dynamic, 1)
            #endif
			  for (int tt = 0 ; tt < sgrid->nrFullGrid() ; tt++)
			  {
				FITOB_OUT_LEVEL1(verb(), "  MultigirdSolver::solvePDE , subspace:" << tt << " out of:" << sgrid->nrFullGrid() << " index:" << sortedFullGrids_[tt]);
				FITOB_OUT_LEVEL1(verb(), "  MultigirdSolver::solvePDE , thread nr:" << FITOB_OPENMP_GET_THREAD_NUM() << " out of:" << FITOB_OPENMP_GET_MAX_THREADS());
				FITOB_OUT_LEVEL1(verb(),  sgrid->getFG(sortedFullGrids_[tt]).toStringFG() );
				// solve for this subspace (full grid)
				solvePDE_with_MG( &(sgrid->getFG(sortedFullGrids_[tt])) , context , models , calc , oneCombiTimeStep , tt);
			  }
			#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_COMBI) )
            }
            #endif
			// -------------------- END of NO MPI parallel section ---------
 }

			// apply combi technique and re compose the full grids
			sgrid->reCompose( combiErrorIndicator , true );

			FITOB_OUT_LEVEL1(verb(), " MultigirdSolver::solvePDE , reCompose() , combiErrorIndicator:" << combiErrorIndicator );
			// here we must add time to the context time
			globalCoordinates_backup[0] += oneCombiTimeStep;
			actualTime = actualTime + oneCombiTimeStep;
			co++;

			if (combiErrorIndicator > 1e-10)
			{nextStep = oneCombiTimeStep*(0.8*combiStepEps_/combiErrorIndicator);}
			else { nextStep = oneCombiTimeStep; }
			// limit the combi time step is necessary
			oneCombiTimeStep = (nextStep > maxCombiDT_) ? maxCombiDT_ : nextStep;
			oneCombiTimeStep = (oneCombiTimeStep > maxCombiDT_) ? maxCombiDT_ : oneCombiTimeStep;
			oneCombiTimeStep = (oneCombiTimeStep > (timeStep-actualTime)) ? (timeStep-actualTime) : oneCombiTimeStep;
			oneCombiTimeStep = (oneCombiTimeStep < minCombiDT_) ? minCombiDT_ : oneCombiTimeStep;

			FITOB_OUT_LEVEL1(verb(), " MultigirdSolver::solvePDE , new combi time step:" << oneCombiTimeStep << " , proposed:" << nextStep);
			FITOB_OUT_LEVEL1(verb(), " MultigirdSolver::solvePDE , maxCombiDT_:" << maxCombiDT_ << " , minCombiDT_:" << minCombiDT_);

		} // end of the while loop

	}
	else
	{
	   // --- call the function to solve with MG ----
	   FITOB_OUT_LEVEL3(verb(), " MultigirdSolver::solvePDE SOLVING FullGrid ")
	   FullGridBase* fg_base = dynamic_cast<FullGridBase*>(grid);
	   FITOB_ERROR_TEST( fg_base != NULL , "MultigirdSolver::solvePDE could not convert MeshBase to FullGridBase ");
	   solvePDE_with_MG( fg_base , context , models , calc , timeStep , 0 );

	}
}

void MultigirdSolver::solvePDE_with_MG(
		              FullGridBase *fgrid , const MeshContext* context ,
                      const ModelCollection* models ,
                      const FitobCalculator* calc ,
                      double timeStep ,
                      int nrFullGrid ){

    // -- now we do a simple GS iteration
	// --- create the highest grid
	double err = 1.0 , err_old = 0.0 , dtime = 0.0 , newDT = 0.0;
	// the pre-smooth and post smooth iteration number
	// (in one iteration there are two alternating iterations )
	presmooth_ = 1;
	postsmooth_ = 1;

	// this grid is the highest in the hierarchy
	boost::ptr_vector<MultigridFGBase> gridContain;
	boost::ptr_vector<SmootherBase> smoothContain;

	int depth = 0;
	MultigridFGBase* actGrid;

	// plot the FG
	if (verb() > 4) { fgrid->plotGrid( "out/FG_bef_MG" , nrFullGrid); }

	// insert the first grid and smoother on the highest level
	FITOB_OUT_LEVEL3(verb(), " MultigirdSolver::solvePDE_with_MG insert first FG_MG , NumThereads = " << FITOB_OPENMP_GET_MAX_THREADS() )
	gridContain.insert( gridContain.end() , MultigridFGBase::createMultigridFG(fgrid) );
	FITOB_OUT_LEVEL3(verb(), " MultigirdSolver::solvePDE_with_MG insert first GS ")
	smoothContain.insert( smoothContain.end() ,  SmootherBase::createSmootherPointer(conf() , &(gridContain[depth]) , models , calc ) );

	actGrid = &(gridContain[depth]);
	// create the hierarchy of grids and smoothers
	while ( (actGrid->caBeRefined()) ){ //|| (depth < 1)
		FITOB_OUT_LEVEL3(verb(), " MultigirdSolver::solvePDE_with_MG insert next FG_MG ")
		// create the FG for the actual level
		gridContain.insert( gridContain.end() , MultigridFGBase::createMultigridFG( actGrid , true ) );
		// create the smoother for the grid
		smoothContain.insert( smoothContain.end() ,  SmootherBase::createSmootherPointer( conf() , &(gridContain[depth+1]) , models , calc ) );

		actGrid = &(gridContain[depth+1]);
		depth++;
	}

	// copy global coords
    DVector globCoord = globalCoordinates_backup;
    DVector unknowns( gridContain[0].u().size() , 0.0 );
    DVector unknowns_backup( gridContain[0].u().size() , 0.0 );
    double predictor_norm_L2 , predictor_norm_Inf , y_max;

	// set the time step
    double rest_time = timeStep;
    int count = 0;
    bool acceptActualStep = true;
	dtime = newDT = minMicroDT_;

    //for (int dt = 0 ; dt < nr_timeStep ; dt++)
    while ( rest_time > FITOB_NUMERICALZERO )
    {
    	// limit the time step to the limits
    	dtime = (dtime < minMicroDT_) ? minMicroDT_ : dtime;
    	dtime = (dtime > maxMicroDT_) ? maxMicroDT_ : dtime;
    	dtime = (dtime > rest_time ) ? rest_time : dtime;

    	for (unsigned int hh = 0 ; hh < smoothContain.size() ; hh++)
    		smoothContain[hh].setTimeStep( dtime );

        // debugging plot of the grid
        if (verb() > 5) { gridContain[0].plotMAT_grad("out/solution_before_MG" , gridContain[0].u() , nrFullGrid ); }

    	// ---------------- TIME STEPPING ---------
        if (usePredictor_) {
        	// calculate the RHS so that we can make one explicit time step
        	smoothContain[0].calculateRHS( globCoord );
			//gridContain[0].plotMAT_grad("unknown_bef_expl" , gridContain[0].u() , 0 );
			smoothContain[0].explicitStep( globCoord , unknowns , predictorUnderrelaxCoef_ );
			//gridContain[0].plotMAT_grad("unknown_aft_expl" , gridContain[0].u() , 0 );
        }
    	// ---------------- END TIME STEPPING ---------

        // backup the actual solution for restart purpose
       	for (unsigned int i = 0 ; i < unknowns.size() ; i++){
       		unknowns_backup[i] = gridContain[0].u()[i];
        }

	    err = err_old = 1.0; count = 0; acceptActualStep = true;
	    // this is the GS cycle
        while (err > solverEpsilon_ ){

           //  call MG recursively
           recursive_V_cycle( gridContain , smoothContain ,
        			depth , 0 , 0 , globCoord);

		   err = smoothContain[0].calcResiduum( globCoord );
		   FITOB_OUT_LEVEL2(verb()," MultigirdSolver::solvePDE_with_MG: , after one V Cycle , err:" << err);

		   count++;
		   // if result is not achieved in a given number of iterations, then smooth aggressively
		   if (count > 15){
				presmooth_ = presmooth_ + 1;
				postsmooth_ = postsmooth_ + 3;
				count = 0;
		   }
		   // if there is divergence then smooth more aggressively
		   if (err > err_old){
				presmooth_ = presmooth_ + 1;
				postsmooth_ = postsmooth_ + 10;
		   }
		   err_old = err;
        }

	    // ---------------- TIME STEPPING -------------
	    // calculate norm, according to the time step size control formula
        if (makeTimeStepControll_){
        	double choosen_norm = 0.0;
        	for (unsigned int i = 0 ; i < unknowns.size() ; i++){
        		unknowns[i] = (2.0/3.0) * ( gridContain[0].u()[i] - unknowns[i] );
        	}
        	predictor_norm_Inf = fitob::inf_norm(&(unknowns));
        	predictor_norm_L2 = fitob::l2_norm(&(unknowns));
        	// calc the maximum of the solution
            y_max = fitob::inf_norm(&(unknowns_backup));
        	//FITOB_OUT_LEVEL3(verb()," MultigirdSolver::solvePDE_with_MG: , time step norm L2: " << predictor_norm_L2 << " Inf:" << predictor_norm_Inf);
        	// choose between the norms
        	if (useInfNormForTimeStepControll_){
        		choosen_norm = predictor_norm_Inf;
        	} else{
        		choosen_norm = predictor_norm_L2;
        	}
        	// apply the form for the next time step
        	newDT = dtime * pow ( ( y_max * timeStepControll_Epsilon_ / choosen_norm) , 1.0/3.0);
        	// test when we accept this step
        	if ( (10.0*newDT < dtime) && ( dtime > 1.2*minMicroDT_) ){
        		// do not accept this step and redo the step
               	for (unsigned int i = 0 ; i < unknowns.size() ; i++){
               		gridContain[0].u()[i] = unknowns_backup[i];
                }
           		acceptActualStep = false;
        	}
        }
	    // ---------------- END TIME STEPPING ---------

        // increase time in the global coordinates, if this step is accepted
        if (acceptActualStep){
        	globCoord[0] += dtime;
        	rest_time = rest_time - dtime;
    	    FITOB_OUT_LEVEL2(verb()," MultigirdSolver::solvePDE_with_MG: , after solve , rest time:" << rest_time << " dt:" << dtime);
        } else {
    	    FITOB_OUT_LEVEL2(verb()," MultigirdSolver::solvePDE_with_MG: , rejecting actual step , redo with dt: " << newDT << " old dt:" <<dtime);
        }

        // set the new time step
        if (makeTimeStepControll_){
        	FITOB_OUT_LEVEL2(verb()," MultigirdSolver::solvePDE_with_MG: , new time step: " << newDT << " , old dt : " << dtime );
        	dtime = newDT;
        }

        // debugging plot of the grid
        if (verb() > 5) { gridContain[0].plotMAT_grad("out/solution_after_MG" , gridContain[0].u() , nrFullGrid ); }

	    // impose constraints , if there are any
		if (calc->getScriptModel().get()->constrainOpSeq().nrOperators() > 0){
			globCoord = globalCoordinates_backup;
			gridContain[0].applyConstraints( &(calc->getScriptModel().get()->constrainOpSeq()) , globCoord );
		}
    }

	// write from mfg the results back
    gridContain[0].writeBackSolution(fgrid);

	// plot the FG
	if (verb() > 4) { fgrid->plotGrid( "out/FG_aft_MG" , nrFullGrid); }

	FITOB_OUT_LEVEL4(verb()," MultigirdSolver::solvePDE: , after create mfg.writeBackSolution");

}

void MultigirdSolver::recursive_V_cycle(
		            boost::ptr_vector<MultigridFGBase> &gridContain ,
					boost::ptr_vector<SmootherBase> &smoothContain ,
					int maxDepth ,
					int actualDepth ,
					int count ,
					DVector &globCoord){

	double err = 1.0;

    // do iteration pre-smoothing
	for (int i = 0 ; i < presmooth_ ; i++){  smoothContain[actualDepth].smoothGrid( globCoord ); }

	if (!(actualDepth >= maxDepth)) {
		// calculate the error
		err = smoothContain[actualDepth].calcResiduum( globCoord );
		//FITOB_OUT_LEVEL3(verb()," MultigirdSolver::solvePDE_with_MG: ,before 2 grid , err:" << err);

		//gridContain[actualDepth].plotMAT_grad("residum_fine_" , mfg.res() , count );

		// this is the classical correction scheme
		ProlongationRestriction::makeDirectRestriction(
				&(gridContain[actualDepth]) , gridContain[actualDepth].res() ,
				&(gridContain[actualDepth+1]) , gridContain[actualDepth+1].rhs() , 1.0 , 0.0 );

		//gridContain[actualDepth+1].plotMAT_grad("restrict" , mfg1.rhs() , count );

		//set the unknowns to zero
		for (unsigned int ii = 0 ; ii < gridContain[actualDepth+1].u().size() ; ii++)
			gridContain[actualDepth+1].u()[ii] = 0.0;

		smoothContain[actualDepth+1].residumIsSet();
		if (actualDepth+1 < maxDepth ){
			recursive_V_cycle( gridContain , smoothContain , maxDepth , (actualDepth+1) , count , globCoord);
		}

		//gridContain[actualDepth+1].plotMAT_grad("restrict_c" , gridContain[actualDepth+1].u() , count );
		//gridContain[actualDepth].plotMAT_grad("restrict_bef_f" , gridContain[actualDepth].u() , count );

		// add the correction to the current unknown vector
		ProlongationRestriction::makeLinearProlongation(
				&(gridContain[actualDepth]) , gridContain[actualDepth].u() ,
				&(gridContain[actualDepth+1]) , gridContain[actualDepth+1].u() , 1.0 , 1.0 );

		//gridContain[actualDepth].plotMAT_grad("restrict_aft_f" , gridContain[actualDepth].u() , count );

		//err = smoothContain[actualDepth].calcResiduum( globCoord );
		//FITOB_OUT_LEVEL3(verb()," MultigirdSolver::solvePDE_with_MG: ,after 2 grid  , err:" << err);
	}

   // do postsmooth_ iteration post smoothing
   for (int i = 0 ; i < postsmooth_ ; i++){  smoothContain[actualDepth].smoothGrid( globCoord ); }

   //err = smoothContain[actualDepth].calcResiduum( globCoord );
   //FITOB_OUT_LEVEL3(verb()," MultigirdSolver::solvePDE_with_MG: , after all gs.smoothGrid , err:" << err);

}


// ===================== OLD DEPREDICATED SOLVER , GAUSS SEIDEL ==============================

// ============ WITH BOUNDARY =============
void MultigirdSolver::solvePDE_with_GS(
		      FullGridBase *fgrid , const MeshContext* context ,
              const ModelCollection* models ,
              const FitobCalculator* calc ,
              double timeStep){

    // -- now we do a simple GS iteration
	// --- create the highest grid
	double err = 1.0 , dtime = 0.0;
	// here we divide the standard time step in 7 small ones
    MultigridFG mfg(fgrid);
    FITOB_OUT_LEVEL3(verb()," MultigirdSolver::solvePDE: , after create MultigridFG");
    GaussSeidel  gs( conf() , &(mfg) , models , calc);
    FITOB_OUT_LEVEL3(verb()," MultigirdSolver::solvePDE: , after create GaussSeidel");
    DVector globCoord = globalCoordinates_backup;
	// todo: make a time step controlling
    int nr_timeStep = ceil( timeStep / minMicroDT_ );
    for (int dt = 0 ; dt < nr_timeStep ; dt++)
    {
    	// set the time step
    	dtime = timeStep/ (double)nr_timeStep;
	    gs.setTimeStep( dtime );
	    err = 1.0;
	    // this is the GS cycle
        while (err > solverEpsilon_ ){
           // do 10 iteration till residuum calculation
		   for (int i = 0 ; i < 10 ; i++){
			  gs.smoothGrid( globCoord );
		   }
		   err = gs.calcResiduum( globCoord );
		   FITOB_OUT_LEVEL3(verb()," MultigirdSolver::solvePDE: , after create gs.smoothGrid , err:" << err);
        }
        // increase time in the global coordinates
        globCoord[0] += dtime;
	    FITOB_OUT_LEVEL3(verb()," MultigirdSolver::solvePDE: , after solve");

	    // impose constraints , if there are any
		if (calc->getScriptModel().get()->constrainOpSeq().nrOperators() > 0){
	          mfg.applyConstraints( &(calc->getScriptModel().get()->constrainOpSeq()) , globCoord );
		}
    }

	// write from mfg the results back
	mfg.writeBackSolution(fgrid);

	FITOB_OUT_LEVEL3(verb()," MultigirdSolver::solvePDE: , after create mfg.writeBackSolution");
}

/*
 * FitobExtrapGrid.cpp
 *
 *  Created on: Aug 27, 2010
 *      Author: benk
 */

#include "FitobExtrapGrid.hpp"
#include "src/mesh/FitobFullGrid.hpp"
#include "src/mesh/FitobFullGrid_WB.hpp"

// the SGpp interface
#include "combigrid.hpp"

using namespace fitob;
using namespace std;

ExtrapGrid::ExtrapGrid(const Domain* dom , int sparseGridType , double diagonalCutOffLevel ,
		bool use_Opticom , const DVector& adaptiveTruncation ,  GridType gridType )
: SparseGrid( dom , " fitob::ExtrapGrid" , use_Opticom) ,
		combinationGrid_(0) , combiScheme_(0), nrFullGrids_(-1) ,
	global_Max_level_(0) , sparseGridType_(sparseGridType) ,
	diagonal_cut_of_level_(1) {

	fullGridLevels_.resize(0);
    setGridType(gridType);
	setVerb(0);

	// get the "global" level (this will be used if we do not have dimension adaptivity)
	global_Max_level_ = domain()->getAxisLevel( domain()->localToGlobalIndex(0));
	IVector adaptLevels(domain()->nrRealAxis(),0 );
	IVector adaptTruncation(domain()->nrRealAxis() , 0 );

	for (int ii = 0 ; ii < domain()->nrRealAxis() ; ii++)
	{
		adaptLevels[ii] = (size_t) domain()->getAxisLevel( domain()->localToGlobalIndex(ii) );
		global_Max_level_ = (global_Max_level_ < domain()->getAxisLevel( domain()->localToGlobalIndex(ii)) ) ?
				domain()->getAxisLevel( domain()->localToGlobalIndex(ii)) : global_Max_level_;
	}

	// choose between sparse grid types
	switch(sparseGridType){
	case 16: {
		// classical sparse grid , with classical combination technique
		FITOB_OUT_LEVEL3(verb(),"ExtrapGrid::ExtrapGrid Ctor create Combi LinearTrapezoidBoundaryGrid , " << global_Max_level_ << " , " << domain()->nrRealAxis() );
		//combinationGrid_ = (new sg::FullGridSet( domain()->nrRealAxis() , global_Max_level_  , "linearTrapezoidBoundary"));
		combiScheme_ = (new combigrid::S_CT( domain()->nrRealAxis() , adaptLevels ));
		combinationGrid_ = (new combigrid::SerialCombiGrid(combiScheme_ , (getGridType() == GRID_WITH_BOUNDARY) ));
		combinationGrid_->createFullGrids();
	} break;
	case 17: {
		// create square root sparse grid
		FITOB_OUT_LEVEL3(verb(),"ExtrapGrid::ExtrapGrid Ctor create Combi SquareRootGrid , " << global_Max_level_ << " , " << domain()->nrRealAxis() );
		combiScheme_ = (new combigrid::TS_CT( domain()->nrRealAxis() , adaptLevels ));
		combinationGrid_ = combinationGrid_ = (new combigrid::SerialCombiGrid(combiScheme_ ,  (getGridType() == GRID_WITH_BOUNDARY) ));
#if defined(FITOB_MPI)
		// the grid should be created only on the 0-th rank processor
		if (fitob::FITOB_MPI_Comm_rank() < 1) { combinationGrid_->createFullGrids(); }
#else
		combinationGrid_->createFullGrids();
#endif
	} break;
	case 18: {
		  if (adaptiveTruncation.size() <= 0){
			// create modified combi technique (cut off the diagonal)
			FITOB_OUT_LEVEL3(verb(),"ExtrapGrid::ExtrapGrid Ctor create Modified Combi Grid , "
					<< global_Max_level_ << " , " << domain()->nrRealAxis() );
			// this level defines which part of the "diagonal" will be cut off
			diagonal_cut_of_level_ = ceil((diagonalCutOffLevel < 0.0) ? (double)global_Max_level_/fabs(diagonalCutOffLevel) : diagonalCutOffLevel) ;
			for (int ij = 0 ; ij < domain()->nrRealAxis() ; ij++){
				FITOB_OUT_LEVEL3(verb(),"NONE , " << ij << " , " << domain()->nrRealAxis() );
				adaptTruncation[ij] = diagonal_cut_of_level_;
				FITOB_OUT_LEVEL3(verb(),"NONE1 , " << ij << " , " << domain()->nrRealAxis() );
				FITOB_OUT_LEVEL3(verb(), "ij=" << ij << ", adaptTruncation[ij]=" << adaptTruncation[ij]);
				FITOB_OUT_LEVEL3(verb(),"DOME , " << ij << " , " << domain()->nrRealAxis() );
			}
		 }
		 else {
			// dimension adaptivity truncation for T-CT
			// create modified combi technique (cut off the diagonal)
			FITOB_OUT_LEVEL3(verb(),"ExtrapGrid::ExtrapGrid Ctor create Modified Combi Grid (adaptive truncation) , "
					<< global_Max_level_ << " , " << domain()->nrRealAxis() );
			FITOB_ERROR_TEST((int)adaptiveTruncation.size() == (int)domain()->nrRealAxis() , " ExtrapGrid::ExtrapGrid Ctor (adaptive truncation) v1:"
					<< adaptiveTruncation.size() << " , v2:" << domain()->nrRealAxis() );
			// this level defines which part of the "diagonal" will be cut off
			for (int ij = 0 ; ij < domain()->nrRealAxis() ; ij++){
				adaptTruncation[ij] = ceil((adaptiveTruncation[ij] < 0.0) ? (double)adaptLevels[ij]/fabs(adaptiveTruncation[ij]) : adaptiveTruncation[ij] ) ;
				FITOB_OUT_LEVEL3(verb(), "ij="<<ij<<" ,  adaptLevels[ij]="<<adaptLevels[ij]<< " , adaptiveTruncation[ij] =" << adaptiveTruncation[ij] <<
						", adaptTruncation[ij]=" << adaptTruncation[ij]);
			}
		 }
		 // here we create the extrapolation grid
	     FITOB_OUT_LEVEL3(verb(),"Before creating stuff" );
		 combiScheme_ = (new combigrid::S_CT( domain()->nrRealAxis() , adaptLevels , adaptTruncation ));
		 combinationGrid_ = (new combigrid::SerialCombiGrid(combiScheme_ , (getGridType() == GRID_WITH_BOUNDARY) ));
#if defined(FITOB_MPI)
		 FITOB_OUT_LEVEL3(verb(),"ExtrapGrid::ExtrapGrid Ctor create only rank 0" );
		// the grid should be created only on the 0-th rank processor
		if (fitob::FITOB_MPI_Comm_rank() < 1) { combinationGrid_->createFullGrids(); }
		 FITOB_OUT_LEVEL3(verb(),"ExtrapGrid::ExtrapGrid Ctor END create only rank 0" );
#else
		combinationGrid_->createFullGrids();
#endif
		} break;
	}

}

ExtrapGrid::~ExtrapGrid() {
	if (combinationGrid_ != 0){
		delete combinationGrid_;
		delete combiScheme_;
	}
}


double ExtrapGrid::eval(const DVector& globalCoords) const {

	DVector localCoords(domain()->nrRealAxis());
	DVector intersect(domain()->nrRealAxis());
	IVector minIndex(domain()->nrRealAxis());
	IVector maxIndex(domain()->nrRealAxis());
	IVector localCoord(domain()->nrRealAxis());
	int middle = 0;
	int nrAxis_ = domain()->nrRealAxis();
	// unit coordinates
	DVector unitCoords(domain()->nrRealAxis(),0.0);

	// get the local coordinates
	domain()->globalToLocal( globalCoords , localCoords);

	// find the position on each axis
	for (int ii = 0 ; ii < nrAxis_ ; ii++)
	{
		int globalIndex = domain()->localToGlobalIndex(ii);
		const DVector& axisGrading = domain()->getGradedAxis( globalIndex );
		minIndex[ii] = 0;
		maxIndex[ii] = axisGrading.size()-1;

		// this should have a complexity of O(logN) or O(L), L-> level
		// bisection search, since we assume that the array is ordered
		//FITOB_OUT_LEVEL6(verb()," ExtrapGrid::eval , globalCoords[globalI]:" << globalCoords[domain()->localToGlobalIndex(ii)]);
		//FITOB_OUT_LEVEL6(verb()," ExtrapGrid::eval , localCoords[ii]:" << localCoords[ii]);
		for (; ;){
		  middle = (minIndex[ii] + maxIndex[ii])/2;
		  //FITOB_OUT_LEVEL6(verb()," ExtrapGrid::eval , ii: " << ii << ", middle:" << middle <<
			//	  ", axisGrading_[ii][middle]:" << axisGrading[middle]);
		  if (localCoords[ii] >= axisGrading[middle] ){
			  minIndex[ii] = middle;
			  //FITOB_OUT_LEVEL6(verb()," minIndex[ii] = " << minIndex[ii]);
		  }
		  else{
			  maxIndex[ii] = middle;
			  //FITOB_OUT_LEVEL6(verb()," maxIndex[ii] = " << maxIndex[ii]);
		  }
          // break the for cycle when we are close enough
		  if ((maxIndex[ii] - minIndex[ii]) < 2) break;
		}

		// calculate the linear
		intersect[ii] = (axisGrading[maxIndex[ii]] - localCoords[ii]) /
				(axisGrading[maxIndex[ii]] - axisGrading[minIndex[ii]]);

		//   - make the transformation to the unit domain -> unitCoords
		unitCoords[ii] = ((double)maxIndex[ii] - intersect[ii]) /
				fitob::powerTwo[ domain()->getAxisLevel(globalIndex) ];

		//FITOB_OUT_LEVEL6(verb()," ExtrapGrid::eval , ii: " << ii << ", level[ii]:" << domain()->getAxisLevel(globalIndex) );
		//FITOB_OUT_LEVEL6(verb()," ExtrapGrid::eval , ii: " << ii << ", localCoords[ii]:" << localCoords[ii]);
		//FITOB_OUT_LEVEL6(verb()," ExtrapGrid::eval , ii: " << ii << ", intersect[ii]:" << intersect[ii]);
		//FITOB_OUT_LEVEL6(verb()," ExtrapGrid::eval , ii: " << ii << ", unitCoords[ii]:" << unitCoords[ii]);
		//FITOB_OUT_LEVEL6(verb()," ExtrapGrid::eval , ii: " << ii << ", axisGrading[minIndex[ii]]:" << axisGrading[minIndex[ii]]);
		//FITOB_OUT_LEVEL6(verb()," ExtrapGrid::eval , ii: " << ii << ", axisGrading[maxIndex[ii]]:" << axisGrading[maxIndex[ii]]);
		//FITOB_OUT_LEVEL6(verb()," ExtrapGrid::eval , ii: " << ii << ", minIndex[ii]:" << minIndex[ii]);
		//FITOB_OUT_LEVEL6(verb()," ExtrapGrid::eval , ii: " << ii << ", maxIndex[ii]:" << maxIndex[ii]);
	}

    // eval the combi grid on the unit coordinates
	double eval_val = combinationGrid_->eval( unitCoords );
	FITOB_OUT_LEVEL5(verb()," ExtrapGrid::eval RET_VAL:" << eval_val );
	return eval_val;
}

void ExtrapGrid::eval(const std::vector<DVector>& globalCoordonates , DVector& resVector) const {
	   for (unsigned int ii = 0 ; ii < globalCoordonates.size() ; ii++){
		   // we just call the simple version of this function
		   // todo: this could be optimized by copy paste the body of the function above
		   resVector[ii] = this->eval(globalCoordonates[ii]);
	   }
}

void ExtrapGrid::setValues(const Evaluable* func , const DVector& globalCoords){

	double dummy = 0.0;
	// get the values from SGpp::FullGrids
	deCompose();

#if defined(FITOB_MPI)
	// Synchronize all processes
	FITOB_MPI_Barrier();
	// if we are not on the 0-th processor then just exit
	if (fitob::FITOB_MPI_Comm_rank() >= 1) {
		FITOB_OUT_LEVEL3(verb()," ExtrapGrid::setValues() , MPI END ");
		return;
	}
#endif

	FITOB_OUT_LEVEL3(verb()," ExtrapGrid::setValues set values");

	// this for loop can be parallelized so that each full grid sets its value separately
#if defined(FITOB_OPENMP)
#pragma omp parallel
{
#pragma omp for schedule(dynamic, 1)
#endif
	for ( int i=0; i < nrFullGrids_ ; i++)
	{
		DVector tmpGlCoo  = globalCoords;
		FullGridBase& fitobFG = fullGridContainer_[i];
	    FITOB_OUT_LEVEL3(verb()," ExtrapGrid::setValues FullGrid nr. i:" << i );
	    fitobFG.setValues( func , tmpGlCoo );
    }
#if defined(FITOB_OPENMP)
}
#endif

	// put the values back to the SGpp::FullGrids
	reCompose(dummy , false);
	FITOB_OUT_LEVEL3(verb()," ExtrapGrid::setValues END");
}

void ExtrapGrid::applyConstraints(const OperatorSequence* constraintOpSeq_ ,
		                      const DVector& globalCoords ) {

	double dummy = 0.0;
	// get the values from SGpp::FullGrids
	deCompose();

#if defined(FITOB_MPI)
	// Synchronize all processes
	FITOB_MPI_Barrier();
	// if we are not on the 0-th processor then just exit
	if (fitob::FITOB_MPI_Comm_rank() >= 1) {
		FITOB_OUT_LEVEL3(verb()," ExtrapGrid::applyConstraints() , MPI END ");
		return;
	}
#endif


	FITOB_OUT_LEVEL3(verb()," ExtrapGrid::applyConstraints apply constraints ");

	// this for loop can be parallelized so that each full grid sets its value separately
#if defined(FITOB_OPENMP)
#pragma omp parallel
{
#pragma omp for schedule(dynamic, 1)
#endif
	  for ( int i=0; i < nrFullGrids_ ; i++)
	  {
		  DVector tmpGlCoo  = globalCoords;
		  FullGridBase& fitobFG = fullGridContainer_[i];
	      FITOB_OUT_LEVEL3(verb()," ExtrapGrid::applyConstraints apply constraints for FullGrid nr. i:" << i );
	      fitobFG.applyConstraints( constraintOpSeq_ , tmpGlCoo );
      }
#if defined(FITOB_OPENMP)
}
#endif

	// put the values back to the SGpp::FullGrids
	reCompose( dummy , false);
	FITOB_OUT_LEVEL3(verb()," ExtrapGrid::applyConstraints END");
}

void ExtrapGrid::deCompose(){

	FITOB_OUT_LEVEL3(verb()," ExtrapGrid::deCompose START");

	// create the full grids if nor existing before
	if (fullGridLevels_.size() < 1)
	{
		// resize the levels
		FITOB_OUT_LEVEL3(verb()," ExtrapGrid::deCompose() , CREATE COMBI GRID ");
		fullGridLevels_.resize( (int)combinationGrid_->getNrFullGrid() );
		// create the grid
		nrFullGrids_ = combinationGrid_->getNrFullGrid();
		for ( int i=0 ; i < nrFullGrids_ ; i++ )
		{
		   	IVector specialLevels( domain()->nrRealAxis() );
		   	combigrid::FullGridD *fg = (combinationGrid_->getFullGrid(i));
		   	// copy the levels
		   	for (int jj = 0 ; jj < domain()->nrRealAxis() ; jj++ ){
		   		FITOB_OUT_LEVEL3(verb()," ExtrapGrid::deCompose , create levels jj:"<<jj<<" fg->getLevel()[jj]: " << fg->getLevels()[jj] );
		   		specialLevels[jj] = fg->getLevels()[jj];
		   	}
		   	// create fullgrid
		   	FITOB_OUT_LEVEL3(verb()," ExtrapGrid::deCompose() , create FULLGRID nrUnkowns:" << fg->getNrElements());

		   	//create FGs in case, which do not have boundary points
		   	switch (getGridType()){
		   	   case (GRID_WITH_BOUNDARY):{
				   	fullGridContainer_.insert( fullGridContainer_.end() ,
				   	    			new FullGrid( domain() , specialLevels , &(fg->getElementVector()) ) );
		   	   break;}
		   	   case (GRID_WITHOUT_BOUNDARY):{
				   	fullGridContainer_.insert( fullGridContainer_.end() ,
				   	    			new FullGrid_WB( domain() , specialLevels , &(fg->getElementVector()) ) );
		   	   break;}
		   	}

		   	// store the levels for this full grid
		   	fullGridLevels_[i] = specialLevels;
		}
	}

#if defined(FITOB_MPI)
	// if we are not on the 0-th processor then just exit
	if (fitob::FITOB_MPI_Comm_rank() >= 1) {
		FITOB_OUT_LEVEL3(verb()," ExtrapGrid::deCompose() , MPI END ");
		return ;
	}
#endif

	// transform the SGpp full grids into fitob::FullGrids
	FITOB_OUT_LEVEL3(verb()," ExtrapGrid::deCompose() , START DE COMBI ");
	for ( int i=0 ; i < combinationGrid_->getNrFullGrid() ; i++ )
	{
		// sg::FullGrid is also the base class for the full grid without boundary points
		combigrid::FullGridD *fg = (combinationGrid_->getFullGrid(i));
	   	int  m = fg->getNrElements();
	   	FullGridBase& fitobFG = fullGridContainer_[i];

	   	// test if the size is not equal
	   	FITOB_ERROR_TEST( (int)fitobFG.unknVect().size() == m , " ExtrapGrid::deCompose , fitobFG.unknVect().size()" <<
	   			fitobFG.unknVect().size() << " , m:" << m );

	   	// copy the vector values
	   	FITOB_OUT_LEVEL3(verb()," ExtrapGrid::deCompose() START i:" << i << " m:"<<m);
	   	//for ( unsigned int j = 0 ; j < m ; j++) {
	   		//fitobFG.unknVect()[j] = (*fg)[j];
	   		//if (verb() > 3) std::cout << ", " << fitobFG.unknVect()[j];
	   	//}
	   	//if (verb() > 3) std::cout << std::endl;
	   	FITOB_OUT_LEVEL3(verb()," ExtrapGrid::deCompose() END i:" << i << " m:"<<m);
	}
	FITOB_OUT_LEVEL3(verb()," ExtrapGrid::deCompose() , END DE COMBI ");
}



void ExtrapGrid::reCompose(double &errorIndicator, bool change_coef){

#if defined(FITOB_MPI)
	// if we are not on the 0-th processor then just exit
	if (fitob::FITOB_MPI_Comm_rank() >= 1) {
		FITOB_OUT_LEVEL3(verb()," ExtrapGrid::reCompose() , MPI END ");
		return ;
	}
#endif

	FITOB_OUT_LEVEL3(verb()," ExtrapGrid::reCompose START ");

	// create the combi object if necessary
	FITOB_ERROR_TEST( combinationGrid_ != 0 , " ExtrapGrid::reCompose , combinationGrid_ is NULL , ");
	FITOB_ERROR_TEST( fullGridLevels_.size() > 0 , " ExtrapGrid::reCompose , fullGridLevels_.size() > 0 ");

	// transform the SGpp full grids into fitob::FullGrids
	FITOB_OUT_LEVEL3(verb()," ExtrapGrid::deCompose() , START COMBI ");
	for ( int i=0 ; i < combinationGrid_->getNrFullGrid() ; i++ )
	{
		// sg::FullGrid is also the base class for the full grid without boundary points
		combigrid::FullGridD *fg = (combinationGrid_->getFullGrid(i));
		int m = fg->getNrElements();
	   	FullGridBase& fitobFG = fullGridContainer_[i];
	   	// test if the size is not equal
	   	FITOB_ERROR_TEST( (int)fitobFG.unknVect().size() == m , " ExtrapGrid::deCompose , fitobFG.unknVect().size()" <<
	   			fitobFG.unknVect().size() << " , m:" << m );
	}
	FITOB_OUT_LEVEL3(verb()," ExtrapGrid::deCompose() , END COMBI ");
}

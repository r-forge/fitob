/*
 * FitobCombiGrid.cpp
 *
 *  Created on: Aug 27, 2010
 *      Author: benk
 */

#include "FitobCombiGrid.hpp"
#include "src/mesh/FitobFullGrid.hpp"
#include "src/mesh/FitobFullGrid_WB.hpp"

// the SGpp interface
#include "combigrid.hpp"

using namespace fitob;
using namespace std;
#ifdef SGPP_DIRECT_SOLVER
using namespace sg::base;
#endif

CombiGrid::CombiGrid(const Domain* dom , int sparseGridType , double diagonalCutOffLevel ,
		bool use_Opticom , const DVector& adaptiveTruncation ,  GridType gridType )
: SparseGrid( dom , " fitob::CombiGrid" , use_Opticom) ,
	combinationGrid_(0) , combiScheme_(0),
#ifdef SGPP_DIRECT_SOLVER
	gridstorageSGpp_(0) , alphas_(0) , minAlpha_(0), maxAlpha_(0),
#endif
	nrFullGrids_(-1) , global_Max_level_(0) , sparseGridType_(sparseGridType) ,
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
	case 6: {
		// classical sparse grid , with classical combination technique
		FITOB_OUT_LEVEL3(verb(),"CombiGrid::CombiGrid Ctor create Combi LinearTrapezoidBoundaryGrid , " << global_Max_level_ << " , " << domain()->nrRealAxis() );
		//combinationGrid_ = (new sg::FullGridSet( domain()->nrRealAxis() , global_Max_level_  , "linearTrapezoidBoundary"));
		combiScheme_ = (new combigrid::S_CT( domain()->nrRealAxis() , adaptLevels ));
		combinationGrid_ = (new combigrid::SerialCombiGrid(combiScheme_ , (getGridType() == GRID_WITH_BOUNDARY) ));
		combinationGrid_->createFullGrids();
		} break;
	case 7: {
		// create square root sparse grid
		FITOB_OUT_LEVEL3(verb(),"CombiGrid::CombiGrid Ctor create Combi SquareRootGrid , " << global_Max_level_ << " , " << domain()->nrRealAxis() );
		combiScheme_ = (new combigrid::TS_CT( domain()->nrRealAxis() , adaptLevels ));
		combinationGrid_ = combinationGrid_ = (new combigrid::SerialCombiGrid(combiScheme_ ,  (getGridType() == GRID_WITH_BOUNDARY) ));
		combinationGrid_->createFullGrids();
		} break;
	case 8: {
		  if (adaptiveTruncation.size() <= 0){
			// create modified combi technique (cut off the diagonal)
			FITOB_OUT_LEVEL3(verb(),"CombiGrid::CombiGrid Ctor create Modified Combi Grid , "
					<< global_Max_level_ << " , " << domain()->nrRealAxis() );
			// this level defines which part of the "diagonal" will be cut off
			diagonal_cut_of_level_ = ceil((diagonalCutOffLevel < 0.0) ? (double)global_Max_level_/fabs(diagonalCutOffLevel) : diagonalCutOffLevel) ;
			for (int ij = 0 ; ij < domain()->nrRealAxis() ; ij++){
				adaptTruncation[ij] = diagonal_cut_of_level_;
				FITOB_OUT_LEVEL3(verb(), "ij="<<ij<<" ,  adaptLevels[ij]="<<adaptLevels[ij]<< " , adaptiveTruncation[ij] =" << adaptiveTruncation[ij] <<
						", adaptTruncation[ij]=" << adaptTruncation[ij]);
			}
		 } else {
				// dimension adaptivity truncation for T-CT
				// create modified combi technique (cut off the diagonal)
				FITOB_OUT_LEVEL3(verb(),"CombiGrid::CombiGrid Ctor create Modified Combi Grid (adaptive truncation) , "
						<< global_Max_level_ << " , " << domain()->nrRealAxis() );
				FITOB_ERROR_TEST((int)adaptiveTruncation.size() == (int)domain()->nrRealAxis() , " CombiGrid::CombiGrid Ctor (adaptive truncation) v1:"
						<< adaptiveTruncation.size() << " , v2:" << domain()->nrRealAxis() );
				// this level defines which part of the "diagonal" will be cut off
				for (int ij = 0 ; ij < domain()->nrRealAxis() ; ij++){
					adaptTruncation[ij] = ceil((adaptiveTruncation[ij] < 0.0) ? (double)adaptLevels[ij]/fabs(adaptiveTruncation[ij]) : adaptiveTruncation[ij] ) ;
					FITOB_OUT_LEVEL3(verb(), "ij="<<ij<<" ,  adaptLevels[ij]="<<adaptLevels[ij]<< " , adaptiveTruncation[ij] =" << adaptiveTruncation[ij] <<
							", adaptTruncation[ij]=" << adaptTruncation[ij]);
				}
		 }
		 // here we create the extrapolation grid
		 combiScheme_ = (new combigrid::S_CT( domain()->nrRealAxis() , adaptLevels , adaptTruncation ));
		 combinationGrid_ = (new combigrid::SerialCombiGrid(combiScheme_ , (getGridType() == GRID_WITH_BOUNDARY) ));
		 combinationGrid_->createFullGrids();
		} break;
	}
}

CombiGrid::~CombiGrid() {
	if (combinationGrid_ != 0){
#ifdef SGPP_DIRECT_SOLVER
		if (gridstorageSGpp_ != 0){
		   	// delete the grid storage
	    	delete alphas_;
	    	delete minAlpha_;
	    	delete maxAlpha_;
			delete gridstorageSGpp_;
		}
#endif
		delete combinationGrid_;
		delete combiScheme_;
	}
}


double CombiGrid::eval(const DVector& globalCoords) const {

	DVector localCoords(domain()->nrRealAxis());
	DVector intersect(domain()->nrRealAxis());
	IVector minIndex(domain()->nrRealAxis());
	IVector maxIndex(domain()->nrRealAxis());
	IVector localCoord(domain()->nrRealAxis());
	int middle = 0;
	int nrAxis_ = domain()->nrRealAxis();

	DVector unitCoords(domain()->nrRealAxis());

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
		//FITOB_OUT_LEVEL6(verb()," CombiGrid::eval , globalCoords[globalI]:" << globalCoords[domain()->localToGlobalIndex(ii)]);
		//FITOB_OUT_LEVEL6(verb()," CombiGrid::eval , localCoords[ii]:" << localCoords[ii]);
		for (; ;){
		  middle = (minIndex[ii] + maxIndex[ii])/2;
		  //FITOB_OUT_LEVEL6(verb()," CombiGrid::eval , ii: " << ii << ", middle:" << middle <<
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

		//FITOB_OUT_LEVEL6(verb()," CombiGrid::eval , ii: " << ii << ", level[ii]:" << domain()->getAxisLevel(globalIndex) );
		//FITOB_OUT_LEVEL6(verb()," CombiGrid::eval , ii: " << ii << ", localCoords[ii]:" << localCoords[ii]);
		//FITOB_OUT_LEVEL6(verb()," CombiGrid::eval , ii: " << ii << ", intersect[ii]:" << intersect[ii]);
		//FITOB_OUT_LEVEL6(verb()," CombiGrid::eval , ii: " << ii << ", unitCoords[ii]:" << unitCoords[ii]);
		//FITOB_OUT_LEVEL6(verb()," CombiGrid::eval , ii: " << ii << ", axisGrading[minIndex[ii]]:" << axisGrading[minIndex[ii]]);
		//FITOB_OUT_LEVEL6(verb()," CombiGrid::eval , ii: " << ii << ", axisGrading[maxIndex[ii]]:" << axisGrading[maxIndex[ii]]);
		//FITOB_OUT_LEVEL6(verb()," CombiGrid::eval , ii: " << ii << ", minIndex[ii]:" << minIndex[ii]);
		//FITOB_OUT_LEVEL6(verb()," CombiGrid::eval , ii: " << ii << ", maxIndex[ii]:" << maxIndex[ii]);
	}

    // eval the combi grid on the unit coordinates
	double eval_val = combinationGrid_->eval( unitCoords );
	FITOB_OUT_LEVEL5(verb()," CombiGrid::eval RET_VAL:" << eval_val );
	return eval_val;
}

void CombiGrid::eval(const std::vector<DVector>& globalCoordonates , DVector& resVector) const {
	   for (unsigned int ii = 0 ; ii < globalCoordonates.size() ; ii++){
		   // we just call the simple version of this function
		   // todo: this could be optimized by copy paste the body of the function above
		   resVector[ii] = this->eval(globalCoordonates[ii]);
	   }
}

void CombiGrid::setValues(const Evaluable* func , const DVector& globalCoords){

	// get the values from SGpp::FullGrids
	deCompose();

	FITOB_OUT_LEVEL3(verb()," CombiGrid::setValues set values");

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
	    FITOB_OUT_LEVEL3(verb()," CombiGrid::setValues FullGrid nr. i:" << i );
	    fitobFG.setValues( func , tmpGlCoo );
    }
#if defined(FITOB_OPENMP)
}
#endif

	// we do not reCompose, because all the values are in the grid
	FITOB_OUT_LEVEL3(verb()," CombiGrid::setValues END");
}

void CombiGrid::applyConstraints(const OperatorSequence* constraintOpSeq_ ,
		                      const DVector& globalCoords ) {
	// get the values from SGpp::FullGrids
	deCompose();

	FITOB_OUT_LEVEL3(verb()," CombiGrid::applyConstraints apply constraints ");

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
	      FITOB_OUT_LEVEL3(verb()," CombiGrid::applyConstraints apply constraints for FullGrid nr. i:" << i );
	      fitobFG.applyConstraints( constraintOpSeq_ , tmpGlCoo );
      }
#if defined(FITOB_OPENMP)
}
#endif

	// we do not reCompose, because all the values are in the grid
	FITOB_OUT_LEVEL3(verb()," CombiGrid::applyConstraints END");
}

void CombiGrid::deCompose(){

	FITOB_OUT_LEVEL3(verb()," CombiGrid::deCompose ");

	// create the full grids if nor existing before
	if (fullGridLevels_.size() < 1)
	{
		// resize the levels
		FITOB_OUT_LEVEL3(verb()," CombiGrid::deCompose() , CREATE COMBI GRID ");
		fullGridLevels_.resize( (int)combinationGrid_->getNrFullGrid() );
		// create the grid
		nrFullGrids_ = combinationGrid_->getNrFullGrid();
		for ( int i=0 ; i < nrFullGrids_ ; i++ )
		{
		   	IVector specialLevels( domain()->nrRealAxis() );
		   	combigrid::FullGridD *fg = (combinationGrid_->getFullGrid(i));
		   	// copy the levels
		   	for (int jj = 0 ; jj < domain()->nrRealAxis() ; jj++ ){
		   		FITOB_OUT_LEVEL3(verb()," CombiGrid::deCompose , create levels jj:"<<jj<<" fg->getLevel()[jj]: " << fg->getLevels()[jj] );
		   		specialLevels[jj] = fg->getLevels()[jj];
		   	}
		   	// create fullgrid
		   	FITOB_OUT_LEVEL3(verb()," CombiGrid::deCompose() , create FULLGRID nrUnkowns:" << fg->getNrElements());

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

	// transform the SGpp full grids into fitob::FullGrids
	FITOB_OUT_LEVEL3(verb()," CombiGrid::deCompose() , START DE COMBI ");
	for ( int i=0 ; i < combinationGrid_->getNrFullGrid() ; i++ )
	{
		// sg::FullGrid is also the base class for the full grid without boundary points
		combigrid::FullGridD *fg = (combinationGrid_->getFullGrid(i));
	   	int  m = fg->getNrElements();
	   	FullGridBase& fitobFG = fullGridContainer_[i];

	   	// test if the size is not equal
	   	FITOB_ERROR_TEST( (int)fitobFG.unknVect().size() == m , " CombiGrid::deCompose , fitobFG.unknVect().size()" <<
	   			fitobFG.unknVect().size() << " , m:" << m );

	   	// copy the vector values
	   	FITOB_OUT_LEVEL3(verb()," CombiGrid::deCompose() START i:" << i << " m:"<<m);
	   	//for ( unsigned int j = 0 ; j < m ; j++) {
	   		//fitobFG.unknVect()[j] = (*fg)[j];
	   		//if (verb() > 3) std::cout << ", " << fitobFG.unknVect()[j];
	   	//}
	   	//if (verb() > 3) std::cout << std::endl;
	   	FITOB_OUT_LEVEL3(verb()," CombiGrid::deCompose() END i:" << i << " m:"<<m);
	}
	FITOB_OUT_LEVEL3(verb()," CombiGrid::deCompose() , END DE COMBI ");
}



void CombiGrid::reCompose(double &errorIndicator, bool change_coef){

	FITOB_OUT_LEVEL3(verb()," CombiGrid::reCompose ");

	// create the combi object if necessary
	FITOB_ERROR_TEST( combinationGrid_ != 0 , " CombiGrid::reCompose , combinationGrid_ is NULL , ");
	FITOB_ERROR_TEST( fullGridLevels_.size() > 0 , " CombiGrid::reCompose , fullGridLevels_.size() > 0 ");

	// transform the SGpp full grids into fitob::FullGrids
	FITOB_OUT_LEVEL3(verb()," CombiGrid::deCompose() , START COMBI ");

#ifdef SGPP_DIRECT_SOLVER
	// create
	if (gridstorageSGpp_ == 0){
    	// create the grid storage
		gridstorageSGpp_ = combinationGrid_->createSGppGridStorage();
    	alphas_ = new DataVector(gridstorageSGpp_->size());
    	minAlpha_ = new DataVector(gridstorageSGpp_->size());
    	maxAlpha_ = new DataVector(gridstorageSGpp_->size());
    	//for (size_t i = 0 ; i < gridstorageSGpp_->size() ; i++){
    	//	alphas_[0] = 0.0; minAlpha_[0] = 0.0; maxAlpha_[0] = 0.0;
    	//}
	}
#endif

#ifdef SGPP_DIRECT_SOLVER
	// recompose the grids
	combinationGrid_->reCompose( gridstorageSGpp_ , alphas_ , minAlpha_ , maxAlpha_ );

	// calculate the error indicators (avoid division by zero)
	// inf( (maxAlpha[i]-minAlpha[i])) / norm_L2(alphas_)
	double tmp, summError = 0.0 , L2NormAlphas = 0 , maxAplha = 0.0;
	double maxPointWiseDiff = 0.0;
	for (size_t i = 0 ; i < alphas_->getSize() ; i++){
		L2NormAlphas = L2NormAlphas + ((*alphas_)[i] * (*alphas_)[i]);
		tmp = ((*maxAlpha_)[i] - (*minAlpha_)[i]); // this is always positive
		maxPointWiseDiff = ( maxPointWiseDiff < tmp ) ? tmp : maxPointWiseDiff;
		summError = summError + tmp*tmp;
		tmp = fabs((*alphas_)[i]);
		maxAplha = (maxAplha < tmp) ? tmp : maxAplha;
	}
	// the error indicator is calculated as specified in the formula,
	if (maxAplha < 1e-14){
		// avoid division by zero
		errorIndicator = 0.0;
	}else{
		errorIndicator = maxPointWiseDiff / maxAplha;
		//errorIndicator = maxPointWiseDiff / (sqrt( L2NormAlphas )/(double)alphas_->getSize());
		//errorIndicator = (sqrt( summError )/(double)alphas_->getSize()) / (sqrt( L2NormAlphas )/(double)alphas_->getSize());
	}

	// decompose the SGpp into the full grids
	combinationGrid_->deCompose( gridstorageSGpp_ , alphas_ );

	for ( int i=0 ; i < combinationGrid_->getNrFullGrid() ; i++ )
	{
		// sg::FullGrid is also the base class for the full grid without boundary points
		combigrid::FullGridD *fg = (combinationGrid_->getFullGrid(i));
		int m = fg->getNrElements();
	   	FullGridBase& fitobFG = fullGridContainer_[i];

	   	// test if the size is not equal
	   	FITOB_ERROR_TEST( (int)fitobFG.unknVect().size() == m , " CombiGrid::deCompose , fitobFG.unknVect().size()" <<
	   			fitobFG.unknVect().size() << " , m:" << m );
	}
#else
	FITOB_OUT_LEVEL3(verb()," CombiGrid::reCompose  NO SGPP active ");
	errorIndicator = 0.0; // nothing to do
#endif
	FITOB_OUT_LEVEL3(verb()," CombiGrid::deCompose() , END COMBI ");
}

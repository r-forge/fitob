/*
 * FitobSmootherBase.cpp
 *
 *  Created on: Jul 6, 2010
 *      Author: benk
 */

#include "FitobSmootherBase.hpp"
#include "src/pdesolver/FitobGaussSeidel.hpp"
#include "src/pdesolver/FitobGaussSeidel_WB.hpp"

using namespace fitob;
using namespace std;

SmootherBase::SmootherBase( const XMLConfiguration* config, MultigridFGBase *grid ,
		 const ModelCollection* models , const FitobCalculator* fitobcalculator
		 ) : configuration_(config) , grid_(grid) , models_(models) ,
		 fitobcalculator_(fitobcalculator) , microTimestep_(0.0) {

	// the number of factor
	nrFactors_ = models_->nrFactors();

    // store the mesh dimension
	gdim_ = grid_->dim();

	// test weather we have any correlations, if not we might save some computational time
	hasCorrelations_ = false;
	for (int i=1 ; i < nrFactors_ ; i++){
		for (int j = 0 ; j < i ; j++)
			if ( fabs(models_->corr(j,i)) > 1e-5 ){
				hasCorrelations_ = true;
			}
	}
}

boost::shared_ptr<SmootherBase> SmootherBase::createSmoother(const XMLConfiguration* config ,
		MultigridFGBase *grid , const ModelCollection* models ,
		const FitobCalculator* fitobcalculator ){

	// later on we might have different smoothers
   switch (grid->getMultigrigFGType()){
   case (MG_FG_WITH_BOUNDARY):
			//return a smoother which works on full grids with boundary points
			return( boost::shared_ptr<SmootherBase>( new GaussSeidel( config , grid , models , fitobcalculator ) ) );
   case (MG_FG_WITHOUT_BOUNDARY):
			//return a smoother which works on full grids without boundary points
			return( boost::shared_ptr<SmootherBase>( new GaussSeidel_WB( config , grid , models , fitobcalculator ) ) );
   default :
	   FITOB_ERROR_EXIT("SmootherBase::createSmoother , no MultigridFG could be created!");
   }
   // return NULL if somehow gets here
   return( boost::shared_ptr<SmootherBase>((GaussSeidel*)NULL));
}

SmootherBase* SmootherBase::createSmootherPointer(const XMLConfiguration* config ,
		MultigridFGBase *grid , const ModelCollection* models ,
		const FitobCalculator* fitobcalculator ){

	// later on we might have different smoothers
   switch (grid->getMultigrigFGType()){
   case (MG_FG_WITH_BOUNDARY):
			//return a smoother which works on full grids with boundary points
			return( new GaussSeidel( config , grid , models , fitobcalculator ) );
   case (MG_FG_WITHOUT_BOUNDARY):
			//return a smoother which works on full grids without boundary points
			return( new GaussSeidel_WB( config , grid , models , fitobcalculator ) );
   default :
			FITOB_ERROR_EXIT("SmootherBase::createSmootherPointer , no Multigrid FG could be created!");
   }
   // return NULL if somehow gets here
   return( (GaussSeidel*)NULL );
}

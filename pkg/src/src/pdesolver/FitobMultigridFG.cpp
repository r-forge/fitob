/*
 * FitobMultigridFG.cpp
 *
 *  Created on: Jul 12, 2010
 *      Author: benk
 */

#include "FitobMultigridFG.hpp"
#include "src/operators/FitobOperatorSequence.hpp"
#include <boost/lexical_cast.hpp>

using namespace fitob;
using namespace std;

MultigridFG::MultigridFG(const FullGridBase* fullgrid) : MultigridFGBase(fullgrid){

	setVerb(0);

	//test if full grid is with boundary points!
	FITOB_ERROR_TEST( fullgrid->getGridType() == GRID_WITH_BOUNDARY , " MultigridFG::MultigridFG , fullgrid must be with boundary points " );

// in this case we create one grid similar to the full grid input
	dim_ = fullgrid->dim();
	innerPoints_ = 1;
	levels_.resize(dim_);
	nrPoints_.resize(dim_);
	offset_.resize(dim_);
	axisScaling_.resize(dim_);
	FITOB_OUT_LEVEL3( verb() , " MultigridFG::MultigridFG , fullgrid->dim():" << fullgrid->dim());
	for (int i = 0 ; i < dim_ ; i++){
		levels_[i] = fullgrid->getAxisLevel(i);
		nrPoints_[i] = fullgrid->nrPointsPerAxis(i);
		innerPoints_ = innerPoints_*(nrPoints_[i]-2);
		offset_[i] = fullgrid->offs(i);
		axisScaling_[i].resize(fullgrid->getScalingAxis(i).size());
		FITOB_OUT_LEVEL3( verb() , " MultigridFG::MultigridFG , i:" << i << " axisScaling_[i].size():"<<axisScaling_[i].size());
		for (unsigned int j=0 ; j < fullgrid->getScalingAxis(i).size() ; j++ ){
			axisScaling_[i][j] = fullgrid->getScalingAxis(i)[j];
		}
	}
	rhs_.resize(fullgrid->totalPoints());
	u_.resize(fullgrid->totalPoints());
	d_.resize(fullgrid->totalPoints());
	residuum_.resize(fullgrid->totalPoints());
// do with OpenMP the initialization, so that by touch first the memory is allocated locally for faster access
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp parallel
{
#endif
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp for schedule(static)
#endif
	for (int i=0 ; i < fullgrid->totalPoints() ; i++){
		//FITOB_OUT_LEVEL3( verb() , " u_["<< i<<"]=" << u_[i] << " out of:" << fullgrid->totalPoints());
		//FITOB_OUT_LEVEL3( verb() , " size " << fullgrid->unknVect().size());
		u_[i] = fullgrid->val(i);
		//FITOB_OUT_LEVEL3( verb() , " AFT ");
		rhs_[i] = 0.0; d_[i] = 0.0;
		residuum_[i] = 0.0;
	}
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
}
#endif

	FITOB_OUT_LEVEL3( verb() , " MultigridFG::MultigridFG , END");
	domain_ = fullgrid->domain();
}

MultigridFG::MultigridFG(const MultigridFG *mgfg , bool coarse) : MultigridFGBase( mgfg , coarse ) {

	setVerb(0);

	//test if full grid is with boundary points!
	FITOB_ERROR_TEST( mgfg->getMultigrigFGType() == MG_FG_WITH_BOUNDARY , " MultigridFG::MultigridFG , MultigridFG must be with boundary points " );

	// store the domain
	domain_ = mgfg->domain_;
	innerPoints_ = 1;

	if (coarse)
	{
		// here is the technique to coarse for the multigrid solution
		// the recommended strategy is to coarse the finest axis
		// set the unknowns to the unknown of the fine mesh
		dim_ = mgfg->dim_;
		levels_.resize(dim_);
		nrPoints_.resize(dim_);
		offset_.resize(dim_);
		axisScaling_.resize(dim_);
		int maxlevel = 0;
		// get the maximum level
		for (int i = 0 ; i < dim_ ; i++){
			maxlevel = (maxlevel < mgfg->levels_[i]) ? (mgfg->levels_[i]): (maxlevel);
		}
		// copy the data and make refinement
		for (int i = 0 ; i < dim_ ; i++){
			// --- HERE WE HAVE THE COARSENING STRATEGY FOR THE MULTIFRID METHOD ---
			// we try to coarse only the axis with the highest level, in this way we might have
			// more points then coarsening all axis, but so we have to hope a more robust MG
			// (similar to line coarsening)
			if ( maxlevel == mgfg->levels_[i]){
			    // MAKE coarsening on this axis
				levels_[i] = mgfg->levels_[i]-1;
				nrPoints_[i] = ((mgfg->nrPoints_[i]+1)/2);
				innerPoints_ = innerPoints_*(nrPoints_[i]-2);
				axisScaling_[i].resize((mgfg->axisScaling_[i].size()+1)/2);
				for (unsigned int j=0 ; j < axisScaling_[i].size() ; j++ )
					axisScaling_[i][j] = mgfg->axisScaling_[i][j*2];
				FITOB_OUT_LEVEL3( verb() , " MultigridFG::MultigridFG COARSE maxlevel:"<<maxlevel<<" nrPoints_[i]"<<nrPoints_[i]);
			}else{
			    // NO coarsening on this axis
				levels_[i] = mgfg->levels_[i];
				nrPoints_[i] = mgfg->nrPoints_[i];
				innerPoints_ = innerPoints_*(nrPoints_[i]-2);
				axisScaling_[i].resize(mgfg->axisScaling_[i].size());
				for (unsigned int j=0 ; j < mgfg->axisScaling_[i].size() ; j++ )
					axisScaling_[i][j] = mgfg->axisScaling_[i][j];
				FITOB_OUT_LEVEL3( verb() , " MultigridFG::MultigridFG NO COARSE maxlevel:"<<maxlevel<<" nrPoints_[i]"<<nrPoints_[i]);
			}
		}
		// set offsets, according to the refinement
		int linOffs = 1;
		for (int i = 0 ; i < dim_ ; i++){
			offset_[i] = linOffs;
			linOffs = linOffs*nrPoints_[i];
		}
		rhs_.resize(linOffs);
		u_.resize(linOffs);
		d_.resize(linOffs);
		residuum_.resize(linOffs);
		// do with OpenMP the initialization, so that by touch first the memory is allocated locally for faster access
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp parallel
{
#endif
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp for schedule(static)
#endif
		for (int i = 0; i < linOffs ; i++){
			rhs_[i] = u_[i] = d_[i] = residuum_[i] = 0.0;
		}
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
}
#endif
	}
	else
	{
	// just clone the mesh , except the diag rhs and res vectors
		dim_ = mgfg->dim_;
		levels_.resize(dim_);
		nrPoints_.resize(dim_);
		offset_.resize(dim_);
		axisScaling_.resize(dim_);
		for (int i = 0 ; i < dim_ ; i++){
			levels_[i] = mgfg->levels_[i];
			nrPoints_[i] = mgfg->nrPoints_[i];
			innerPoints_ = innerPoints_*(nrPoints_[i]-2);
			offset_[i] = mgfg->offset_[i];
			axisScaling_[i].resize(mgfg->axisScaling_[i].size());
			for (unsigned int j=0 ; j < mgfg->axisScaling_[i].size() ; j++ )
				axisScaling_[i][j] = mgfg->axisScaling_[i][j];
		}
		rhs_.resize(mgfg->rhs_.size());
		u_.resize(mgfg->u_.size());
		d_.resize(mgfg->d_.size());
		residuum_.resize(mgfg->residuum_.size());
	// do with OpenMP the initialization, so that by touch first the memory is allocated locally for faster access
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp parallel
{
#endif
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp for schedule(static)
#endif
		for ( unsigned int i=0 ; i < mgfg->u_.size() ; i++){
			u_[i] = mgfg->u_[i]; rhs_[i] = 0.0; d_[i] = 0.0;
			residuum_[i] = 0.0;
		}
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
}
#endif
	}
}

void MultigridFG::applyConstraints(
		const OperatorSequence* constraintOpSeq_ ,  const DVector& globalCoords ){
#if defined(FITOB_OPENMP)
#pragma omp parallel
{
#endif
	DVector globalc = globalCoords;

	// apply the restriction for each point
	int tmp_I , linearIndex;
	IVector axisIndex(dim_);
#if defined(FITOB_OPENMP)
#pragma omp for schedule(static)
#endif
	for (unsigned int ind=0 ; ind < u_.size() ; ind++){
		   tmp_I = (int)ind;
		   linearIndex = 0;
		   for (int i = dim_-1 ; i >= 0 ; i--){
			   axisIndex[i] = tmp_I / (offset_[i]);
			   tmp_I = tmp_I % offset_[i];
			   linearIndex = linearIndex + (axisIndex[i])*offset_[i];
		       // set the global coordinates
			   //FITOB_OUT_LEVEL3(5," GaussSeidl::updateStage: ind:" << ind << " , tmp_I:" << tmp_I <<
				//	   ", nrPoints_[i]:" << nrPoints_[i] << " , u_.size():" << u_.size());
			   //FITOB_OUT_LEVEL3(5," GaussSeidl::updateStage: i:"<<i<<" domain_->localToGlobalIndex(i):" << domain_->localToGlobalIndex(i)
				//	   << " , axisIndex[i]:"<<axisIndex[i]);
			   //FITOB_OUT_LEVEL3(5," GaussSeidl::updateStage: globalc.size():"<<globalc.size()<<" axisScaling_[i].size():" << axisScaling_[i].size());
			   globalc[ domain_->localToGlobalIndex(i) ] = axisScaling_[i][axisIndex[i]];
		   }

		   // call the apply constraint for one global set of coordinates
		   globalc[1] = u_[ind];
		   constraintOpSeq_->applyConstraints(globalc);
		   u_[ind] = globalc[1];
	}
#if defined(FITOB_OPENMP)
}
#endif
}

void MultigridFG::plotMAT_grad(const string& filename , DVector& vect , int count ) const {

   const string completeName = filename + boost::lexical_cast<std::string>(count)+".m";
   DVector globalCoord = domain_->getAverage();
   DVector result(0);

   // --- PLOTTING --- 1D
   if (axisScaling_.size() == 1)
   {
      result.resize( axisScaling_[0].size() );
      // loop and evaluate points
      for ( unsigned int ii = 0 ; ii < axisScaling_[0].size() ; ii++){
    	  result[ii] = vect[ii];
      }

       // writing file
       ofstream myfile;
       myfile.open(completeName.c_str());
       myfile << "X = [ " << axisScaling_[0][0];
       for ( unsigned int ii = 1 ; ii < axisScaling_[0].size() ; ii++){
    	   myfile << " , " << axisScaling_[0][ii];
       }
       myfile << "]; \n ";
       myfile << "res = [ " << result[0];
       for ( unsigned int ii = 1 ; ii < axisScaling_[0].size() ; ii++) {
    	   myfile << " , " << result[ii];
       }
       myfile << "]; \n ";
       myfile << " plot(X,res); \n ";
       myfile.close();

   }

   // -- PLOTTING ---  >= 2D
   if (axisScaling_.size() >= 2)
   {
	  result.resize( axisScaling_[0].size() * axisScaling_[1].size() );
	  // loop and evaluate points
      for ( unsigned int ii = 0 ; ii < axisScaling_[0].size() ; ii++){
    	  for ( unsigned int jj = 0 ; jj < axisScaling_[1].size() ; jj++){
    		  // todo: this will work only in 2D
    	     result[ii*axisScaling_[1].size() + jj] = vect[jj*axisScaling_[0].size() + ii];
    	  }
      }

      ofstream myfile;
      myfile.open(completeName.c_str());
      myfile << "X = [ " << axisScaling_[0][0];
      for ( unsigned int ii = 1 ; ii < axisScaling_[0].size() ; ii++){
   	     myfile << " , " << axisScaling_[0][ii];
      }
      myfile << "]; \n ";
      myfile << "Y = [ " << axisScaling_[1][0];
      for ( unsigned int ii = 1 ; ii < axisScaling_[1].size() ; ii++){
   	     myfile << " , " << axisScaling_[1][ii];
      }
      myfile << "]; \n ";
      myfile << "res = [ " << result[0];
      for ( unsigned int ii = 1 ; ii < axisScaling_[0].size()*axisScaling_[1].size() ; ii++) {
    	  if ( (ii % axisScaling_[1].size()) == 0){
    		  myfile << " ; " << result[ii];
    	  }
    	  else{
   	         myfile << " , " << result[ii];
    	  }
      }
      myfile << "]; \n ";
      myfile << "[x,y]=meshgrid(Y,X);\n";
      myfile << " surf(x,y,res); \n ";
      myfile.close();

   }
}

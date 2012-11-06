/*
 * FitobProlongationRestriction.cpp
 *
 *  Created on: Jul 6, 2010
 *      Author: benk
 */

#include "FitobProlongationRestriction.hpp"

using namespace fitob;
using namespace std;

ProlongationRestriction::ProlongationRestriction() {
   // nothing to do here
}


void ProlongationRestriction::makeLinearProlongation( MultigridFGBase* gFine , DVector& vFine ,
		                            const MultigridFGBase* gCoarse , const DVector& vCoarse ,
		                            double coef , double coefOLD ){
	// based on the grid type we have to do different types of prolongations
	switch (gFine->getMultigrigFGType()){
		case (MG_FG_WITH_BOUNDARY):{
			makeLinearProlongation_withBoundary( gFine , vFine ,gCoarse , vCoarse , coef , coefOLD);
			break;
		}
		case (MG_FG_WITHOUT_BOUNDARY):{
			makeLinearProlongation_withoutBoundary( gFine , vFine ,gCoarse , vCoarse , coef , coefOLD);
			break;
		}
	}
}

void ProlongationRestriction::makeDirectRestriction( const MultigridFGBase* gFine , const DVector& vFine ,
		                           MultigridFGBase* gCoarse , DVector& vCoarse ,
		                           double coef , double coefOLD ){
	// based on the grid type we have to do different types of restrictions
	switch (gFine->getMultigrigFGType()){
		case (MG_FG_WITH_BOUNDARY):{
			makeDirectRestriction_withBoundary( gFine , vFine ,gCoarse , vCoarse , coef , coefOLD);
			break;
		}
		case (MG_FG_WITHOUT_BOUNDARY):{
			makeDirectRestriction_withoutBoundary( gFine , vFine ,gCoarse , vCoarse , coef , coefOLD);
			break;
		}
	}
}


void ProlongationRestriction::makeLinearProlongation_withBoundary(
		MultigridFGBase* gFine , DVector& vFine ,
		const MultigridFGBase* gCoarse , const DVector& vCoarse ,
		double coef , double coefOLD){
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp parallel
{
#endif
	// variable declaration
	//const int verb = 0;
	const int dim = gFine->dim();
	const int nrFinePoints = gFine->u().size();
	const int nrCellPoints = fitob::powerTwo[dim];
	int tmp_I , linearIndexC , tmp_I2 , axis , i ;
	double val = 0.0 , tmp_D;
	int tOffs = 0;
    IVector axisIndexC(dim);
    IVector axisIndexF(dim);
    IVector mulFactors(dim,1);
    DVector axisDiv(dim);

    // first see which axis index (on the fine mesh) have to be multiplied by two
    for (i = 0 ; i < dim ; i++){
    	mulFactors[i] = (gFine->nrPoints(i)+1)/gCoarse->nrPoints(i);
    }
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp for schedule(static)
#endif
    // here we iterate on all fine grid points
    for (int pInd = 0 ; pInd < nrFinePoints ; pInd++)
    {
       tmp_I = pInd;
       //       - locate the fine point on the coarse mesh (fine the N-dim cell)
       //       - make interpolation (eval of the basis function) on the given coord (caz stetching this is not in the middle)
       //            - take the min per axis
       //       - set the value of the fine mesh point

       // calc fine grid point index
       tmp_I = pInd;
       linearIndexC = 0;
	   for (i = dim-1 ; i >= 0 ; i--){
		   axisIndexF[i] = tmp_I / (gFine->offset()[i]);
		   tmp_I = tmp_I % gFine->offset()[i];
		   // locate the lower corner point of the coarse grid // todo: optimize this
		   if (gCoarse->level(i) == gFine->level(i)){
			   // NO LEVEL reduction
			   axisIndexC[i] = (axisIndexF[i] >= (gCoarse->nrPoints(i)-1)) ? (gCoarse->nrPoints(i)-2) : (axisIndexF[i]) ;
		   } else {
			   // LEVEL reduction
			   axisIndexC[i] = (axisIndexF[i]/2 >= (gCoarse->nrPoints(i)-1)) ? (gCoarse->nrPoints(i)-2) : (axisIndexF[i]/2) ;
		   }
		   linearIndexC = linearIndexC + axisIndexC[i]*gCoarse->offset()[i];
		   //FITOB_OUT_LEVEL3( verb , "makeLinearProlongation LOC i:" << i << ", axisIndexF[i]:" << axisIndexF[i]
		   //                 << ", axisIndexC[i]:" << axisIndexC[i] << " Lc:" << gCoarse->level(i) << " Lf:" << gFine->level(i));
		   //FITOB_OUT_LEVEL3( verb , "makeLinearProlongation LOC i:" << i << " sizeC:" << gCoarse->axisScaling(i).size() <<
			//	                     " sizeF:" << gFine->axisScaling(i).size() );
	   }

	   // get the intersection per axis
	   for (i = 0 ; i < dim ; i++){
		   axisDiv[i] = gFine->axisScaling(i)[axisIndexF[i]] - gCoarse->axisScaling(i)[axisIndexC[i]];
		   //FITOB_OUT_LEVEL3( verb , "makeLinearProlongation INT i:" << i << " F:" << gFine->axisScaling(i)[axisIndexF[i]]
		     //          << " C1:" << gCoarse->axisScaling(i)[axisIndexC[i]] << " C2:" << gCoarse->axisScaling(i)[axisIndexC[i]+1]);
		   axisDiv[i] = axisDiv[i]/(gCoarse->axisScaling(i)[axisIndexC[i]+1] - gCoarse->axisScaling(i)[axisIndexC[i]]);
		   // this value must be between 0.0 and 1.0 !!!
		   //FITOB_ERROR_TEST( (axisDiv[i] > -1e-10) && (axisDiv[i] < 1.0+1e-10),
			//	   " ProlongationRestriction::makeLinearProlongation COND: (axisDiv[i] > -1e-8) && (axisDiv[i] < 1.0+1e-8) : " << axisDiv[i]);
	   }

	   // loop over 2^D points and calculate the value of the interpolation
	   val = 0.0 ; tOffs = 0;

	   //todo: this part can be done much smarter, not to evaluate the basis function, pre-compute the possibilities
	   for (i = 0 ; i < nrCellPoints ; i++){
		   tmp_D = 1.0;
		   tOffs = 0;
		   // calculate the basis function value in N-dim
		   for (axis = 0 ; axis < dim ; axis++){
			   // calculate the index on this axis
			   tmp_I2 = (i / powerTwo[axis]) % 2;
			   // N-linear basis function
			   tmp_D =  tmp_D*((double)(tmp_I2)*axisDiv[axis] + (double)(1-tmp_I2)*(1.0-axisDiv[axis]));
			   tOffs = tOffs + tmp_I2*gCoarse->offset()[axis];
			   //FITOB_OUT_LEVEL3( verb , " Calc basis function, axis:"<< axis << " , tmp_I2:" << tmp_I2 << ",tmp_D:"<<tmp_D);
		   }
		   // add the contribution from this point
		   //FITOB_OUT_LEVEL3( verb , "pInd:"<<pInd<<", basis value:" << tmp_D << ", linearIndexC:" << linearIndexC << ", tOffs:" << tOffs);
		   val = val + tmp_D*vCoarse[ linearIndexC + tOffs ];
	   }//
	   vFine[pInd] = coefOLD*vFine[pInd] + coef*val;
       //FITOB_OUT_LEVEL3( verb , " ProlongationRestriction::makeLinearProlongation , indexFine:"
    	//	   << pInd << " , val:"<<val);
    }// end loop over fine grid points
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
}
#endif
}


void ProlongationRestriction::makeDirectRestriction_withBoundary(
		const MultigridFGBase* gFine , const DVector& vFine ,
		MultigridFGBase* gCoarse , DVector& vCoarse ,
		double coef , double coefOLD){
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp parallel
{
#endif
	// variable declaration
	//const int verb = 0;
	const int dim = gFine->dim();
	const int nrCoarsePoints = gCoarse->u().size();
	int tmp_I , linearIndexC , linearIndexF , i;
    IVector axisIndex(dim);
    IVector mulFactors(dim,1);

    // first see which axis index (on the fine mesh) have to be multiplied by two
    for (int i = 0 ; i < dim ; i++){
    	mulFactors[i] = (gFine->nrPoints(i)+1)/gCoarse->nrPoints(i);
        //FITOB_OUT_LEVEL3( verb , " ProlongationRestriction::makeDirectRestriction , mulFactors[i]:" << mulFactors[i] );
    }

#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp for schedule(static)
#endif
    // here we iterate on all coarse grid points
    for (int pInd = 0 ; pInd < nrCoarsePoints ; pInd++){
       tmp_I = pInd;
       // calc coarse grid point index
       linearIndexC = 0;
	   linearIndexF = 0;
	   for (i = dim-1 ; i >= 0 ; i--){
		   axisIndex[i] = tmp_I / (gCoarse->offset()[i]);
		   tmp_I = tmp_I % gCoarse->offset()[i];
		   // calculate the corresponding fine grid point index
		   linearIndexF = linearIndexF + (mulFactors[i]*axisIndex[i])*gFine->offset()[i];
		   //linearIndexC = linearIndexC + (axisIndex[i])*gCoarse->offset()[i];
	   }

	   // make the direct restriction by adding the value with a given coefficient
	   linearIndexC = pInd;
	   vCoarse[linearIndexC] = coefOLD*vCoarse[linearIndexC] + coef*vFine[linearIndexF];
       //FITOB_OUT_LEVEL3( verb , " ProlongationRestriction::makeDirectRestriction , linearIndexC:"
    	//	   << linearIndexC << " , linearIndexF:"<<linearIndexF);
    }// end loop over coarse grid points
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
}
#endif
}


void ProlongationRestriction::makeLinearProlongation_withoutBoundary(
		                           MultigridFGBase* gFine , DVector& vFine ,
		                           const MultigridFGBase* gCoarse , const DVector& vCoarse ,
		                           double coef , double coefOLD ){
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp parallel
{
#endif
	// variable declaration
	//const int verb = 0;
	const int dim = gFine->dim();
	const int nrFinePoints = gFine->u().size();
	const int nrCellPoints = fitob::powerTwo[dim];
	int tmp_I , linearIndexC , tmp_I2 , axis , i;
	double val = 0.0 , tmp_D;
	int tOffs = 0;
    IVector axisIndexC(dim);
    IVector axisIndexF(dim);
    IVector mulFactors(dim,1);
    DVector axisDiv(dim);

    // first see which axis index (on the fine mesh) have to be multiplied by two
    for (i = 0 ; i < dim ; i++){
    	mulFactors[i] = (gFine->nrPoints(i)+1)/gCoarse->nrPoints(i);
    }

#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp for schedule(static)
#endif
    // here we iterate on all fine grid points
    for (int pInd = 0 ; pInd < nrFinePoints ; pInd++)
    {
       tmp_I = pInd;
       //       - locate the fine point on the coarse mesh (fine the N-dim cell)
       //       - make interpolation (eval of the basis function) on the given coord (caz stetching this is not in the middle)
       //            - take the min per axis
       //       - set the value of the fine mesh point

       // calc fine grid point index
       tmp_I = pInd;
       linearIndexC = 0;
	   for (i = dim-1 ; i >= 0 ; i--)
	   {
		   axisIndexF[i] = tmp_I / (gFine->offset()[i]);
		   tmp_I = tmp_I % gFine->offset()[i];
		   // locate the lower corner point of the coarse grid
		   if (gCoarse->level(i) == gFine->level(i)){
			   // NO LEVEL reduction
			   axisIndexC[i] = (axisIndexF[i] >= (gCoarse->nrPoints(i)-1)) ? (gCoarse->nrPoints(i)-2) : (axisIndexF[i]) ;
		   } else {
			   // LEVEL reduction , we subtract 1 because of the boundary
			   axisIndexC[i] = (axisIndexF[i]/2 >= (gCoarse->nrPoints(i)-1)) ? (gCoarse->nrPoints(i)-2) : ( FITOB_IMAX(axisIndexF[i]-1,0)/2) ;
		   }
		   linearIndexC = linearIndexC + axisIndexC[i]*gCoarse->offset()[i];
		   //FITOB_OUT_LEVEL3( verb , "makeLinearProlongation LOC i:" << i << ", axisIndexF[i]:" << axisIndexF[i]
		    //                << ", axisIndexC[i]:" << axisIndexC[i] << " Lc:" << gCoarse->level(i) << " Lf:" << gFine->level(i));
		   //FITOB_OUT_LEVEL3( verb , "makeLinearProlongation LOC i:" << i << " sizeC:" << gCoarse->axisScaling(i).size() <<
			//	                     " sizeF:" << gFine->axisScaling(i).size() );
	   }

	   // get the intersection per axis
	   for (i = 0 ; i < dim ; i++){
		   // the scaling starts at position 1 not at position 0
		   axisDiv[i] = gFine->axisScaling(i)[1+axisIndexF[i]] - gCoarse->axisScaling(i)[1+axisIndexC[i]];
		   //FITOB_OUT_LEVEL3( verb , "makeLinearProlongation INT i:" << i << " F:" << gFine->axisScaling(i)[1+axisIndexF[i]]
		   //            << " C1:" << gCoarse->axisScaling(i)[1+axisIndexC[i]] << " C2:" << gCoarse->axisScaling(i)[1+axisIndexC[i]+1]);
		   axisDiv[i] = axisDiv[i]/(gCoarse->axisScaling(i)[1+axisIndexC[i]+1] - gCoarse->axisScaling(i)[1+axisIndexC[i]]);
		   // since the finer grid boundary points can be outside the cell
		   axisDiv[i] = FITOB_DMAX( FITOB_DMIN(axisDiv[i],1.0) , 0.0);
	   }

	   // loop over 2^D points and calculate the value of the interpolation
	   val = 0.0 ; tOffs = 0;

	   //todo: this part can be done much smarter, not to evaluate the basis function, pre-compute the possibilities
	   for (i = 0 ; i < nrCellPoints ; i++){
		   tmp_D = 1.0;
		   tOffs = 0;
		   // calculate the basis function value in N-dim
		   for (axis = 0 ; axis < dim ; axis++){
			   // calculate the index on this axis
			   tmp_I2 = (i / powerTwo[axis]) % 2;
			   // N-linear basis function
			   tmp_D =  tmp_D*((double)(tmp_I2)*axisDiv[axis] + (double)(1-tmp_I2)*(1.0-axisDiv[axis]));
			   tOffs = tOffs + tmp_I2*gCoarse->offset()[axis];
			   //FITOB_OUT_LEVEL3( verb , " Calc basis function, axis:"<< axis << " , tmp_I2:" << tmp_I2 << ",tmp_D:"<<tmp_D);
		   }
		   // add the contribution from this point
		   //FITOB_OUT_LEVEL3( verb , "pInd:"<<pInd<<", basis value:" << tmp_D << ", linearIndexC:" << linearIndexC << ", tOffs:" << tOffs);
		   val = val + tmp_D*vCoarse[ linearIndexC + tOffs ];
	   }//
	   vFine[pInd] = coefOLD*vFine[pInd] + coef*val;
       //FITOB_OUT_LEVEL3( verb , " ProlongationRestriction::makeLinearProlongation , indexFine:"
    	//	   << pInd << " , val:"<<val);
    }// end loop over fine grid points
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
}
#endif
}

void ProlongationRestriction::makeDirectRestriction_withoutBoundary(
		                           const MultigridFGBase* gFine , const DVector& vFine ,
		                           MultigridFGBase* gCoarse , DVector& vCoarse ,
		                           double coef , double coefOLD ){
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp parallel
{
#endif
	// variable declaration
	//const int verb = 6;
	const int dim = gFine->dim();
	const int nrCoarsePoints = gCoarse->u().size();
	int tmp_I , linearIndexC , linearIndexF , i;
    IVector axisIndex(dim);
    IVector mulFactors(dim,1);
    IVector start_offs(dim,0);

    // first see which axis index (on the fine mesh) have to be multiplied by two
    for (int i = 0 ; i < dim ; i++){
    	mulFactors[i] = (gFine->nrPoints(i)+1)/gCoarse->nrPoints(i); // this is either 1 or 2
    	start_offs[i] = (mulFactors[i]>1) ? 1 : 0;
        //FITOB_OUT_LEVEL3( verb , " ProlongationRestriction::makeDirectRestriction , mulFactors[i]:" << mulFactors[i] );
    }

    //FITOB_OUT_LEVEL3( verb , " nrCoarsePoints:" << nrCoarsePoints );
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp for schedule(static)
#endif
    // here we iterate on all coarse grid points
    for (int pInd = 0 ; pInd < nrCoarsePoints ; pInd++){
       tmp_I = pInd;
       // calc coarse grid point index
	   linearIndexF = 0;
	   for (i = dim-1 ; i >= 0 ; i--){
		   axisIndex[i] = (tmp_I / (gCoarse->offset()[i]));
		   tmp_I = tmp_I % gCoarse->offset()[i];
		   // calculate the corresponding fine grid point index , restriction starts with point 1
		   linearIndexF = linearIndexF + (start_offs[i] + mulFactors[i]*axisIndex[i])*gFine->offset()[i];
		   //FITOB_OUT_LEVEL3( verb , " i:" << i << " axisIndex[i]:" << axisIndex[i] << " , mulFactors[i]:" << mulFactors[i] );
		   //FITOB_OUT_LEVEL3( verb , " i:" << i << " linearIndexF:" << linearIndexF << " , gFine->offset()[i]:" << gFine->offset()[i] );
	   }

	   // make the direct restriction by adding the value with a given coefficient
	   //FITOB_OUT_LEVEL3( verb , " pInd:" << pInd );
	   linearIndexC = pInd;
	   vCoarse[linearIndexC] = coefOLD*vCoarse[linearIndexC] + coef*vFine[linearIndexF];
       //FITOB_OUT_LEVEL3( verb , " ProlongationRestriction::makeDirectRestriction , linearIndexC:"
    	//	   << linearIndexC << " , linearIndexF:"<<linearIndexF);
    }// end loop over coarse grid points
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
}
#endif
}

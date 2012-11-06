/*
 * FideumGaussSeidelWB.cpp
 *
 *  Created on: Jan 3, 2011
 *      Author: benk
 */

#include "FitobGaussSeidel_WB.hpp"

using namespace fitob;
using namespace std;

GaussSeidel_WB::GaussSeidel_WB( const XMLConfiguration* config , MultigridFGBase *grid ,
	     const ModelCollection* models ,
	     const FitobCalculator* fitobcalculator) : SmootherBase(config,grid, models, fitobcalculator){

	setVerb(4);

	isDiagonalCalculated_ = false;
	calculateRHS_ = false;

    factorIndex_to_local.resize(nrFactors_);
    innerOffset.resize(gdim_);

	innerPoints = 1; // number of points where we'll iterate
	// the local index of each factor

    // copy the constant global coordinates
    // todo: update time in the global coordinates (or this should be done somewhere above)
    const ModelCollection* models_tmp = models;

    int linearIndex = 1;
	for (int i = 0 ; i < gdim_ ; i++){
		// calculate the inner offset
		innerOffset[i] = linearIndex;
		linearIndex = linearIndex*(grid->nrPoints(i));
	    // the inner points are
		innerPoints = innerPoints*(grid->nrPoints(i));
		FITOB_OUT_LEVEL3(verb()," GaussSeidel_WB::GaussSeidel_WB: i:"<<i<<" , linearIndex:"<<linearIndex<<" , innerOffset[i]:"
				             << innerOffset[i] << " , innerPoints:"<<innerPoints);
	}
	for (int i = 0 ; i < nrFactors_ ; i++){
	    // get the local index of each factor
		const FactorModelBase& fm = models_tmp->getModel(i);
		int glIndex = fm.getGlobalIndex();
		factorIndex_to_local[i] = grid->domain()->globalToLocalIndex( glIndex );
	}

}

void GaussSeidel_WB::smoothGrid( const DVector& globalCoords) {
	// the function which iterates one over the mesh

	// -- setup phase for the iteration
    DVector globCoord = globalCoords;

	// compute the diagonal if it is not computed
	calculateDiag( globCoord );

	// compute the right hand side if not computed
	calculateRHS( globCoord );

	// - iterate differently through the mesh
	oneSmoothGrid( globCoord );
}


void GaussSeidel_WB::calculateRHS( DVector& globCoord){

	  //int verb = 6;
	  // ====================== calculate right hand side =====================
		// in case of the multigrid only at the highest level we should do this in correction scheme
		// in FAS case on each level we should do this and similar things on each level
		if (!calculateRHS_){
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp parallel
{
#endif
			   int linearIndex;
			   IVector axisIndex(gdim_,0);
			   DVector h0(nrFactors_);
			   DVector h1(nrFactors_);
			   DVector a_x_ij(nrFactors_);
			   DVector a_xx_ij(nrFactors_);
			   DVector gridV( 3 * nrFactors_  , 0.0);
			   IVector tmpMark( 3 * nrFactors_ , 0 );
			   DVector gridV_xy( 9 , 0.0 );
			   IVector tmpMark_xy( 9 , 0.0 );
			   DVector globCoord_tmp = globCoord;

			   // variables to sum up the formula
			   double s_mu,s_sig,s_coor,s_r;

#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp for schedule(static)
#endif
		       for (int pI = 0 ; pI < innerPoints ; pI++){

	               // update the status
		    	   updateStage( pI , globCoord_tmp , linearIndex , axisIndex , gridV ,
		    			   tmpMark , h0 , h1 , a_x_ij,  a_xx_ij );

		    	   // calculate diagonal according to 1.33
		    	   s_mu = 0.0; s_sig = 0.0; s_coor = 0.0; s_r = 0.0;
		    	   // convection
                           #pragma simd reduction(+:s_mu)
		    	   for (int i = 0 ; i < nrFactors_ ; i++){
		    		   s_mu = s_mu + models()->getModel(i).convectionCoef(globCoord_tmp) *
		    				   V_ij_x( i , h0[i] , h1[i] , a_x_ij[i] , gridV );
		    	   }
		    	   // diffusion
                           #pragma simd reduction(+:s_sig)
		    	   for (int i = 0 ; i < nrFactors_ ; i++){
		    		   s_sig = s_sig + models()->getModel(i).diffusionCoef(globCoord_tmp) * models()->getModel(i).diffusionCoef(globCoord_tmp)
		    				   * V_ij_xx( i , h0[i] , h1[i] , a_xx_ij[i] , gridV );
		    	   }
		    	   // correlation
		    	   if ( this->hasCorrelations_){
		    		   for (int i = 0 ; i < nrFactors_ ; i++)
		    		   {
                                           #pragma simd reduction(+:s_coor)
		    			   for (int j = i+1 ; j < nrFactors_ ; j++)
		    			   {
		    				   updateStage_xy( i , j , gridV_xy , gridV , linearIndex , axisIndex , tmpMark_xy);
		    				   s_coor = s_coor + models()->corr(i,j) *
		    						   models()->getModel(i).diffusionCoef(globCoord_tmp) * models()->getModel(j).diffusionCoef(globCoord_tmp)
		    						 * V_ij_xy( i , j , h0[i] , h1[i] , h0[j] , h1[j] , a_x_ij[i] , a_x_ij[j] , gridV_xy );
		    			   }
		    		   }
		    	   }
		    	   // interest rate
		    	   s_r = - models()->r(globCoord_tmp) * grid()->u()[linearIndex];

		    	   // sum up all the terms
		    	   grid()->rhs()[linearIndex] = grid()->u()[linearIndex] + 0.5*microTimestep_*(s_mu + 0.5*s_sig + s_coor + s_r);
		  		   //FITOB_OUT_LEVEL3( verb , "calculateRHS linearIndex:" << linearIndex << " , rhs():" << grid()->rhs()[linearIndex]);
		       }
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
}
#endif
		}
		calculateRHS_ = true;
	  // ====================== end of calculating right hand side =====================
}

double GaussSeidel_WB::calcResiduum( const DVector& globalCoords ){

	double errorL2 = 0.0;
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp parallel reduction(+ : errorL2)
{
#endif
	// -- setup phase for the iteration
    DVector globCoord = globalCoords;

    //int verb = 6;
	// variables to sum up the formula
	double s_mu,s_sig,s_coor,s_r;
	int linearIndex;
	IVector axisIndex(gdim_,0);
	DVector h0(nrFactors_);
	DVector h1(nrFactors_);
	DVector a_x_ij(nrFactors_);
	DVector a_xx_ij(nrFactors_);
	DVector gridV( 3 * nrFactors_  , 0.0);
	IVector tmpMark( 3 * nrFactors_ , 0 );
	DVector gridV_xy( 9 , 0.0 );
	IVector tmpMark_xy( 9 , 0.0 );
	DVector globCoord_tmp = globCoord;

#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp for schedule(static)
#endif
    for (int pI = 0 ; pI < innerPoints ; pI++){

       // update the status
	   updateStage( pI , globCoord_tmp , linearIndex , axisIndex , gridV ,
			   tmpMark , h0 , h1 , a_x_ij,  a_xx_ij);

	   // calculate diagonal according to 1.33
	   s_mu = 0.0; s_sig = 0.0; s_coor = 0.0; s_r = 0.0;
	   // convection
           #pragma simd reduction(+:s_mu)
	   for (int i = 0 ; i < nrFactors_ ; i++){
		   s_mu = s_mu + models()->getModel(i).convectionCoef(globCoord_tmp) *
				   V_ij_x_sp( i , h0[i] , h1[i] , gridV);
	   }
	   // diffusion
           #pragma simd reduction(+:s_sig)
	   for (int i = 0 ; i < nrFactors_ ; i++){
		   s_sig = s_sig + models()->getModel(i).diffusionCoef(globCoord_tmp) * models()->getModel(i).diffusionCoef(globCoord_tmp)
				   * V_ij_xx_sp( i , h0[i] , h1[i] , gridV);
	   }
	   // correlation
	   if ( this->hasCorrelations_){
		   for (int i = 0 ; i < nrFactors_ ; i++){
                           #pragma simd reduction(+:s_coor)
			   for (int j = i+1 ; j < nrFactors_ ; j++)
			   {
				   updateStage_xy( i , j , gridV_xy , gridV , linearIndex , axisIndex , tmpMark_xy);
				   s_coor = s_coor + models()->corr(i,j) *
						   models()->getModel(i).diffusionCoef(globCoord_tmp) * models()->getModel(j).diffusionCoef(globCoord_tmp)
						* ( V_ij_xy_sp( i , j , h0[i] , h1[i] , h0[j] , h1[j] , a_x_ij[i] , a_x_ij[j] , gridV_xy )
						  + a_x_ij[i] * V_ij_x_sp( j , h0[j] , h1[j] , gridV )	    );
			   }
		   }
	   }
	   // sum up all the terms
	   //grid()->res()[linearIndex] = grid()->u()[linearIndex] - (grid()->rhs()[linearIndex] + s_mu + s_sig + s_coor + s_r) / grid()->d()[linearIndex];
	   grid()->res()[linearIndex] = grid()->rhs()[linearIndex]
	               + ( 0.5*microTimestep_*(s_mu + 0.5*s_sig + s_coor) - grid()->u()[linearIndex] * grid()->d()[linearIndex] );
	   //FITOB_OUT_LEVEL3( verb , "calcResiduum() linearIndex:" << linearIndex << " , res():" << grid()->res()[linearIndex]);
	   // square the error
	   errorL2 += grid()->res()[linearIndex]*grid()->res()[linearIndex];
   }
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
}
#endif
    // - return the square root
    return sqrt(errorL2/(double)innerPoints);
}

void GaussSeidel_WB::calculateDiag( DVector& globCoord) {
	//int verb = 6;
	// ================ if diagonal is not calculated then calculate it =============================
		if (!isDiagonalCalculated_){
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp parallel
{
#endif
		   // variables to sum up the formula
		   double s_mu,s_sig,s_coor,s_r;
		   int linearIndex;
		   IVector axisIndex(gdim_,0);
		   DVector h0(nrFactors_);
		   DVector h1(nrFactors_);
		   DVector a_x_ij(nrFactors_);
		   DVector a_xx_ij(nrFactors_);
		   DVector gridV( 3 * nrFactors_  , 0.0);
		   IVector tmpMark( 3 * nrFactors_ , 0 );
		   DVector gridV_xy( 9 , 0.0 );
		   IVector tmpMark_xy( 9 , 0.0 );
		   DVector globCoord_tmp = globCoord;

#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp for schedule(static)
#endif
	       for (int pI = 0 ; pI < innerPoints ; pI++){

               // update the status
	    	   updateStage( pI , globCoord_tmp , linearIndex , axisIndex , gridV ,
	    			   tmpMark , h0 , h1 , a_x_ij,  a_xx_ij);

	    	   // calculate diagonal according to 1.33
	    	   s_mu = 0.0; s_sig = 0.0; s_coor = 0.0; s_r = 0.0;
	    	   // convection
                   #pragma simd reduction(+:s_mu)
	    	   for (int i = 0 ; i < nrFactors_ ; i++){
	    		   s_mu = s_mu - models()->getModel(i).convectionCoef(globCoord_tmp)*a_x_ij[i];
	    	   }
	    	   // diffusion
                   #pragma simd reduction(+:s_sig)
	    	   for (int i = 0 ; i < nrFactors_ ; i++){
	    		   s_sig = s_sig - a_xx_ij[i] * models()->getModel(i).diffusionCoef(globCoord_tmp) * models()->getModel(i).diffusionCoef(globCoord_tmp);
	    	   }
	    	   // correlation
	    	   if ( this->hasCorrelations_){
	    		   for (int i = 0 ; i < nrFactors_ ; i++){
	    			   for (int j = i+1 ; j < nrFactors_ ; j++)
	    			   {
	    				   // we do not need to update the buffer for the correlation
	    				   s_coor = s_coor - models()->corr(i,j) *
	    						   models()->getModel(i).diffusionCoef(globCoord_tmp) * models()->getModel(j).diffusionCoef(globCoord_tmp)
                                   * a_x_ij[i] * a_x_ij[j];
	    			   }
	    		   }
	    	   }
	    	   // interest rate
	    	   s_r = models()->r(globCoord_tmp);
	    	   // sum up all the terms
	    	   grid()->d()[linearIndex] = 1 + 0.5*microTimestep_*(s_mu + 0.5*s_sig + s_coor + s_r);
	  		   //FITOB_OUT_LEVEL3( verb , "calculateDiag linearIndex:" << linearIndex << " , d():" << grid()->d()[linearIndex]);
	       }
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
}
#endif
		}
	// ================ end diagonal calculation =============================
	isDiagonalCalculated_ = true;
}

void GaussSeidel_WB::oneSmoothGrid( DVector& globCoord){

#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp parallel
{
#endif
	// variables to sum up the formula
	double s_mu,s_sig,s_coor,s_r;
	int linearIndex;
	IVector axisIndex(gdim_,0);
	DVector h0(nrFactors_);
	DVector h1(nrFactors_);
	DVector a_x_ij(nrFactors_);
	DVector a_xx_ij(nrFactors_);
	DVector gridV( 3 * nrFactors_  , 0.0);
	IVector tmpMark( 3 * nrFactors_ , 0 );
	DVector gridV_xy( 9 , 0.0 );
	IVector tmpMark_xy( 9 , 0.0 );
	DVector globCoord_tmp = globCoord;
   //int verb = 6;

#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp for schedule(static)
#endif
   for (int pI = 0 ; pI < innerPoints ; pI++){

       // update the status
	   updateStage( pI , globCoord_tmp , linearIndex , axisIndex , gridV ,
			   tmpMark , h0 , h1 , a_x_ij,  a_xx_ij);

	   // calculate diagonal according to 1.33
	   s_mu = 0.0; s_sig = 0.0; s_coor = 0.0; s_r = 0.0;
	   // convection
           #pragma simd reduction(+:s_mu)
	   for (int i = 0 ; i < nrFactors_ ; i++){
		   s_mu = s_mu + models()->getModel(i).convectionCoef(globCoord_tmp) *
				   V_ij_x_sp( i , h0[i] , h1[i] , gridV );
	   }
	   // diffusion
           #pragma simd reduction(+:s_sig)
	   for (int i = 0 ; i < nrFactors_ ; i++){
		   //FITOB_OUT_LEVEL3(4, " factor [" << i << "] diffcoff =" << models()->getModel(i).diffusionCoef(globCoord) );
		   s_sig = s_sig + models()->getModel(i).diffusionCoef(globCoord_tmp) * models()->getModel(i).diffusionCoef(globCoord_tmp)
				   * V_ij_xx_sp( i , h0[i] , h1[i] , gridV );
	   }
	   // correlation
	   if ( this->hasCorrelations_){
		   for (int i = 0 ; i < nrFactors_ ; i++){
                           #pragma simd reduction(+:s_coor)
			   for (int j = i+1 ; j < nrFactors_ ; j++)
			   {
				   updateStage_xy( i , j , gridV_xy , gridV , linearIndex , axisIndex , tmpMark_xy);
				   s_coor = s_coor + models()->corr(i,j) *
						   models()->getModel(i).diffusionCoef(globCoord_tmp) * models()->getModel(j).diffusionCoef(globCoord_tmp)
						* ( V_ij_xy_sp( i , j , h0[i] , h1[i] , h0[j] , h1[j] , a_x_ij[i] , a_x_ij[j] , gridV_xy )
						  + a_x_ij[i] * V_ij_x_sp( j , h0[j] , h1[j] , gridV )	    );
			   }
		   }
	   }
	   // sum up all the terms
	   //FITOB_OUT_LEVEL3( verb , "oneSmoothGrid OLD linearIndex:" << linearIndex << " , u():" << grid()->u()[linearIndex]);
	   grid()->u()[linearIndex] = (grid()->rhs()[linearIndex] + 0.5*microTimestep_*(s_mu + 0.5*s_sig + s_coor)) / grid()->d()[linearIndex];
	   //FITOB_OUT_LEVEL3( verb , "oneSmoothGrid NEW linearIndex:" << linearIndex << " , u():" << grid()->u()[linearIndex]);
   }

#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp for schedule(static)
#endif
   // ---------- second iteration , backward iteration ------
   for (int pI = innerPoints-1 ; pI >= 0 ; pI--){
    // update the status
	   updateStage( pI , globCoord_tmp , linearIndex , axisIndex , gridV ,
			   tmpMark , h0 , h1 , a_x_ij,  a_xx_ij);

	   // calculate diagonal according to 1.33
	   s_mu = 0.0; s_sig = 0.0; s_coor = 0.0; s_r = 0.0;
	   // convection
           #pragma simd reduction(+:s_mu)
	   for (int i = 0 ; i < nrFactors_ ; i++){
		   s_mu = s_mu + models()->getModel(i).convectionCoef(globCoord_tmp) *
				   V_ij_x_sp( i , h0[i] , h1[i] , gridV );
	   }
	   // diffusion
           #pragma simd reduction(+:s_sig)
	   for (int i = 0 ; i < nrFactors_ ; i++){
		   s_sig = s_sig + models()->getModel(i).diffusionCoef(globCoord_tmp) * models()->getModel(i).diffusionCoef(globCoord_tmp)
				   * V_ij_xx_sp( i , h0[i] , h1[i] , gridV );
	   }
	   // correlation
	   if ( this->hasCorrelations_){
		   for (int i = 0 ; i < nrFactors_ ; i++){
                           #pragma simd reduction(+:s_coor)
			   for (int j = i+1 ; j < nrFactors_ ; j++)
			   {
				   updateStage_xy( i , j , gridV_xy , gridV , linearIndex , axisIndex , tmpMark_xy);
				   s_coor = s_coor + models()->corr(i,j) *
						   models()->getModel(i).diffusionCoef(globCoord_tmp) * models()->getModel(j).diffusionCoef(globCoord_tmp)
						* ( V_ij_xy_sp( i , j , h0[i] , h1[i] , h0[j] , h1[j] , a_x_ij[i] , a_x_ij[j] , gridV_xy)
						  + a_x_ij[i] * V_ij_x_sp( j , h0[j] , h1[j] , gridV )	    );
			   }
		   }
	   }
	   // sum up all the terms
	   //FITOB_OUT_LEVEL3( verb , "oneSmoothGrid OLD linearIndex:" << linearIndex << " , u():" << grid()->u()[linearIndex]);
	   grid()->u()[linearIndex] = (grid()->rhs()[linearIndex] + 0.5*microTimestep_*(s_mu + 0.5*s_sig + s_coor)) / grid()->d()[linearIndex];
	   //FITOB_OUT_LEVEL3( verb , "oneSmoothGrid NEW linearIndex:" << linearIndex << " , u():" << grid()->u()[linearIndex]);
   }
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
}
#endif
}

void GaussSeidel_WB::explicitStep( DVector& globCoord , DVector& unknowns , double underrelaxFactor ){

#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp parallel
{
#endif
	// test if the sizes match
	// variables to sum up the formula
	double s_mu,s_sig,s_coor,s_r;
	int linearIndex;
	IVector axisIndex(gdim_,0);
	DVector h0(nrFactors_);
	DVector h1(nrFactors_);
	DVector a_x_ij(nrFactors_);
	DVector a_xx_ij(nrFactors_);
	DVector gridV( 3 * nrFactors_  , 0.0);
	IVector tmpMark( 3 * nrFactors_ , 0 );
	DVector gridV_xy( 9 , 0.0 );
	IVector tmpMark_xy( 9 , 0.0 );
	DVector globCoord_tmp = globCoord;
	//int verb = 6;

	double tmp_val ;
	//underrelaxFactor = 1e-2;

	FITOB_ERROR_TEST(unknowns.size() == grid()->u().size() , " GaussSeidel_WB::explicitStep , unknowns.size():" << unknowns.size()
			        << " , grid()->u().size():" << grid()->u().size() );

#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp for schedule(static)
#endif
	// ------ overwrite the unknown vector first -----
    for (unsigned int pI = 0 ; pI < grid()->u().size() ; pI++){
    	unknowns[pI] = grid()->u()[pI];
    }

#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp for schedule(static)
#endif
    // -------- make the first step , till the mid point, from t till t+0.5*dt -----
    for (int pI = 0 ; pI < innerPoints ; pI++){
        // update the status
 	   updateStage( pI , globCoord_tmp , linearIndex , axisIndex , gridV ,
			   tmpMark , h0 , h1 , a_x_ij,  a_xx_ij);
 	   // calculate diagonal according to 1.33
 	   s_mu = 0.0; s_sig = 0.0; s_coor = 0.0; s_r = 0.0;
 	   // convection
           #pragma simd reduction(+:s_mu)
 	   for (int i = 0 ; i < nrFactors_ ; i++){
 		   s_mu = s_mu + models()->getModel(i).convectionCoef(globCoord_tmp) *
 				   V_ij_x( i , h0[i] , h1[i] , a_x_ij[i] , gridV);
 	   }
 	   // diffusion
           #pragma simd reduction(+:s_sig)
 	   for (int i = 0 ; i < nrFactors_ ; i++){
 		   s_sig = s_sig + models()->getModel(i).diffusionCoef(globCoord_tmp) * models()->getModel(i).diffusionCoef(globCoord_tmp)
 				   * V_ij_xx( i , h0[i] , h1[i] , a_xx_ij[i] , gridV );
 	   }
 	   // correlation
 	   if ( this->hasCorrelations_){
 		   for (int i = 0 ; i < nrFactors_ ; i++){
                           #pragma simd reduction(+:s_coor)
 			   for (int j = i+1 ; j < nrFactors_ ; j++)
 			   {
 				   updateStage_xy( i , j , gridV_xy , gridV , linearIndex , axisIndex , tmpMark_xy);
 				   s_coor = s_coor + models()->corr(i,j) *
 						   models()->getModel(i).diffusionCoef(globCoord_tmp) * models()->getModel(j).diffusionCoef(globCoord_tmp)
 						 * V_ij_xy( i , j , h0[i] , h1[i] , h0[j] , h1[j] , a_x_ij[i] , a_x_ij[j] , gridV_xy );
 			   }
 		   }
 	   }
 	   // interest rate
 	   s_r = - models()->r(globCoord_tmp) * grid()->u()[linearIndex];
 	   // sum up all the terms
 	   //FITOB_OUT_LEVEL3( verb , " Explicit step linearIndex:" << linearIndex << " , s_mu:" << s_mu << " , s_sig:" << s_sig <<
 		//	                 " , s_coor:" << s_coor << " , s_r:" << s_r );
 	   unknowns[linearIndex] = grid()->u()[linearIndex] + underrelaxFactor * 0.5 * microTimestep_ * (s_mu + 0.5*s_sig + s_coor + s_r);
 	   //FITOB_OUT_LEVEL3( verb , " old:" << grid()->u()[linearIndex] << " , new:" << unknowns[linearIndex]);
    }

#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp for schedule(static)
#endif
    // ------ swap "unknowns" and "grid()->u()" -----
    for (unsigned int pI = 0 ; pI < grid()->u().size() ; pI++){
    	tmp_val = unknowns[pI];
    	unknowns[pI] = grid()->u()[pI];
    	grid()->u()[pI] = tmp_val;
    }

#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp for schedule(static)
#endif
    // ----- make the second step, complete step from t till t+dt --------
    for (int pI = 0 ; pI < innerPoints ; pI++){
        // update the status
 	   updateStage( pI , globCoord_tmp , linearIndex , axisIndex , gridV ,
			   tmpMark , h0 , h1 , a_x_ij,  a_xx_ij);
 	   // calculate diagonal according to 1.33
 	   s_mu = 0.0; s_sig = 0.0; s_coor = 0.0; s_r = 0.0;
 	   // convection
           #pragma simd reduction(+:s_mu)
 	   for (int i = 0 ; i < nrFactors_ ; i++){
 		   s_mu = s_mu + models()->getModel(i).convectionCoef(globCoord_tmp) *
 				   V_ij_x( i , h0[i] , h1[i] , a_x_ij[i] , gridV );
 	   }
 	   // diffusion
           #pragma simd reduction(+:s_sig)
 	   for (int i = 0 ; i < nrFactors_ ; i++){
 		   s_sig = s_sig + models()->getModel(i).diffusionCoef(globCoord_tmp) * models()->getModel(i).diffusionCoef(globCoord_tmp)
 				   * V_ij_xx( i , h0[i] , h1[i] , a_xx_ij[i] , gridV );
 	   }
 	   // correlation
 	   if ( this->hasCorrelations_){
 		   for (int i = 0 ; i < nrFactors_ ; i++){
                           #pragma simd reduction(+:s_coor)
 			   for (int j = i+1 ; j < nrFactors_ ; j++)
 			   {
 				   updateStage_xy( i , j , gridV_xy , gridV , linearIndex , axisIndex , tmpMark_xy);
 				   s_coor = s_coor + models()->corr(i,j) *
 						   models()->getModel(i).diffusionCoef(globCoord_tmp) * models()->getModel(j).diffusionCoef(globCoord_tmp)
 						 * V_ij_xy( i , j , h0[i] , h1[i] , h0[j] , h1[j] , a_x_ij[i] , a_x_ij[j] , gridV_xy );
 			   }
 		   }
 	   }
 	   // interest rate
 	   s_r = - models()->r(globCoord_tmp) * grid()->u()[linearIndex];
 	   // sum up all the terms
 	   // here we make a complete step
 	   unknowns[linearIndex] = grid()->u()[linearIndex] + underrelaxFactor*microTimestep_*(s_mu + 0.5*s_sig + s_coor + s_r);
    }

#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp for schedule(static)
#endif
    // copy the solution back
    for (unsigned int pI = 0 ; pI < grid()->u().size() ; pI++){
    	grid()->u()[pI] = unknowns[pI];
    }
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
}
#endif
}

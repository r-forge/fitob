/*
 * FitobGaussSeidel.cpp
 *
 *  Created on: Jul 15, 2010
 *      Author: benk
 */

#include "FitobGaussSeidel.hpp"

using namespace fitob;
using namespace std;

GaussSeidel::GaussSeidel( const XMLConfiguration* config , MultigridFGBase *grid ,
		const ModelCollection* models , const FitobCalculator* fitobcalculator)
:SmootherBase(config,grid, models, fitobcalculator){

	setVerb(4);

	isDiagonalCalculated_ = false;
	calculateRHS_ = false;

    factorIndex_to_local.resize(nrFactors_);
    innerOffset.resize(gdim_);

	innerPoints = 1; // number of points where we'll iterate
	int linearIndex = 0;
	// the local index of each factor

    // copy the constant global coordinates
    // todo: update time in the global coordinates (or this should be done somewhere above)
    const ModelCollection* models_tmp = models;

    linearIndex = 1;
	for (int i = 0 ; i < gdim_ ; i++){
		// calculate the inner offset
		innerOffset[i] = linearIndex;
		linearIndex = linearIndex*(grid->nrPoints(i)-2);
	    // the inner points are
		innerPoints = innerPoints*(grid->nrPoints(i)-2);
		FITOB_OUT_LEVEL3(verb()," GaussSeidl::GaussSeidl: i:"<<i<<" , linearIndex:"<<linearIndex<<" , innerOffset[i]:"
				             << innerOffset[i] << " , innerPoints:"<<innerPoints);
	}
	for (int i = 0 ; i < nrFactors_ ; i++){
	    // get the local index of each factor
		const FactorModelBase& fm = models_tmp->getModel(i);
		int glIndex = fm.getGlobalIndex();
		factorIndex_to_local[i] = grid->domain()->globalToLocalIndex( glIndex );
	}
}


void GaussSeidel::smoothGrid( const DVector& globalCoords) {
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


void GaussSeidel::calculateDiag( DVector& globCoord){

	// ================ if diagonal is not calculated then calculate it =============================
		if (!isDiagonalCalculated_){
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp parallel
{
#endif
		   // variables to sum up the formula
		   double s_mu,s_sig,s_coor,s_r;
		   int linearIndex;
		   IVector axisIndex(gdim_);
		   DVector h0(nrFactors_);
		   DVector h1(nrFactors_);
		   DVector a_x_ij(nrFactors_);
		   DVector a_xx_ij(nrFactors_);
		   DVector globCoord_tmp = globCoord;

#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp for schedule(static)
#endif
	       for (int pI = 0 ; pI < innerPoints ; pI++){

               // update the status
	    	   updateStage( pI , globCoord_tmp , linearIndex, axisIndex, h0 , h1 , a_x_ij , a_xx_ij);

	    	   // calculate diagonal according to 1.33
	    	   s_mu = 0.0; s_sig = 0.0; s_coor = 0.0; s_r = 0.0;
	    	   // convection
	    	   for (int i = 0 ; i < nrFactors_ ; i++){
	    		   s_mu = s_mu - models()->getModel(i).convectionCoef(globCoord_tmp)*a_x_ij[i];
	    	   }
	    	   // diffusion
	    	   for (int i = 0 ; i < nrFactors_ ; i++){
	    		   s_sig = s_sig - a_xx_ij[i] * models()->getModel(i).diffusionCoef(globCoord_tmp) * models()->getModel(i).diffusionCoef(globCoord_tmp);
	    	   }
	    	   // correlation
	    	   if ( this->hasCorrelations_){
	    		   for (int i = 0 ; i < nrFactors_ ; i++){
	    			   for (int j = i+1 ; j < nrFactors_ ; j++)
	    			   {
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
	       }
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
}
#endif
		}
	// ================ end diagonal calculation =============================
	isDiagonalCalculated_ = true;
}

void GaussSeidel::calculateRHS( DVector& globCoord){
  // ====================== calculate right hand side =====================
	// in case of the multigrid only at the highest level we should do this in correction scheme
	// in FAS case on each level we should do this and similar things on each level
	if (!calculateRHS_){
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp parallel
{
#endif
		   // variables to sum up the formula
		   double s_mu,s_sig,s_coor,s_r;
		   int linearIndex;
		   IVector axisIndex(gdim_);
		   DVector h0(nrFactors_);
		   DVector h1(nrFactors_);
		   DVector a_x_ij(nrFactors_);
		   DVector a_xx_ij(nrFactors_);
		   DVector globCoord_tmp = globCoord;

#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp for schedule(static)
#endif
	       for (int pI = 0 ; pI < innerPoints ; pI++){

               // update the status
	    	   updateStage( pI , globCoord_tmp , linearIndex, axisIndex, h0 , h1 , a_x_ij , a_xx_ij );

	    	   // calculate diagonal according to 1.33
	    	   s_mu = 0.0; s_sig = 0.0; s_coor = 0.0; s_r = 0.0;
	    	   // convection
	    	   for (int i = 0 ; i < nrFactors_ ; i++){
	    		   s_mu = s_mu + models()->getModel(i).convectionCoef(globCoord_tmp) *
	    				   V_ij_x( linearIndex , grid()->offset()[i] , h0[i] , h1[i] , a_x_ij[i] );
	    	   }
	    	   // diffusion
	    	   for (int i = 0 ; i < nrFactors_ ; i++){
	    		   s_sig = s_sig + models()->getModel(i).diffusionCoef(globCoord_tmp) * models()->getModel(i).diffusionCoef(globCoord_tmp)
	    				   * V_ij_xx( linearIndex , grid()->offset()[i] , h0[i] , h1[i] , a_xx_ij[i] );
	    	   }
	    	   // correlation
	    	   if ( this->hasCorrelations_){
	    		   for (int i = 0 ; i < nrFactors_ ; i++)
	    		   {
	    			   for (int j = i+1 ; j < nrFactors_ ; j++)
	    			   {
	    				   s_coor = s_coor + models()->corr(i,j) *
	    						   models()->getModel(i).diffusionCoef(globCoord_tmp) * models()->getModel(j).diffusionCoef(globCoord_tmp)
	    						 * V_ij_xy( linearIndex , grid()->offset()[i] , grid()->offset()[j] ,
	    								    h0[i] , h1[i] , h0[j] , h1[j] , a_x_ij[i] , a_x_ij[j]     );
	    			   }
	    		   }
	    	   }
	    	   // interest rate
	    	   s_r = - models()->r(globCoord_tmp) * grid()->u()[linearIndex];

	    	   // sum up all the terms
	    	   grid()->rhs()[linearIndex] = grid()->u()[linearIndex] + 0.5*microTimestep_*(s_mu + 0.5*s_sig + s_coor + s_r);
	       }
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
}
#endif
	}
	calculateRHS_ = true;
  // ====================== end of calculating right hand side =====================
}


void GaussSeidel::oneSmoothGrid( DVector& globCoord){

#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp parallel
{
#endif
	// variables to sum up the formula
	double s_mu,s_sig,s_coor,s_r;
	int linearIndex;
	IVector axisIndex(gdim_);
	DVector h0(nrFactors_);
	DVector h1(nrFactors_);
	DVector a_x_ij(nrFactors_);
	DVector a_xx_ij(nrFactors_);
	DVector globCoord_tmp = globCoord;

#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp for schedule(static)
#endif
    for (int pI = 0 ; pI < innerPoints ; pI++){

       // update the status
 	   updateStage( pI , globCoord_tmp , linearIndex, axisIndex, h0 , h1 , a_x_ij , a_xx_ij );

 	   // calculate diagonal according to 1.33
 	   s_mu = 0.0; s_sig = 0.0; s_coor = 0.0; s_r = 0.0;
 	   // convection
 	   for (int i = 0 ; i < nrFactors_ ; i++){
 		   s_mu = s_mu + models()->getModel(i).convectionCoef(globCoord_tmp) *
 				   V_ij_x_sp( linearIndex , grid()->offset()[i] , h0[i] , h1[i] );
 	   }
 	   // diffusion
 	   for (int i = 0 ; i < nrFactors_ ; i++){
  		   //FITOB_OUT_LEVEL3(4, " factor [" << i << "] diffcoff =" << models()->getModel(i).diffusionCoef(globCoord) );
 		   s_sig = s_sig + models()->getModel(i).diffusionCoef(globCoord_tmp) * models()->getModel(i).diffusionCoef(globCoord_tmp)
 				   * V_ij_xx_sp( linearIndex , grid()->offset()[i] , h0[i] , h1[i] );
 	   }
 	   // correlation
 	   if ( this->hasCorrelations_){
 		   for (int i = 0 ; i < nrFactors_ ; i++){
 			   for (int j = i+1 ; j < nrFactors_ ; j++)
 			   {
 				   s_coor = s_coor + models()->corr(i,j) *
 						   models()->getModel(i).diffusionCoef(globCoord_tmp) * models()->getModel(j).diffusionCoef(globCoord_tmp)
 						* ( V_ij_xy_sp( linearIndex , grid()->offset()[i] , grid()->offset()[j] ,
 								        h0[i] , h1[i] , h0[j] , h1[j] , a_x_ij[i] , a_x_ij[j] )
 						  + a_x_ij[i] * V_ij_x_sp( linearIndex , grid()->offset()[j] , h0[j] , h1[j] )	    );
 			   }
 		   }
 	   }
 	   // sum up all the terms
 	   grid()->u()[linearIndex] = (grid()->rhs()[linearIndex] + 0.5*microTimestep_*(s_mu + 0.5*s_sig + s_coor)) / grid()->d()[linearIndex];
    }

#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp for schedule(static)
#endif
    // ---------- second iteration , backward iteration ------
    for (int pI = innerPoints-1 ; pI >= 0 ; pI--){
       // update the status
 	   updateStage( pI , globCoord_tmp , linearIndex, axisIndex, h0 , h1 , a_x_ij , a_xx_ij );

 	   // calculate diagonal according to 1.33
 	   s_mu = 0.0; s_sig = 0.0; s_coor = 0.0; s_r = 0.0;
 	   // convection
 	   for (int i = 0 ; i < nrFactors_ ; i++){
 		   s_mu = s_mu + models()->getModel(i).convectionCoef(globCoord_tmp) *
 				   V_ij_x_sp( linearIndex , grid()->offset()[i] , h0[i] , h1[i] );
 	   }
 	   // diffusion
 	   for (int i = 0 ; i < nrFactors_ ; i++){
 		   s_sig = s_sig + models()->getModel(i).diffusionCoef(globCoord_tmp) * models()->getModel(i).diffusionCoef(globCoord_tmp)
 				   * V_ij_xx_sp( linearIndex , grid()->offset()[i] , h0[i] , h1[i] );
 	   }
 	   // correlation
 	   if ( this->hasCorrelations_){
 		   for (int i = 0 ; i < nrFactors_ ; i++){
 			   for (int j = i+1 ; j < nrFactors_ ; j++)
 			   {
 				   s_coor = s_coor + models()->corr(i,j) *
 						   models()->getModel(i).diffusionCoef(globCoord_tmp) * models()->getModel(j).diffusionCoef(globCoord_tmp)
 						* ( V_ij_xy_sp( linearIndex , grid()->offset()[i] , grid()->offset()[j] ,
 								        h0[i] , h1[i] , h0[j] , h1[j] , a_x_ij[i] , a_x_ij[j] )
 						  + a_x_ij[i] * V_ij_x_sp( linearIndex , grid()->offset()[j] , h0[j] , h1[j] )	    );
 			   }
 		   }
 	   }
 	   // sum up all the terms
 	   grid()->u()[linearIndex] = (grid()->rhs()[linearIndex] + 0.5*microTimestep_*(s_mu + 0.5*s_sig + s_coor)) / grid()->d()[linearIndex];
    }
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
}
#endif
}

double GaussSeidel::calcResiduum( const DVector& globalCoords ){

	double errorL2 = 0.0;
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp parallel reduction(+ : errorL2)
{
#endif
	// -- setup phase for the iteration
    DVector globCoord = globalCoords;
	// variables to sum up the formula
	double s_mu,s_sig,s_coor,s_r;
	int linearIndex;
	IVector axisIndex(gdim_);
	DVector h0(nrFactors_);
	DVector h1(nrFactors_);
	DVector a_x_ij(nrFactors_);
	DVector a_xx_ij(nrFactors_);
	DVector globCoord_tmp = globCoord;

#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp for schedule(static)
#endif
    for (int pI = 0 ; pI < innerPoints ; pI++){

       // update the status
	   updateStage( pI , globCoord_tmp , linearIndex, axisIndex, h0 , h1 , a_x_ij , a_xx_ij );

	   // calculate diagonal according to 1.33
	   s_mu = 0.0; s_sig = 0.0; s_coor = 0.0; s_r = 0.0;
	   // convection
	   for (int i = 0 ; i < nrFactors_ ; i++){
		   s_mu = s_mu + models()->getModel(i).convectionCoef(globCoord_tmp) *
				   V_ij_x_sp( linearIndex , grid()->offset()[i] , h0[i] , h1[i] );
	   }
	   // diffusion
	   for (int i = 0 ; i < nrFactors_ ; i++){
		   s_sig = s_sig + models()->getModel(i).diffusionCoef(globCoord_tmp) * models()->getModel(i).diffusionCoef(globCoord_tmp)
				   * V_ij_xx_sp( linearIndex , grid()->offset()[i] , h0[i] , h1[i] );
	   }
	   // correlation
	   if ( this->hasCorrelations_){
		   for (int i = 0 ; i < nrFactors_ ; i++){
			   for (int j = i+1 ; j < nrFactors_ ; j++)
			   {
				   s_coor = s_coor + models()->corr(i,j) *
						   models()->getModel(i).diffusionCoef(globCoord_tmp) * models()->getModel(j).diffusionCoef(globCoord_tmp)
						* ( V_ij_xy_sp( linearIndex , grid()->offset()[i] , grid()->offset()[j] ,
								        h0[i] , h1[i] , h0[j] , h1[j] , a_x_ij[i] , a_x_ij[j] )
						  + a_x_ij[i] * V_ij_x_sp( linearIndex , grid()->offset()[j] , h0[j] , h1[j] )	    );
			   }
		   }
	   }
	   // sum up all the terms
	   //grid()->res()[linearIndex] = grid()->u()[linearIndex] - (grid()->rhs()[linearIndex] + s_mu + s_sig + s_coor + s_r) / grid()->d()[linearIndex];
	   grid()->res()[linearIndex] = grid()->rhs()[linearIndex] + ( 0.5*microTimestep_*(s_mu + 0.5*s_sig + s_coor) - grid()->u()[linearIndex] * grid()->d()[linearIndex] );
	   // square the error
	   errorL2 += grid()->res()[linearIndex]*grid()->res()[linearIndex];
   }
#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
}
#endif
    // - return the square root
    return sqrt(errorL2/(double)innerPoints);
}


// ===================== EXPLICIT PREDICTOR TIME STEP =======================
// The code below is the midpoint rule code

void GaussSeidel::explicitStep( DVector& globCoord , DVector& unknowns , double underrelaxFactor ) {

#if ( defined(FITOB_OPENMP) && defined(FITOB_OPENMP_SOLVER) )
#pragma omp parallel
{
#endif
	// test if the sizes match
	// variables to sum up the formula
	double s_mu,s_sig,s_coor,s_r;
	int linearIndex;
	IVector axisIndex(gdim_);
	DVector h0(nrFactors_);
	DVector h1(nrFactors_);
	DVector a_x_ij(nrFactors_);
	DVector a_xx_ij(nrFactors_);
	DVector globCoord_tmp = globCoord;
	//int verb = 0;

	double tmp_val ;
	//underrelaxFactor = 1e-2;

	FITOB_ERROR_TEST(unknowns.size() == grid()->u().size() , " GaussSeidel::explicitStep , unknowns.size():" << unknowns.size()
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
 	   updateStage( pI , globCoord_tmp , linearIndex, axisIndex, h0 , h1 , a_x_ij , a_xx_ij );
 	   // calculate diagonal according to 1.33
 	   s_mu = 0.0; s_sig = 0.0; s_coor = 0.0; s_r = 0.0;
 	   // convection
 	   for (int i = 0 ; i < nrFactors_ ; i++){
 		   s_mu = s_mu + models()->getModel(i).convectionCoef(globCoord_tmp) *
 				   V_ij_x( linearIndex , grid()->offset()[i] , h0[i] , h1[i] , a_x_ij[i] );
 		  /*FITOB_OUT_LEVEL3( verb , " convection: " << models()->getModel(i).convectionCoef(globCoord) );
 		  FITOB_OUT_LEVEL3( verb , " V_ij_x:" << V_ij_x( linearIndex , grid()->offset()[i] , h0[i] , h1[i] , a_x_ij[i] ) );
 		  FITOB_OUT_LEVEL3( verb , " V_ij_x_sp:" << V_ij_x_sp( linearIndex , grid()->offset()[i] , h0[i] , h1[i] ) );
 		  FITOB_OUT_LEVEL3( verb , "(hh0/hh1):" << h0[i]/h1[i] << " , (hh1/hh0):" << h1[i]/h0[i] << " , (hh0+hh1):" << h0[i]+h1[i]);
 		  FITOB_OUT_LEVEL3( verb , "grid()->u()[linearIndex-grid()->offset()[i]]:" << grid()->u()[linearIndex-grid()->offset()[i]] );
 		  FITOB_OUT_LEVEL3( verb , "grid()->u()[linearIndex]:" << grid()->u()[linearIndex] );
 		  FITOB_OUT_LEVEL3( verb , "grid()->u()[linearIndex+grid()->offset()[i]]:" << grid()->u()[linearIndex+grid()->offset()[i]] ); */
 	   }
 	   // diffusion
 	   for (int i = 0 ; i < nrFactors_ ; i++){
 		   s_sig = s_sig + models()->getModel(i).diffusionCoef(globCoord_tmp) * models()->getModel(i).diffusionCoef(globCoord_tmp)
 				   * V_ij_xx( linearIndex , grid()->offset()[i] , h0[i] , h1[i] , a_xx_ij[i] );
 		  /*FITOB_OUT_LEVEL3( verb , " diffusion: " << models()->getModel(i).diffusionCoef(globCoord) );
 		  FITOB_OUT_LEVEL3( verb , " V_ij_xx:" << V_ij_xx( linearIndex , grid()->offset()[i] , h0[i] , h1[i] , a_xx_ij[i] )  << " , a_xx_ij[i]:" << a_xx_ij[i]);
 		  FITOB_OUT_LEVEL3( verb , " V_ij_xx_sp:" << V_ij_xx_sp(linearIndex , grid()->offset()[i] , h0[i] , h1[i]) );
 		  FITOB_OUT_LEVEL3( verb , "0.5*hh0*hh1*(hh1+hh0)" << 0.5*h0[i]*h1[i]*(h1[i]+h0[i]) );*/
 	   }
 	   // correlation
 	   if ( this->hasCorrelations_){
 		   for (int i = 0 ; i < nrFactors_ ; i++){
 			   for (int j = i+1 ; j < nrFactors_ ; j++)
 			   {
 				   s_coor = s_coor + models()->corr(i,j) *
 						   models()->getModel(i).diffusionCoef(globCoord_tmp) * models()->getModel(j).diffusionCoef(globCoord_tmp)
 						 * V_ij_xy( linearIndex , grid()->offset()[i] , grid()->offset()[j] ,
 								    h0[i] , h1[i] , h0[j] , h1[j] , a_x_ij[i] , a_x_ij[j]     );
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
 	   updateStage( pI , globCoord_tmp , linearIndex, axisIndex, h0 , h1 , a_x_ij , a_xx_ij );
 	   // calculate diagonal according to 1.33
 	   s_mu = 0.0; s_sig = 0.0; s_coor = 0.0; s_r = 0.0;
 	   // convection
 	   for (int i = 0 ; i < nrFactors_ ; i++){
 		   s_mu = s_mu + models()->getModel(i).convectionCoef(globCoord_tmp) *
 				   V_ij_x( linearIndex , grid()->offset()[i] , h0[i] , h1[i] , a_x_ij[i] );
 	   }
 	   // diffusion
 	   for (int i = 0 ; i < nrFactors_ ; i++){
 		   s_sig = s_sig + models()->getModel(i).diffusionCoef(globCoord_tmp) * models()->getModel(i).diffusionCoef(globCoord_tmp)
 				   * V_ij_xx( linearIndex , grid()->offset()[i] , h0[i] , h1[i] , a_xx_ij[i] );
 	   }
 	   // correlation
 	   if ( this->hasCorrelations_){
 		   for (int i = 0 ; i < nrFactors_ ; i++){
 			   for (int j = i+1 ; j < nrFactors_ ; j++)
 			   {
 				   s_coor = s_coor + models()->corr(i,j) *
 						   models()->getModel(i).diffusionCoef(globCoord_tmp) * models()->getModel(j).diffusionCoef(globCoord_tmp)
 						 * V_ij_xy( linearIndex , grid()->offset()[i] , grid()->offset()[j] ,
 								    h0[i] , h1[i] , h0[j] , h1[j] , a_x_ij[i] , a_x_ij[j]     );
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

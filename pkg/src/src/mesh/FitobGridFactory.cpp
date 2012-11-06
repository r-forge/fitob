/*
 * FitobGridFactory.cpp
 *
 *  Created on: Jul 5, 2010
 *      Author: benk
 */

#include "src/mesh/FitobGridFactory.hpp"
#include "src/mesh/FitobConstantMesh.hpp"
#include "src/mesh/FitobFullGrid.hpp"
#include "src/mesh/FitobFullGrid_WB.hpp"
#include "src/mesh/FitobSparseGrid.hpp"
#include "src/mesh/FitobSGppSparseGrid.hpp"
#include "src/mesh/FitobSGppSparseGridAdaptive.hpp"
#include "src/mesh/FitobCombiGrid.hpp"
#include "src/mesh/FitobExtrapGrid.hpp"

#include "src/scripteval/FitobCalculator.hpp"

using namespace fitob;
using namespace std;

GridFactory::GridFactory(boost::shared_ptr<XMLConfiguration>& XMLConfig) :
XMLConfig_(XMLConfig.get()) , dimensionAdaptiveTruncationLevels_(0) {

	setVerb(0);

	string gridName = XMLConfig->getStringConfiguration("thetaconfigurations.gridproperties.GRID_TYPE.<xmlattr>.value");

	typeOfGrid_ = GRID_WITH_BOUNDARY;
	gridtype_ = -1;
	special_diagonal_cut_off_level_ = 0.0;
	use_OptiCom_ = false;
	use_Scaling_ = true;

// ------------------ SGpp grids --------------------------
	// build in the grids the grids without boundary points
	if (gridName == "sparsegrid"){
		FITOB_OUT_LEVEL3(verb(),"GridFactory set to sparsegrid : 1" );
		gridtype_ = 1;
	}
	if (gridName == "sparsegrid_WB"){
		FITOB_OUT_LEVEL3(verb(),"GridFactory WB set to sparsegrid : 1" );
		gridtype_ = 1; typeOfGrid_ = GRID_WITHOUT_BOUNDARY;
	}
	if (gridName == "SQRTgrid"){
		FITOB_OUT_LEVEL3(verb(),"GridFactory set to sparsegrid (SQRT): 2" );
		gridtype_ = 2;
	}
	if (gridName == "SQRTgrid_WB"){
		FITOB_OUT_LEVEL3(verb(),"GridFactory WB set to sparsegrid (SQRT): 2" );
		gridtype_ = 2; typeOfGrid_ = GRID_WITHOUT_BOUNDARY;
	}
	if (gridName == "truncated_sparsegrid"){
		FITOB_OUT_LEVEL3(verb(),"GridFactory set to sparsegrid (truncated_sparsegrid) : 3" );
		special_diagonal_cut_off_level_ =
			XMLConfig->getDoubleConfiguration("thetaconfigurations.gridproperties.combigrid.COMBI_DIAGONAL_CUT_LEVEL.<xmlattr>.value");
		gridtype_ = 3;
	}
	if (gridName == "truncated_sparsegrid_WB"){
		FITOB_OUT_LEVEL3(verb(),"GridFactory WB set to sparsegrid (truncated_sparsegrid) : 3" );
		special_diagonal_cut_off_level_ =
			XMLConfig->getDoubleConfiguration("thetaconfigurations.gridproperties.combigrid.COMBI_DIAGONAL_CUT_LEVEL.<xmlattr>.value");
		gridtype_ = 3; typeOfGrid_ = GRID_WITHOUT_BOUNDARY;
	}
	if (gridName == "sparsegrid_stretched"){
		FITOB_OUT_LEVEL3(verb(),"GridFactory WB set to sparsegrid (sparsegrid_stretched) : 4" );
		gridtype_ = 4;
	}

// ----------------------- COMBINATION GRID --------------------
	if (gridName == "combigrid"){
		FITOB_OUT_LEVEL3(verb(),"GridFactory set to combigrid : 1" );
		gridtype_ = 6;
	}
	if (gridName == "combigrid_WB"){
		FITOB_OUT_LEVEL3(verb(),"GridFactory WB set to combigrid : 1" );
		gridtype_ = 6; typeOfGrid_ = GRID_WITHOUT_BOUNDARY;
	}
	if (gridName == "combiSQRTgrid"){
		FITOB_OUT_LEVEL3(verb(),"GridFactory set to combigrid (SQRT): 2" );
		gridtype_ = 7;
	}
	if (gridName == "combiSQRTgrid_WB"){
		FITOB_OUT_LEVEL3(verb(),"GridFactory WB set to combigrid (SQRT): 2" );
		gridtype_ = 7; typeOfGrid_ = GRID_WITHOUT_BOUNDARY;
	}
	if (gridName == "truncated_combigrid"){
		dimensionAdaptiveTruncationLevels_.resize(0);
		FITOB_OUT_LEVEL3(verb(),"GridFactory set to combigrid (truncated_sparsegrid) : 3" );
		special_diagonal_cut_off_level_ =
		XMLConfig->getDoubleConfiguration("thetaconfigurations.gridproperties.combigrid.COMBI_DIAGONAL_CUT_LEVEL.<xmlattr>.value");
		XMLConfig->getDoubleVectorConfiguration(
				"thetaconfigurations.gridproperties.combigrid.COMBI_DIAGONAL_CUT_LEVEL.<xmlattr>.vector",
                ',' , dimensionAdaptiveTruncationLevels_);
		gridtype_ = 8;
	}
	if (gridName == "truncated_combigrid_WB"){
		dimensionAdaptiveTruncationLevels_.resize(0);
		FITOB_OUT_LEVEL3(verb(),"GridFactory WB set to combigrid (truncated_sparsegrid) : 3" );
		special_diagonal_cut_off_level_ =
		XMLConfig->getDoubleConfiguration("thetaconfigurations.gridproperties.combigrid.COMBI_DIAGONAL_CUT_LEVEL.<xmlattr>.value");
		XMLConfig->getDoubleVectorConfiguration(
				"thetaconfigurations.gridproperties.combigrid.COMBI_DIAGONAL_CUT_LEVEL.<xmlattr>.vector",
                ',' , dimensionAdaptiveTruncationLevels_);
		gridtype_ = 8; typeOfGrid_ = GRID_WITHOUT_BOUNDARY;
	}


// --------------------- EXTRAPOLATION GRID -----------------
	if (gridName == "extrap_combigrid"){
		FITOB_OUT_LEVEL3(verb(),"GridFactory set to extrap_combigrid : 1" );
		gridtype_ = 16;
	}
	if (gridName == "extrap_combigrid_WB"){
		FITOB_OUT_LEVEL3(verb(),"GridFactory WB set to extrap_combigrid : 1" );
		gridtype_ = 16; typeOfGrid_ = GRID_WITHOUT_BOUNDARY;
	}
	if (gridName == "extrap_combiSQRTgrid"){
		FITOB_OUT_LEVEL3(verb(),"GridFactory set to extrap_combigrid (SQRT): 2" );
		gridtype_ = 17;
	}
	if (gridName == "extrap_combiSQRTgrid_WB"){
		FITOB_OUT_LEVEL3(verb(),"GridFactory WB set to extrap_combigrid (SQRT): 2" );
		gridtype_ = 17; typeOfGrid_ = GRID_WITHOUT_BOUNDARY;
	}
	if (gridName == "extrap_truncated_combigrid"){
		dimensionAdaptiveTruncationLevels_.resize(0);
		FITOB_OUT_LEVEL3(verb(),"GridFactory set to extrap_truncated_combigrid : 3" );
		special_diagonal_cut_off_level_ =
		XMLConfig->getDoubleConfiguration("thetaconfigurations.gridproperties.combigrid.COMBI_DIAGONAL_CUT_LEVEL.<xmlattr>.value");
		XMLConfig->getDoubleVectorConfiguration(
				"thetaconfigurations.gridproperties.combigrid.COMBI_DIAGONAL_CUT_LEVEL.<xmlattr>.vector",
                ',' , dimensionAdaptiveTruncationLevels_);
		gridtype_ = 18;
	}
	if (gridName == "extrap_truncated_combigrid_WB"){
		dimensionAdaptiveTruncationLevels_.resize(0);
		FITOB_OUT_LEVEL3(verb(),"GridFactory WB set to extrap_combigrid (extrap_truncated_combigrid_WB) : 3" );
		special_diagonal_cut_off_level_ =
		XMLConfig->getDoubleConfiguration("thetaconfigurations.gridproperties.combigrid.COMBI_DIAGONAL_CUT_LEVEL.<xmlattr>.value");
		XMLConfig->getDoubleVectorConfiguration(
				"thetaconfigurations.gridproperties.combigrid.COMBI_DIAGONAL_CUT_LEVEL.<xmlattr>.vector",
                ',' , dimensionAdaptiveTruncationLevels_);
		gridtype_ = 18; typeOfGrid_ = GRID_WITHOUT_BOUNDARY;
	}


// ----------------- ADAPTIVE SGpp GRIDS ---------------------
	if (gridName == "sparsegrid_adaptive"){
		FITOB_OUT_LEVEL3(verb(),"GridFactory set to adaptive sparsegrid: 9" );
		gridtype_ = 9;
	}


// ------------------- FULL GRIDS -----------------------------
	if (gridName == "fullgrid"){
		FITOB_OUT_LEVEL3(verb(),"GridFactory set to full grid : 10" );
		gridtype_ = 10;
	}
	if (gridName == "fullgrid_WB"){
		FITOB_OUT_LEVEL3(verb(),"GridFactory WB set to full grid : 10" );
		gridtype_ = 10; typeOfGrid_ = GRID_WITHOUT_BOUNDARY;
	}

	// set the Opticom flag
	if ( ((gridtype_ > 0) && (gridtype_ < 10)) || ((gridtype_ > 15)) ){
		 // weather to use opticom or not
		 if ( "true" == XMLConfig->getStringConfiguration("thetaconfigurations.gridproperties.combigrid.USE_OPTICOM.<xmlattr>.value") ){
			 use_OptiCom_ = true;
		 } else {
			 use_OptiCom_ = false;
		 }
		 FITOB_OUT_LEVEL3(verb(),"GridFactory using OPTICOM : " << use_OptiCom_);
	}

	// weather to use scaling factors
	if ( "false" == XMLConfig->getStringConfiguration("thetaconfigurations.gridproperties.USE_SCALING.<xmlattr>.value") ){
		use_Scaling_ = false;
	} else {
		use_Scaling_ = true;
	}
	FITOB_OUT_LEVEL3(verb(),"GridFactory using SCALING : " << use_Scaling_);

}

boost::shared_ptr<MeshBase> GridFactory::createMesh(
		const Domain* gridDomain ,
		const FitobCalculator* calc ) const {

	bool useNewdomain = false;
	Domain* gridDomain_str = 0;

	if (!use_Scaling_){
		gridDomain_str = gridDomain->getNonScaledDomain();
		FITOB_OUT_LEVEL2( verb() , "NON SCALED DOMAIN:" << gridDomain_str->toString() );
		// we add to one container so that later will be deleted
		domainDepo_.push_back(gridDomain_str);
		useNewdomain = true;
	}
	// here we apply the dimension adaptive levels
	if (calc->dimensionAdaptiveComp()){
		useNewdomain = true;
		if (gridDomain_str == 0)
		{
			gridDomain_str = new Domain(gridDomain);
			// apply the lower dimension adaptivity
			gridDomain_str->applyDimAdaptLevels(calc->dimAdaptiveLevels());
			domainDepo_.push_back(gridDomain_str);
		}
		else
		{
			// apply the lower dimension adaptivity
			gridDomain_str->applyDimAdaptLevels(calc->dimAdaptiveLevels());
		}
	}

	// set the domain which will be used
	const Domain* gridDomain_tmp = (!useNewdomain) ? (gridDomain) : (gridDomain_str);

	// if there is no real dimensions then just create a constant mesh (the mesh is only a point)
	if (gridDomain_tmp->nrRealAxis() == 0) {
		// the constant grid will be set by default to 0.0
	   return( boost::shared_ptr<MeshBase>( new ConstantMesh(gridDomain_tmp, 0.0 ) ) );
	}
	// there are real dimensions so create the corresponding
	if ( (gridtype_ <= 4) && ( gridtype_ >= 1 ) /*&& ( gridDomain_tmp->nrRealAxis() > 1 )*/ ) {
        // create sparse grid
		FITOB_OUT_LEVEL3(verb(),"GridFactory::createMesh create SGppSparseGrid with domain:" << gridDomain_tmp->toString() );
#ifdef SGPP_DIRECT_SOLVER
		return( boost::shared_ptr<MeshBase>( new SGppSparseGrid( gridDomain_tmp , gridtype_ ,
				special_diagonal_cut_off_level_ , use_OptiCom_ , typeOfGrid_) ) );
#else
		FITOB_ERROR_EXIT(" use \"direct_solver=yes\" to use SGpp grids, in this configuration they are not available ");
		return( boost::shared_ptr<MeshBase>( new CombiGrid( gridDomain_tmp , gridtype_ ,
				special_diagonal_cut_off_level_ ,
				use_OptiCom_ , dimensionAdaptiveTruncationLevels_ , typeOfGrid_ ) ) );
#endif
	}
	if ( (gridtype_ >= 6) && ( gridtype_ < 9 ) && ( gridDomain_tmp->nrRealAxis() > 1 ) ) {
        // create sparse grid
		FITOB_OUT_LEVEL3(verb(),"GridFactory::createMesh create CombiGrid with domain:" << gridDomain_tmp->toString() );
		if (dimensionAdaptiveTruncationLevels_.size() < 1){
			// no dimension adaptive truncation
			return( boost::shared_ptr<MeshBase>( new CombiGrid( gridDomain_tmp , gridtype_ ,
					special_diagonal_cut_off_level_ ,
					use_OptiCom_ , dimensionAdaptiveTruncationLevels_ , typeOfGrid_ ) ) );
		} else {
			// dimension adaptive truncation
			return( boost::shared_ptr<MeshBase>( new CombiGrid( gridDomain_tmp , gridtype_ ,
					special_diagonal_cut_off_level_ ,
					use_OptiCom_ , dimensionAdaptiveTruncationLevels_ , typeOfGrid_ ) ) );
		}
	}
	if ( (gridtype_ >= 16) && ( gridtype_ < 19 ) && ( gridDomain_tmp->nrRealAxis() > 1 ) ) {
        // create sparse grid
		FITOB_OUT_LEVEL3(verb(),"GridFactory::createMesh create CombiGrid with domain:" << gridDomain_tmp->toString() );
		if (dimensionAdaptiveTruncationLevels_.size() < 1){
			// no dimension adaptive truncation
			return( boost::shared_ptr<MeshBase>( new ExtrapGrid( gridDomain_tmp , gridtype_ ,
					special_diagonal_cut_off_level_ ,
					use_OptiCom_ , dimensionAdaptiveTruncationLevels_ , typeOfGrid_ ) ) );
		} else {
			// dimension adaptive truncation
			return( boost::shared_ptr<MeshBase>( new ExtrapGrid( gridDomain_tmp , gridtype_ ,
					special_diagonal_cut_off_level_ ,
					use_OptiCom_ , dimensionAdaptiveTruncationLevels_ , typeOfGrid_ ) ) );
		}
	}
	if ( (gridtype_ == 9)) {
		// create adaptive sparse grid
		FITOB_OUT_LEVEL3(verb(),"GridFactory::createMesh create SGppSparseGridAdaptive with domain:" << gridDomain_tmp->toString() );
#ifdef SGPP_DIRECT_SOLVER
		return( boost::shared_ptr<MeshBase>( new SGppSparseGridAdaptive( gridDomain_tmp , gridtype_ ,
				special_diagonal_cut_off_level_ , use_OptiCom_ , typeOfGrid_, XMLConfig_) ) );
#else
		FITOB_ERROR_EXIT(" use \"direct_solver=yes\" to use SGpp adpative grids, in this configuration they are not available ");
		return( boost::shared_ptr<MeshBase>( new CombiGrid( gridDomain_tmp , gridtype_ ,
				special_diagonal_cut_off_level_ ,
				use_OptiCom_ , dimensionAdaptiveTruncationLevels_ , typeOfGrid_ ) ) );
#endif
	}
	if ( (gridtype_ == 10) || ( gridDomain_tmp->nrRealAxis() == 1 ) ) {
        // create full grid
		FITOB_OUT_LEVEL3(verb(),"GridFactory::createMesh create FullGrid with domain:" << gridDomain_tmp->toString() );
		if (typeOfGrid_ == GRID_WITH_BOUNDARY) {
			return( boost::shared_ptr<MeshBase>( new FullGrid( gridDomain_tmp ) ) );
		}
		if (typeOfGrid_ == GRID_WITHOUT_BOUNDARY) {
			return( boost::shared_ptr<MeshBase>( new FullGrid_WB( gridDomain_tmp ) ) );
		}
	}
	if ( gridDomain_tmp->nrRealAxis() <  1){
        // create constant grid (only one point)
		FITOB_OUT_LEVEL3(verb(),"GridFactory::createMesh create ConstantMesh with domain:" << gridDomain_tmp->toString() );
		return( boost::shared_ptr<MeshBase>( new ConstantMesh( gridDomain_tmp , 0.0) ) );
	}

	// Wall
	if (gridtype_ < 0){
		FITOB_ERROR_EXIT( "GridFactory::createMesh , negative configuration parameter");
	}
	FITOB_ERROR_EXIT( "GridFactory::createMesh , wrong grid name");
	return( boost::shared_ptr<MeshBase>( new ConstantMesh(gridDomain_tmp, 0.0 ) ) );
}

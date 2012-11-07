/*
 * FitobGridPlotter.cpp
 *
 *  Created on: Jul 7, 2010
 *      Author: benk
 */

#include "FitobGridPlotter.hpp"
#include "src/expressions/FitobExpressionBasis.hpp"
#include "src/scripteval/FitobCalculator.hpp"
#include <boost/lexical_cast.hpp>
//#include <boost/filesystem.hpp>

#include <iostream>
//#include <fstream>
//#include <sys/stat.h>
//#include <sys/types.h>

using namespace fitob;
using namespace std;

const string GridPlotter::outDir_ = "out";

GridPlotter::GridPlotter(const boost::shared_ptr<XMLConfiguration>& config)
 : couter_(0) , couter_gsp_(0) , res_x_(-1) , res_y_(-1) , activePlotting_(true) , makeGridSpecificPlotting_(false){
	setVerb(0);

	//  debugging flag to plot or not
	if ( config->getStringConfiguration("thetaconfigurations.gridproperties.DEBUG_PLOTTING.<xmlattr>.value") == "true" ){
		activePlotting_ = true;
	} else {
		activePlotting_ = false;
	}

	if ( config->getStringConfiguration("thetaconfigurations.gridproperties.DEBUG_PLOTTING.<xmlattr>.gridspecificplotting") == "true" ){
		makeGridSpecificPlotting_ = true;
	} else {
		makeGridSpecificPlotting_ = false;
	}

	res_x_ = config->getIntConfiguration("thetaconfigurations.gridproperties.DEBUG_PLOTTING.<xmlattr>.resolution");

	axisIndex0_ = config->getIntConfiguration("thetaconfigurations.gridproperties.DEBUG_PLOTTING.<xmlattr>.axis0");
	axisIndex0_ = (axisIndex0_ < 0) ? 0 : axisIndex0_;
	axisIndex1_ = config->getIntConfiguration("thetaconfigurations.gridproperties.DEBUG_PLOTTING.<xmlattr>.axis1");
	axisIndex1_ = (axisIndex1_ < 0) ? 1 : axisIndex1_;

	if ( res_x_ > 0 )
	{
		res_y_ = res_x_;
	} else {
		res_y_ = res_x_ = -1;
	}

	// create the "GridPlotter::outDir_" directory
	FITOB_OUT_LEVEL3(verb(),"GridPlotter create out directory");
	//__mode_t md = S_IRWXU | S_IRWXG | S_IRWXO;
	//mkdir( GridPlotter::outDir_.c_str() , md );
	//boost::filesystem::create_directories(GridPlotter::outDir_.c_str());
	//system("mkdir out");

	FITOB_OUT_LEVEL3(verb(),"GridPlotter created");
}

void GridPlotter::plot(const MeshContext* meshContext , const string& filename ,
		const ExpressionBasis* expr , const FitobCalculator* calc ){


	FITOB_OUT_LEVEL2(verb(),"GridPlotter::plot");

#if defined(FITOB_MPI)
	// if we are not on the 0-th processor then just exit ,since the corresponding grids just do not exist
	if (fitob::FITOB_MPI_Comm_rank() >= 1) {
		FITOB_OUT_LEVEL2(verb()," GridPlotter::plot , MPI END ");
		return ;
	}
#endif

	// test if we should do grid specific plotting
	if (makeGridSpecificPlotting_){
		plotGridSpecific( meshContext , filename );
	}

	// if plotting is not active then just return
	if (!activePlotting_) return;

	// todo: choose the correct plotting method based on the configuration

	// if the resolution is given then plot a cartesian grid (eqvidistant)
	if ( (res_y_ > 0 ) &&  ( res_x_ > 0 )){
		plotMAT( meshContext , filename , expr , calc );
	} else {
		plotMAT_grad( meshContext , filename , expr , calc);
	}

	// increment the counter
	couter_++;
}

void GridPlotter::plotVTK(const MeshContext* meshContext , const string& filename ,
		const ExpressionBasis* expr , const FitobCalculator* calc ){
  // todo: implement this
}

void GridPlotter::plotMAT(const MeshContext* meshContext , const string& filename ,
		const ExpressionBasis* expr , const FitobCalculator* calc ){

   FITOB_OUT_LEVEL2(verb(),"GridPlotter::plotMAT");

   //const string completeName = GridPlotter::outDir_+"/"+filename + boost::lexical_cast<std::string>(couter_)+".m";
   const string completeName = filename + boost::lexical_cast<std::string>(couter_)+".m";
   // it is important to get the domain of the mesh and not of the context
   // in case of dimension adaptivity these might be different
   const Domain& dom = meshContext->getMesh()->domain();
   DVector globalCoord = dom.getAverage();
   DVector result(0);

   // --- PLOTTING --- 1D
   if (dom.nrRealAxis() == 1)
   {
      double minX = dom.getGradedAxisMin(dom.localToGlobalIndex(axisIndex0_));
      double maxX = dom.getGradedAxisMax(dom.localToGlobalIndex(axisIndex0_));
      result.resize(res_x_);
      // loop and evaluate points
      for (int ii = 0 ; ii < res_x_ ; ii++){
    	  globalCoord[dom.localToGlobalIndex(axisIndex0_)] = minX + (maxX - minX) * (double(ii)/double(res_x_-1));
    	  result[ii] = meshContext->eval(globalCoord);
    	  if (expr != 0) {
        	  // todo: here could be different export variables
        	  globalCoord[dom.exportToGlobal(0)] = result[ii];
    		  result[ii] = expr->eval(globalCoord);
    	  }
      }

       // writing file
       ofstream myfile;
       myfile.open(completeName.c_str());
       myfile << "X = [ " << minX;
       for (int ii = 1 ; ii < res_x_ ; ii++){
    	   myfile << " , " << (minX + (maxX - minX) * (double(ii)/double(res_x_-1)));
       }
       myfile << "]; \n ";
       myfile << "res = [ " << result[0];
       for (int ii = 1 ; ii < res_x_ ; ii++) {
    	   myfile << " , " << result[ii];
       }
       myfile << "]; \n ";
       myfile << " plot(X,res); \n " << result[0];
       if ( calc != 0 ) {
    	   const boost::shared_ptr<ScriptModel> &scriptModel = calc->getScriptModel();
    	   myfile << "xlabel(\'" <<
    		(scriptModel->getVariable(dom.localToGlobalIndex(axisIndex0_)))->getVariableName() <<  "\'); \n ";
       }
       myfile.close();

   }

   // -- PLOTTING ---  >= 2D
   if (dom.nrRealAxis() >= 2)
   {
	  double minX = dom.getGradedAxisMin(dom.localToGlobalIndex(axisIndex0_));
	  double maxX = dom.getGradedAxisMax(dom.localToGlobalIndex(axisIndex0_));
	  double minY = dom.getGradedAxisMin(dom.localToGlobalIndex(axisIndex1_));
	  double maxY = dom.getGradedAxisMax(dom.localToGlobalIndex(axisIndex1_));
	  result.resize(res_x_*res_y_);
	  // loop and evaluate points
      for (int ii = 0 ; ii < res_x_ ; ii++){
    	  for (int jj = 0 ; jj < res_y_ ; jj++){
    	     globalCoord[dom.localToGlobalIndex(axisIndex0_)] = minX + (maxX - minX) * (double(ii)/double(res_x_-1));
    	     globalCoord[dom.localToGlobalIndex(axisIndex1_)] = minY + (maxY - minY) * (double(jj)/double(res_y_-1));
    	     result[ii*res_y_ + jj] = meshContext->eval(globalCoord);
       	     if (expr != 0) {
           	    // todo: here could be different export variables
           	    globalCoord[dom.exportToGlobal(0)] = result[ii*res_y_ + jj];
           	    result[ii*res_y_ + jj] = expr->eval(globalCoord);
       	     }
    	  }
      }

      ofstream myfile;
      myfile.open(completeName.c_str());
      myfile << "X = [ " << minX;
      for (int ii = 1 ; ii < res_x_ ; ii++){
   	     myfile << " , " << (minX + (maxX - minX) * (double(ii)/double(res_x_-1)));
      }
      myfile << "]; \n ";
      myfile << "Y = [ " << minY;
      for (int ii = 1 ; ii < res_y_ ; ii++){
   	     myfile << " , " << (minY + (maxY - minY) * (double(ii)/double(res_y_-1)));
      }
      myfile << "]; \n ";
      myfile << "res = [ " << result[0];
      for (int ii = 1 ; ii < res_x_*res_y_ ; ii++) {
       	  if ( (ii % res_y_ ) == 0){
        	  myfile << " ; " << result[ii];
          }
          else{
       	      myfile << " , " << result[ii];
          }
      }
      myfile << "]; \n ";
      myfile << "[x,y]=meshgrid(Y,X);\n";
      myfile << " surf(x,y,res); \n ";
      if ( calc != 0 ) {
   	   const boost::shared_ptr<ScriptModel> &scriptModel = calc->getScriptModel();
   	   myfile << "xlabel(\'" <<
   		(scriptModel->getVariable(dom.localToGlobalIndex(axisIndex1_)))->getVariableName() <<  "\'); \n ";
   	   myfile << "ylabel(\'" <<
   		(scriptModel->getVariable(dom.localToGlobalIndex(axisIndex0_)))->getVariableName() <<  "\'); \n ";
      }
      myfile.close();

   }
}


void GridPlotter::plotMAT_grad(const MeshContext* meshContext , const string& filename ,
		const ExpressionBasis* expr , const FitobCalculator* calc ){

   FITOB_OUT_LEVEL2(verb(),"GridPlotter::plotMAT_grad");

   //const string completeName = GridPlotter::outDir_+"/"+filename + boost::lexical_cast<std::string>(couter_)+".m";
   const string completeName = filename + boost::lexical_cast<std::string>(couter_)+".m";
   // it is important to get the domain of the mesh and not of the context
   // in case of dimension adaptivity these might be different
   const MeshBase* mesh = meshContext->getMesh();
   const Domain& dom = mesh->domain();
   DVector globalCoord = dom.getAverage();
   DVector result(0);

   // --- PLOTTING --- 1D
   if (dom.nrRealAxis() == 1)
   {
      result.resize( dom.getGradedAxis(dom.localToGlobalIndex(axisIndex0_)).size() );
      // loop and evaluate points
      for ( unsigned int ii = 0 ; ii < dom.getGradedAxis(dom.localToGlobalIndex(axisIndex0_)).size() ; ii++){
    	  globalCoord[ dom.localToGlobalIndex(axisIndex0_)] = dom.getGradedAxis(dom.localToGlobalIndex(axisIndex0_))[ii];
    	  result[ii] = meshContext->eval(globalCoord);
    	  if (expr != 0) {
        	  // todo: here could be different export variables
        	  globalCoord[dom.exportToGlobal(0)] = result[ii];
    		  result[ii] = expr->eval(globalCoord);
    	  }
      }

       // writing file
       ofstream myfile;
       myfile.open(completeName.c_str());
       myfile << "X = [ " << dom.getGradedAxis(dom.localToGlobalIndex(axisIndex0_))[0];
       for ( unsigned int ii = 1 ; ii < dom.getGradedAxis(dom.localToGlobalIndex(axisIndex0_)).size() ; ii++){
    	   myfile << " , " << dom.getGradedAxis(dom.localToGlobalIndex(axisIndex0_))[ii];
       }
       myfile << "]; \n ";
       myfile << "res = [ " << result[0];
       for ( unsigned int ii = 1 ; ii < dom.getGradedAxis(dom.localToGlobalIndex(axisIndex0_)).size() ; ii++) {
    	   myfile << " , " << result[ii];
       }
       myfile << "]; \n ";
       myfile << " plot(X,res); \n " << result[0];
       if ( calc != 0 ) {
    	   const boost::shared_ptr<ScriptModel> &scriptModel = calc->getScriptModel();
    	   myfile << "xlabel(\'" <<
    		(scriptModel->getVariable(dom.localToGlobalIndex(axisIndex0_)))->getVariableName() <<  "\'); \n ";
       }
       myfile.close();

   }

   // -- PLOTTING ---  >= 2D
   if (dom.nrRealAxis() >= 2)
   {
	  result.resize(dom.getGradedAxis(dom.localToGlobalIndex(axisIndex0_)).size() *
			  dom.getGradedAxis(dom.localToGlobalIndex(axisIndex1_)).size());
	  // loop and evaluate points
	  double res = 0.0;
      for ( unsigned int ii = 0 ; ii < dom.getGradedAxis(dom.localToGlobalIndex(axisIndex0_)).size() ; ii++){
    	  for ( unsigned int jj = 0 ; jj < dom.getGradedAxis(dom.localToGlobalIndex(axisIndex1_)).size() ; jj++){
    	     globalCoord[dom.localToGlobalIndex(axisIndex0_)] = dom.getGradedAxis(dom.localToGlobalIndex(axisIndex0_))[ii];
    	     globalCoord[dom.localToGlobalIndex(axisIndex1_)] = dom.getGradedAxis(dom.localToGlobalIndex(axisIndex1_))[jj];
     	     FITOB_OUT_LEVEL3(verb(),"GridPlotter::plotMAT_grad ii:"<<ii<<",i0: " << dom.localToGlobalIndex(axisIndex0_)
     	    		 << ", x0:" << globalCoord[dom.localToGlobalIndex(axisIndex0_)]);
     	     FITOB_OUT_LEVEL3(verb(),"GridPlotter::plotMAT_grad jj:"<<jj<<",i1: " << dom.localToGlobalIndex(axisIndex1_)
     	    		 << ", x1:" << globalCoord[dom.localToGlobalIndex(axisIndex1_)]);
    	     res= meshContext->eval(globalCoord);
    	     FITOB_OUT_LEVEL3(verb(),"GridPlotter::plotMAT_grad , result: " << res);
    	     result[ii*dom.getGradedAxis(dom.localToGlobalIndex(axisIndex1_)).size() + jj] = res;
       	     if (expr != 0) {
           	    // todo: here could be different export variables
           	    globalCoord[dom.exportToGlobal(0)] = result[ii*dom.getGradedAxis(dom.localToGlobalIndex(axisIndex1_)).size() + jj];
           	    result[ii*dom.getGradedAxis(dom.localToGlobalIndex(axisIndex1_)).size() + jj] = expr->eval(globalCoord);
       	     }
    	  }
      }

      ofstream myfile;
      myfile.open(completeName.c_str());
      myfile << "X = [ " << dom.getGradedAxis(dom.localToGlobalIndex(axisIndex0_))[0];
      for ( unsigned int ii = 1 ; ii < dom.getGradedAxis(dom.localToGlobalIndex(axisIndex0_)).size() ; ii++){
   	     myfile << " , " << dom.getGradedAxis(dom.localToGlobalIndex(axisIndex0_))[ii];
      }
      myfile << "]; \n ";
      myfile << "Y = [ " << dom.getGradedAxis(dom.localToGlobalIndex(axisIndex1_))[0];
      for ( unsigned int ii = 1 ; ii < dom.getGradedAxis(dom.localToGlobalIndex(axisIndex1_)).size() ; ii++){
   	     myfile << " , " << dom.getGradedAxis(dom.localToGlobalIndex(axisIndex1_))[ii];
      }
      myfile << "]; \n ";
      myfile << "res = [ " << result[0];
      for ( unsigned int ii = 1 ; ii < dom.getGradedAxis(dom.localToGlobalIndex(axisIndex0_)).size()*
                             dom.getGradedAxis(dom.localToGlobalIndex(axisIndex1_)).size() ; ii++) {
    	  if ( (ii % dom.getGradedAxis(dom.localToGlobalIndex(axisIndex1_)).size()) == 0){
    		  myfile << " ; " << result[ii];
    	  }
    	  else{
   	         myfile << " , " << result[ii];
    	  }
      }
      myfile << "]; \n ";
      myfile << "[x,y]=meshgrid(Y,X);\n";
      myfile << " surf(x,y,res); \n ";
      if ( calc != 0 ) {
   	   const boost::shared_ptr<ScriptModel> &scriptModel = calc->getScriptModel();
   	   myfile << "xlabel(\'" <<
   		(scriptModel->getVariable(dom.localToGlobalIndex(axisIndex1_)))->getVariableName() <<  "\'); \n ";
   	   myfile << "ylabel(\'" <<
   		(scriptModel->getVariable(dom.localToGlobalIndex(axisIndex0_)))->getVariableName() <<  "\'); \n ";
      }
      myfile.close();

   }
}


// grid specific plotting
void GridPlotter::plotGridSpecific(const MeshContext* meshContext , const string& filename){
	// here we create the grid specific plotting filename
	/*const string completeName = GridPlotter::outDir_+"/"+filename +"_gsp_"
			+boost::lexical_cast<std::string>(couter_gsp_);
	meshContext->getMesh()->gridSpecificPlot(completeName);
	couter_gsp_++;*/
}

void GridPlotter::plotGNU(const MeshContext* meshContext , const string& filename ,
		const ExpressionBasis* expr , const FitobCalculator* calc ){
  // todo: implement this
}


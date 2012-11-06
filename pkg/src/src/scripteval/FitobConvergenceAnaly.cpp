/*
 * FitobConvergenceAnaly.cpp
 *
 *  Created on: Aug 7, 2010
 *      Author: benk
 */

#include "src/scripteval/FitobConvergenceAnaly.hpp"
#include "src/utils/FitobXMLConfiguration.hpp"
#include "src/mesh/FitobFullGrid.hpp"
#include "src/evalcontext/FitobMeshContext.hpp"

#include <boost/lexical_cast.hpp>

using namespace fitob;
using namespace std;

ConvergenceAnaly::ConvergenceAnaly(
		const string& configurationXMLFile,
        const string& scriptFile,
        int lowLevel , int highLevel) :
  lowLevel_(lowLevel) ,  highLevel_(highLevel) ,
  confXMLFile_(configurationXMLFile) , sciprtFile_(scriptFile) , dimensionAdaptive_(false)
{
   setVerb(6);

   XMLConfiguration conf(confXMLFile_);

   FITOB_OUT_LEVEL3( verb() , "ConvergenceAnaly::ConvergenceAnaly , lowLevel:" << lowLevel_ << " , highLevel:" << highLevel_);
   // if we have wrong(default) input, then read them from the XML file
   if ( (lowLevel_ < 0) || (highLevel_ < 0)){
	   lowLevel_ = conf.getIntConfiguration("thetaconfigurations.gridproperties.convergence-tool.CONV_LOWEST_LEVEL.<xmlattr>.value");
	   highLevel_ = conf.getIntConfiguration("thetaconfigurations.gridproperties.convergence-tool.CONV_REF_LEVEL.<xmlattr>.value");
	   FITOB_OUT_LEVEL3( verb() , "ConvergenceAnaly::ConvergenceAnaly XML , lowLevel:" << lowLevel_ << " , highLevel:"<<highLevel_);
   }

	// if we have dimension adaptive convergence analysis
	if ("true" == conf.getStringConfiguration("thetaconfigurations.gridproperties.convergence-tool.<xmlattr>.dimAdaptive")){
        int nrRuns = conf.nrXMLNodes("thetaconfigurations.gridproperties.convergence-tool.CONV_DIM_ADAPT");
        lowLevel_ = 1; highLevel_ = nrRuns;
        dimensionAdaptive_ = true;
	}
}

void ConvergenceAnaly::doAnalysisConvergence(
		DVector& prices , DVector& conv,
		DVector& norm_L2 , DVector& norm_Inf,
		DVector& rate_L2 , DVector& rate_Linf, DVector& rate_pointw ){

	// resize the result vectors
	prices.resize(highLevel_-lowLevel_+1);
	conv.resize(highLevel_-lowLevel_+1);
	norm_L2.resize(highLevel_-lowLevel_+1);
	norm_Inf.resize(highLevel_-lowLevel_+1);
	rate_pointw.resize(highLevel_-lowLevel_+1);
	rate_L2.resize(highLevel_-lowLevel_+1);
	rate_Linf.resize(highLevel_-lowLevel_+1);

	std::vector< boost::shared_ptr<FullGrid> > fullgridarray( highLevel_-lowLevel_+1 );

	XMLConfiguration conf(confXMLFile_);

	for (int run = 0 ; run <= highLevel_-lowLevel_ ; run++){
		// create the calculator with the prescribed level
    	FitobCalculator calc( confXMLFile_ , sciprtFile_ , lowLevel_ + run , run);
    	// get the pice and the standard mesh
    	FITOB_OUT_LEVEL3( verb() , "ConvergenceAnaly::doAnalysisConvergence() , start calculating level:" << lowLevel_ + run );
    	fullgridarray[run] = calc.evaluateScript_FG( prices[run] );
    	FITOB_OUT_LEVEL3( verb() , "ConvergenceAnaly::doAnalysisConvergence() , calc level:" << lowLevel_ + run <<::std::setprecision( 12 ) << " , Price:" << prices[run]);


//    in case of MPI only rank 0 should do this
#if defined(FITOB_MPI)
	// evaluate the grid on the evaluation point
	if ( FITOB_MPI_Comm_rank() == 0 ) {
#endif
    	// calculate the convergence values
    	FITOB_OUT( " ====== ConvergenceAnaly::doAnalysisConvergence , printing results ====== " );
    	// for new output
    	FITOB_OUT( "  ");
		FITOB_OUT( "level, price, L-inf err, L2 err, L-inf rate, L2 rate, pointw. err, pointw. rate");
		FITOB_OUT( "-------------------------------------------------------------------------------");
    	for (int ii = 0 ; ii <= run ; ii++){
    		norm_L2[ii] = fitob::l2_norm( &(fullgridarray[ii]->unknVect()) , &(fullgridarray[run]->unknVect()) );
    		norm_Inf[ii] = fitob::inf_norm( &(fullgridarray[ii]->unknVect()) , &(fullgridarray[run]->unknVect()) );
    		conv[ii] = (prices[ii] - prices[run]) / prices[run];

    		rate_L2[ii] = log(norm_L2[0]/norm_L2[ii])/log(2.0)/ii;
    		rate_Linf[ii] = log(norm_Inf[0]/norm_Inf[ii])/log(2.0)/ii;
    		rate_pointw[ii] = log(fabs(conv[0]/conv[ii]))/log(2.0)/ii;


    		if (dimensionAdaptive_){
    			// todo: print the level vector
    			IVector dimensionAdaptiveLevels;
    			conf.getIntVectorConfiguration("thetaconfigurations.gridproperties.convergence-tool.CONV_DIM_ADAPT",
            	        ii , "<xmlattr>.level_vect" , ',' , dimensionAdaptiveLevels);
            	std::cout << "(";
            	for (int jj = 0 ; jj < (int)dimensionAdaptiveLevels.size() ; jj++){
            		std::cout << dimensionAdaptiveLevels[jj] << ",";
            	}
            	std::cout << ")";
    		}
    		// old output
    		//FITOB_OUT( "L"<<lowLevel_+ii<<" Price: " << prices[ii] << " , conv(rel):" << conv[ii] << " , L2:" << norm_L2[ii] << " , Inf:" << norm_Inf[ii] );
    		// new output
    		FITOB_OUT( (lowLevel_+ii) << ", " << prices[ii] << ", " << norm_Inf[ii] << ", " << norm_L2[ii] << ", " << rate_Linf[ii] << ", " << rate_L2[ii] << ", " << conv[ii] << ", " << rate_pointw[ii] );


    		// ===== plott the error of the grid ====
    		/*
    		FullGrid* fg = (fullgridarray[ii].get());
    		DVector copyVect;
    		copyVect.resize(fg->unknVect().size());
    		for (int tmp_ii = 0 ; tmp_ii < (int)fg->unknVect().size() ; tmp_ii++ ) { copyVect[tmp_ii] = fg->unknVect()[tmp_ii]; }

    		GridPlotter* plotter = calc.getPlotter().get();
    		Domain dom = Domain(calc.getNormEvalDomain());
    		MeshContext meshcontex = MeshContext( dom , 0.0 );
    		IVector specialLevels = dom.getLevelVector();
    		boost::shared_ptr<MeshBase> fgnew = boost::shared_ptr<MeshBase>(new FullGrid( &dom , specialLevels , &copyVect ));
    		meshcontex.setMesh( fgnew );
    		fitob::vect_diff( &copyVect , &(fullgridarray[run]->unknVect()));
    		string filename = "convergence_" + boost::lexical_cast<std::string>(ii) + "_";
    		// plot the error mesh
    		plotter->plot( &meshcontex , filename );*/
    	}
#if defined(FITOB_MPI)
	  } // the condition from the MPI so that this section is executed only by the 0-th rank
#endif
	}
	FITOB_OUT( " ====== ConvergenceAnaly::doAnalysisConvergence , END ====== " );
}

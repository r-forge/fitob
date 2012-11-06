/*
 * FitobPlotOperator.cpp
 *
 *  Created on: Aug 24, 2011
 *      Author: benk
 */

#include "FitobPlotOperator.hpp"
#include "src/mesh/FitobGridPlotter.hpp"
#include "src/scripteval/FitobCalculator.hpp"

using namespace fitob;


int PlotOperator::plotCounter_ = 0;


PlotOperator::PlotOperator(ExpressionBasis* plotExpression  )  :
	    		  OperatorBasis("PlotOperator" , (PlOperator) ) ,
 plotExpression_(plotExpression) {
}


PlotOperator::~PlotOperator() {
	// Nothing to do
}


void PlotOperator::backwardCalculation( boost::ptr_vector<MeshContext>& contextStack ,
		                          FitobCalculator* calc ,
			                      int& stackIndex ,
		                          double& timeStamp ){
	// get the plotter
	const boost::shared_ptr<GridPlotter> &plotter = calc->getPlotter();
	// plot the expression , regardless if it is turned off
	bool tmp_b = plotter->getActivePlotting();
	plotter->setActivePlotting(true);
	plotter->plot( &(contextStack[stackIndex]) , "plot" , plotExpression_ , calc );
	plotter->setActivePlotting(tmp_b);
	PlotOperator::plotCounter_++;
}

void PlotOperator::forwardMCSimulation( boost::ptr_vector<MeshContext>& contextStack ,
								  FitobCalculator* calc , MCMachine* MC ,
								  int& stackIndex , double& timeStamp ) {
	//todo: this could be implemented later
}

void PlotOperator::forwardMCEvaluation( boost::ptr_vector<MeshContext>& contextStack ,
								  FitobCalculator* calc , MCMachine* MC ,
								  int& stackIndex , double& timeStamp ) {
	//todo: this could be implemented later
}

void PlotOperator::backwardMCEvaluation( boost::ptr_vector<MeshContext>& contextStack ,
								  FitobCalculator* calc , MCMachine* MC ,
								  int& stackIndex , double& timeStamp ) {
	//todo: this could be implemented later
}

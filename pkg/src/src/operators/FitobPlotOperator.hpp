/* @HEADER@ */
// ************************************************************************
//
//                              Fitob
//            Copyright (2012) Janos Benk (benkjanos@gmail.com)
//
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Janos Benk (benkjanos@gmail.com),
//
// ************************************************************************
/* @HEADER@ */

/*
 * FitobPlotOperator.hpp
 *
 *  Created on: Aug 24, 2011
 *      Author: benk
 */

#ifndef FITOBPLOTOPERATOR_HPP_
#define FITOBPLOTOPERATOR_HPP_

#include "FitobOperatorBasis.hpp"
#include "src/expressions/FitobExpressionBasis.hpp"
#include "src/expressions/FitobConstantExpression.hpp"

namespace fitob {

using namespace std;

/** operator to plot one the actual mesh <br>
 * PLOT(P,0,1); -> this is by default the plotting <br>
 * but you can also do PLOT(MAX(S,P),1,2); in 3D case <br>
 * In 1D you also have to provide the third argument PLOT(P,0,0);*/
class PlotOperator : public fitob::OperatorBasis {
public:

	/** Ctor, which axes will be plotted depends on the Grid plotter object
	 * @param plotExpression the expression which should be plotted at each mesh point */
	PlotOperator( ExpressionBasis* plotExpression );

	virtual ~PlotOperator();

	/** function to do forward estimation */
	virtual void forwardEstimation( boost::ptr_vector<MeshContext>& contextStack ,
			                        FitobCalculator* calc ,
			                        int& stackIndex ,
			                        double& timeStamp ) { /* NOP */ }

	/** function to do backward calculation */
	virtual void backwardCalculation( boost::ptr_vector<MeshContext>& contextStack ,
			                          FitobCalculator* calc ,
				                      int& stackIndex ,
			                          double& timeStamp );

	/** function for simulate the scenarios forward */
	virtual void forwardMCSimulation( boost::ptr_vector<MeshContext>& contextStack ,
									  FitobCalculator* calc , MCMachine* MC ,
									  int& stackIndex , double& timeStamp );

	/** function for evaluate the scenarios forward */
	virtual void forwardMCEvaluation( boost::ptr_vector<MeshContext>& contextStack ,
									  FitobCalculator* calc , MCMachine* MC ,
									  int& stackIndex , double& timeStamp );

	/** function for evaluate the scenarios backward */
	virtual void backwardMCEvaluation( boost::ptr_vector<MeshContext>& contextStack ,
									  FitobCalculator* calc , MCMachine* MC ,
									  int& stackIndex , double& timeStamp );

	/** return the name of the operator */
	string toString() { return string(" Plot ");}

private:

	/** the expresion which should be plotted  on each grid point */
	ExpressionBasis* plotExpression_;

	static int plotCounter_;
};

}

#endif /* FITOBPLOTOPERATOR_HPP_ */

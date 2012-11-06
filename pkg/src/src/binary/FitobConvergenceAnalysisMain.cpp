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
 * FitobConvergenceAnalysisMain.cpp
 *
 *  Created on: Aug 18, 2010
 *      Author: benk
 */

#include "src/scripteval/FitobConvergenceAnaly.hpp"
#include "src/mesh/FitobFullGrid.hpp"
#include <iostream>

#include <R.h>
#include <Rdefines.h>

/** binary to study the convergence of the pricing method*/
double mainScriptPrice(int argc, char* argv1, char argv2, int level){

	// todo: implement this later
	/*
	// test if we have all the input parameters
	if (argc < 3){
		std::cout << " ERROR: Program needs at least two parameters , XMLFile and ScriptFile " << std::endl;
		return -1;
	}

    // MPI init in case
    fitob::FITOB_MPI_Init( &argc , NULL );
    std::string XMLFileName(argv1);
    std::string ScriptFileName(argv2);

    // todo: 
	fitob::ConvergenceAnaly conv_analysis(XMLFileName , ScriptFileName );
	DVector prices;
	DVector conv;
	DVector norm_L2;
	DVector norm_Inf;
	DVector rate_L2;
	DVector rate_Inf;
	DVector rate_pointw;

	// do the convergence analysis
	conv_analysis.doAnalysisConvergence(prices , conv , norm_L2 , norm_Inf, rate_L2, rate_Inf, rate_pointw );

   	// MPI finalize in case
	fitob::FITOB_MPI_Finalize( );*/
}

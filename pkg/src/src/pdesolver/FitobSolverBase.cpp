/*
 * FitobSolver.cpp
 *
 *  Created on: Jul 6, 2010
 *      Author: benk
 */

#include "FitobSolverBase.hpp"
#include "FitobMultigirdSolver.hpp"
#include "FitobDirectSolver.hpp"

using namespace fitob;
using namespace std;

SolverBase::SolverBase(const XMLConfiguration* config) : configuration_(config){

}

boost::shared_ptr<SolverBase> SolverBase::generateSolver(const XMLConfiguration* config){

	string solverType = config->getStringConfiguration("thetaconfigurations.solver.<xmlattr>.type");

	// return multigrid combi solver
	if ( solverType == "multigrid" )
	{
		return boost::shared_ptr<SolverBase>( new MultigirdSolver(config) );
	}
	// if the extern word is found then creat external (direct) solver
	if ( solverType.find("extern") != string::npos )
	{
		return boost::shared_ptr<SolverBase>( new DirectSolver(config) );
	}

	// throw error if no solver has been created
	FITOB_ERROR_TEST( false, " SolverBase::generateSolver , could not find matching solver")

	return boost::shared_ptr<SolverBase>( (SolverBase*)0 );
}


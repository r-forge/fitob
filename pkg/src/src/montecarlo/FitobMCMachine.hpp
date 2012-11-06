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
 * FitobMCMachine.hpp
 *
 *  Created on: Mar 14, 2011
 *      Author: benk
 */

#ifndef FITOBMCMACHINE_HPP_
#define FITOBMCMACHINE_HPP_

#include "src/utils/fitobdefs.hpp"
#include "src/montecarlo/FitobMCStep.hpp"
#include "src/utils/FitobXMLConfiguration.hpp"
#include "src/diffusionmodel/FitobModelCollection.hpp"
#include "src/montecarlo/FitobRandom.hpp"

using namespace std;

namespace fitob {

/** Class to contain the data for the Monte-Carlo simulation */
class MCMachine : public VerbClass {
public:

	MCMachine(const XMLConfiguration* xmlConfig ,
			  boost::shared_ptr<ModelCollection>& modelcollection );

	virtual ~MCMachine() {;}

	/** return the number of MC steps */
	inline int nrMCSteps() const { return nrMCSteps_; }

	/** add one Monte-Carlo step */
	void addMCStep(const Domain& domain , double Acttime , int evalContextIndex);

	/** return the requested Monte-Carlo step */
	MCStep& getMCStep(int stepI) { return MCSimulations_[stepI]; }

	const MCStep& getMCStep(int stepI) const { return MCSimulations_[stepI]; }

	/** nr of evaluation cycles*/
	int getNrEvaluation() const { return evaluationModes_.size(); }

	/** get the evaluation mode (1 forward , -1 backward) */
	int getEvaluationMode( int ind ) const { return evaluationModes_[ind]; }

private:

	/** makes the */
	void makeForwardMCStep( double actTime, double dt , MCStep& prevMCStep , MCStep& actMCStep );

	/** XMl configuration*/
	const XMLConfiguration* xmlConfig_;

	/** number of MC scenarios */
	int nrMCScenarios_;

	/** number of MC steps */
	int nrMCSteps_;

	/** the array of MC steps where all the simulations is stored */
	boost::ptr_vector<MCStep> MCSimulations_;

	/** the model collection object which contains all the diffusion models */
	boost::shared_ptr<ModelCollection> modelcollection_;

	/** random number generator*/
	Random randomNrGenerator_;

	/** the decomposed correlation matrix */
	DVector cholevskyMatrix_;

	/** vector of evaluation modes , usually there will be one, either
	 * forward or backward, but the user could specify more than that*/
	IVector evaluationModes_;

	/** the standard deviation of the initial diffusion (normal diffusion)*/
	double initialDiffusionSTD_;

	/** the prescribed time step */
	double maxDT_;

};

}

#endif /* FITOBMCMACHINE_HPP_ */

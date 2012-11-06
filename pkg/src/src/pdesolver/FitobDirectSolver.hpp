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
 * FitobDirectSolver.hpp
 *
 *  Created on: Aug 15, 2010
 *      Author: benk
 */

#ifndef FITOBDIRECTSOLVER_HPP_
#define FITOBDIRECTSOLVER_HPP_

#include "src/pdesolver/FitobSolverBase.hpp"

namespace fitob {

enum DirectSolverType{ ONLY_BS=0 , BS_HW=1 , NO_SOLVER = 2, BS_PARTIAL=3};


using namespace std;

/** Class for direct solving on the mesh ()*/
class DirectSolver : public SolverBase {
public:

	/** standard Ctor*/
	DirectSolver(const XMLConfiguration* config);

	virtual ~DirectSolver() {;}

	/** solves the PDE on the full mesh
	 * @param grid
	 * @param context
	 * @param models
	 * @param timeStep */
	void solvePDE(MeshBase *grid , const MeshContext* context ,
	              const ModelCollection* models ,
	              const FitobCalculator* calc ,
	              double timeStep);

private:

	/**
	 * calculate the theta value for the Hull White model
	 *
	 * @param a is the mean reversion rate
	 * @param sigma is the volatility
	 * @param T is the maturity
	 * @param t is the current time
	 *
	 * @return returns 0 if the file was successfully read, otherwise -1
	 */

	double calculatetheta4HW(double a, double sigma, double T, int t);

};

}

#endif /* FITOBDIRECTSOLVER_HPP_ */

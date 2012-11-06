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
 * FitobDirectSolver.cpp
 *
 *  Created on: Aug 15, 2010
 *      Author: benk
 */

#include "FitobDirectSolver.hpp"
#include "src/mesh/FitobSGppSparseGrid.hpp"
#include "src/diffusionmodel/FitobFactorModelHullWhite.hpp"
#include "src/diffusionmodel/FitobFactorModelGeomBr.hpp"
#include "src/diffusionmodel/FitobFactorModelNormBr.hpp"
#include <boost/lexical_cast.hpp>

//#define SGPP_DIRECT_SOLVER

using namespace fitob;
using namespace std;


#define CRNIC_IMEUL_STEPS 3

DirectSolver::DirectSolver(const XMLConfiguration* config) : SolverBase(config) {

	setVerb(6);
}


void DirectSolver::solvePDE(
		      MeshBase *grid , const MeshContext* context ,
              const ModelCollection* models ,
              const FitobCalculator* calc ,
              double timeStep){
}


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

double DirectSolver::calculatetheta4HW(double a, double sigma, double T, int t)
{
	double theta=0;
	return theta=0.04*a + pow(sigma,2.0)*(1-exp(-2*a*(T-t)))/(2*a);
}


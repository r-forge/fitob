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
 * FitobSolver.hpp
 *
 *  Created on: Jul 6, 2010
 *      Author: benk
 */

#ifndef FITOBSOLVER_HPP_
#define FITOBSOLVER_HPP_

#include "src/utils/fitobdefs.hpp"
#include "src/diffusionmodel/FitobModelCollection.hpp"
#include "src/evalcontext/FitobMeshContext.hpp"
#include "src/utils/FitobXMLConfiguration.hpp"

namespace fitob {

  // forward declaration of this class
  class FitobCalculator;

  using namespace std;

  class SolverBase : public VerbClass{
  public:

	  // input is the XML config
	SolverBase(const XMLConfiguration* config);

	virtual ~SolverBase() {;}

	/** solves the PDE on the full mesh
	 * @param grid
	 * @param context
	 * @param models
	 * @param timeStep */
	virtual void solvePDE(MeshBase *grid , const MeshContext* context ,
			              const ModelCollection* models ,
			              const FitobCalculator* calc ,
			              double timeStep) = 0;

	/** get the configuration object */
	inline const XMLConfiguration* conf() const { return configuration_; }

	/** Factory model and generate Solver based on the configuration */
	static boost::shared_ptr<SolverBase> generateSolver(const XMLConfiguration* config);

  private:

	/** pointer to configuration object */
	const XMLConfiguration* configuration_;

};

}

#endif /* FITOBSOLVER_HPP_ */

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
 * FitobRegularization.hpp
 *
 *  Created on: Jul 7, 2011
 *      Author: benk
 */

#ifndef FITOBREGULARIZATION_HPP_
#define FITOBREGULARIZATION_HPP_

#include  "src/utils/fitobdefs.hpp"
#include  "src/utils/FitobXMLConfiguration.hpp"
#include  "src/expressions/FitobExpectedExpression.hpp"


namespace combigrid{
   class GridDomain;
}

namespace sg{
	namespace datadriven{
		class LearnerBase;
		class LearnerBaseSP;
	}
}

namespace fitob {

  using namespace std;

  /** forward declaration of the Fitob calculator*/
  class MeshBase;

/** The class which */
class Regularization {

public:

	/** Ctor */
	Regularization( const FitobCalculator* calc , DVector& XCoords , DVector& YCoords ,
			const Domain& dom , const ExpectedExpression* expectExpr);

	/** Dtor delete the created objects (Mesh)*/
	virtual ~Regularization();

	/** return the mesh which has been generated */
	inline const MeshBase* getMesh() const { return actualMesh_.get(); }

private:

	/** special function */
	void configureSGppSolver(const FitobCalculator* calc , DVector& XCoords , DVector& YCoords , const Domain& dom );

	/** the mesh on which the regularization (Tikhonov) will be done */
	boost::shared_ptr<MeshBase> actualMesh_;

	/** domain of the grid */
	boost::shared_ptr<Domain> fitobDomain_;

	/** The expected expression which trigered the regression*/
	const ExpectedExpression* expectedExpression_;

	combigrid::GridDomain* domain_;
#ifdef COMBI_REGRESSION_SOLVER
	sg::datadriven::LearnerBase* myLearner_;

	sg::datadriven::LearnerBaseSP* myLearnerSP_;
#endif
};

}

#endif /* FITOBREGULARIZATION_HPP_ */

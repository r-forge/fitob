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
 * FitobDaOperatorEvaluable.hpp
 *
 *  Created on: Jul 7, 2010
 *      Author: benk
 */

#ifndef FITOBDAOPERATOREVALUABLE_HPP_
#define FITOBDAOPERATOREVALUABLE_HPP_

#include "src/evalcontext/FitobEvaluable.hpp"
#include "src/evalcontext/FitobMeshContext.hpp"
#include "src/operators/FitobDaOperator.hpp"

namespace fitob{

  using namespace std;

  /** The base class for all grids (mesh) */
  class DaOperatorEvaluable : public Evaluable{
  public:

	 DaOperatorEvaluable(const MeshContext* meshcontext
	   , const DaOperator* daop) : meshcontext_(meshcontext) , daop_(daop){
	 }

	/** see super class */
	double eval(const DVector& globalCoordonates) const {
		DVector tmp_vect_ = globalCoordonates;
		// first shift the global coordinates with the Da operator
		tmp_vect_[daop_->getVariable()->getGlobalIndex()] = daop_->getExpression()->eval(tmp_vect_);
		// return the evaluation result of the grid
		return meshcontext_->eval(tmp_vect_);
	}

	/** see super class */
	void eval(const std::vector<DVector>& globalCoordonates , DVector& resVector) const {
	   std::vector<DVector> globalCoordonates_tmp;
	   DVector tmp_res(globalCoordonates.size());
       // first copy the coordinates
	   for (unsigned int ii = 0 ; ii < globalCoordonates.size() ; ii++){
		   globalCoordonates_tmp[ii] = globalCoordonates[ii];
	   }
	   // eval the da expression
	   daop_->getExpression()->eval(globalCoordonates_tmp,tmp_res);
	   // set the global coordinates accordingly
	   for (unsigned int ii = 0 ; ii < globalCoordonates.size() ; ii++){
		   globalCoordonates_tmp[ii][daop_->getVariable()->getGlobalIndex()] = tmp_res[ii];
	   }
       // eval the mesh on those points
	   meshcontext_->eval(globalCoordonates_tmp,resVector);
	}

  private:

	/** Mesh context */
	const MeshContext* meshcontext_;

	/** Da operator (used only to get its expression )*/
	const DaOperator* daop_;

	/** used to have a permanent vector*/
	//DVector tmp_vect_;
  };
}

#endif /* FITOBDAOPERATOREVALUABLE_HPP_ */

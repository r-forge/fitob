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
 * FitobPowExpression.hpp
 *
 *  Created on: Feb 27, 2010
 *      Author: benk
 */

#ifndef FITOBPOWEXPRESSION_HPP_
#define FITOBPOWEXPRESSION_HPP_

#include "FitobExpressionBasis.hpp"

namespace fitob{

  using namespace std;

  /** */
  class PowExpression: public fitob::ExpressionBasis {
  public:

	 /** */
	  PowExpression(ExpressionBasis* left , ExpressionBasis* right);

	 /** */
	 virtual ~PowExpression();

	/** see superclass*/
	inline bool isConstantExpression(const MeshContext& context) const {
		return (left_->isConstantExpression(context) && right_->isConstantExpression(context));
	}

	/** see superclass*/
	inline double eval(const DVector& globalCoordonates) const {
		return pow(left_->eval(globalCoordonates),right_->eval(globalCoordonates));
	}

	/** see superclass*/
	void eval(const std::vector<DVector>& globalCoordonates , DVector& resVector) const {
		DVector left_res(resVector.size());
		left_->eval(globalCoordonates,left_res);
		right_->eval(globalCoordonates,resVector);
	    for (unsigned int ii = 0 ; ii < resVector.size() ; ii++){
	    	resVector[ii] = pow(left_res[ii] , resVector[ii]);
	    }
	}

    inline string toString() const { return "PowExpression(" +
    	left_->toString() + ") ^ ("+ right_->toString() + ")"; }

  private:

	 /** */
	 ExpressionBasis* left_;

	 /** */
	 ExpressionBasis* right_;
  };

}

#endif /* FITOBPOWEXPRESSION_HPP_ */

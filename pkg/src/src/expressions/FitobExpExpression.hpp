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
 * FitobExpExpression.hpp
 *
 *  Created on: Feb 27, 2010
 *      Author: benk
 */

#ifndef FITOBEXPEXPRESSION_HPP_
#define FITOBEXPEXPRESSION_HPP_

#include "FitobExpressionBasis.hpp"
#include <math.h>

namespace fitob{

  using namespace std;

  /** */
  class ExpExpression: public fitob::ExpressionBasis {
  public:

	 /** */
	  ExpExpression(ExpressionBasis* expr);

	 /** */
	 virtual ~ExpExpression();

	/** see superclass*/
	inline bool isConstantExpression(const MeshContext& context) const {
		return exp_->isConstantExpression(context);
	}

	/** see superclass*/
	inline double eval(const DVector& globalCoordonates) const {
		return exp(exp_->eval(globalCoordonates));
	}

	/** see superclass*/
	void eval(const std::vector<DVector>& globalCoordonates , DVector& resVector) const {
	   exp_->eval(globalCoordonates , resVector );
	   for (unsigned int ii = 0 ; ii < resVector.size() ; ii++){
		   resVector[ii] = exp(resVector[ii]);
	   }
	}

    inline string toString() const {
    	return "ExpExpression(" +
    	exp_->toString() + ") ";
    }

  private:

	 /** */
	 ExpressionBasis* exp_;
  };

}

#endif /* FITOBEXPEXPRESSION_HPP_ */

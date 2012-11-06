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
 * FitobSqrtExpression.hpp
 *
 *  Created on: Feb 27, 2010
 *      Author: benk
 */

#ifndef FITOBSQRTEXPRESSION_HPP_
#define FITOBSQRTEXPRESSION_HPP_

#include "FitobExpressionBasis.hpp"

namespace fitob{

  using namespace std;

  /** */
  class SqrtExpression: public fitob::ExpressionBasis {
  public:

	 /** */
	  SqrtExpression(ExpressionBasis* expr);

	 /** */
	 virtual ~SqrtExpression();

	/** see superclass */
	inline bool isConstantExpression(const MeshContext& context) const {
		return exp_->isConstantExpression(context);
	}

	/** see superclass */
	inline double eval(const DVector& globalCoordonates) const {
		const double val = exp_->eval(globalCoordonates);
		FITOB_ERROR_TEST(val >= 0 , " SqrtExpression::eval , negative value ");
		return sqrt(val);
	}

	/** see superclass */
	void eval(const std::vector<DVector>& globalCoordonates , DVector& resVector) const {
	   exp_->eval(globalCoordonates,resVector);
	   for (unsigned int ii = 0 ; ii < resVector.size() ; ii++){
			FITOB_ERROR_TEST(resVector[ii] >= 0 , " SqrtExpression::eval 2, negative value ");
			resVector[ii] = sqrt(resVector[ii]);
	   }
	}

    inline string toString() const { return "SqrtExpression(" + exp_->toString() + ")"; }

  private:

	 /** */
	 ExpressionBasis* exp_;
  };

}

#endif /* FITOBSQRTEXPRESSION_HPP_ */

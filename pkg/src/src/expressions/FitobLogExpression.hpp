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
 * FitobLogExpression.hpp
 *
 *  Created on: Feb 27, 2010
 *      Author: benk
 */

#ifndef FITOBLOGEXPRESSION_HPP_
#define FITOBLOGEXPRESSION_HPP_

#include "FitobExpressionBasis.hpp"

namespace fitob{

  using namespace std;

  /** */
  class LogExpression: public fitob::ExpressionBasis {
  public:

	 /** */
	  LogExpression(ExpressionBasis* expr);

	 /** */
	 virtual ~LogExpression();

	/** see superclass*/
	inline bool isConstantExpression(const MeshContext& context) const {
		return exp_->isConstantExpression(context);
	}

	/** see superclass*/
	inline double eval(const DVector& globalCoordonates) const {
		double val = exp_->eval(globalCoordonates);
		FITOB_ERROR_TEST(val >= 0 , " LogExpression:eval , negative value");
		return log(val);
	}

	/** see superclass*/
	void eval(const std::vector<DVector>& globalCoordonates , DVector& resVector) const {
		exp_->eval(globalCoordonates, resVector );
	    for (unsigned int ii = 0 ; ii < resVector.size() ; ii++){
	    	FITOB_ERROR_TEST( resVector[ii] >= 0 , " LogExpression:eval 2 , negative value");
	    	resVector[ii] = log(resVector[ii]);
	    }
	}

    inline string toString() const { return "LogExpression(" + exp_->toString() + ")"; }

  private:

	 /** */
	 ExpressionBasis* exp_;
  };

}

#endif /* FITOBLOGEXPRESSION_HPP_ */

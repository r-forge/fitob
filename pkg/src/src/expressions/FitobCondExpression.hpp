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
 * FitobCondExpression.hpp
 *
 *  Created on: Aug 23, 2011
 *      Author: benk
 */

#ifndef FITOBCONDEXPRESSION_HPP_
#define FITOBCONDEXPRESSION_HPP_

#include "FitobExpressionBasis.hpp"

namespace fitob {

/** this expression corresponds to the  () ? : construct in C/C++ <br>
 * the fist expression is a boolean expression, if it is true then the value of
 * the first expression will be returned, otherwise the value of the second
 * expression will be returned. */
class CondExpression : public fitob::ExpressionBasis {
public:

	CondExpression(ExpressionBasis* cond , ExpressionBasis* trueExpr , ExpressionBasis* falseExpr);

	virtual ~CondExpression();

	/** see superclass*/
	inline bool isConstantExpression(const MeshContext& context) const {
		return (cond_->isConstantExpression(context) && trueExpr_->isConstantExpression(context)
				&& falseExpr_->isConstantExpression(context));
	}

	/** see superclass*/
	inline double eval(const DVector& globalCoordonates) const {
		return (cond_->eval(globalCoordonates) > FITOB_NUMERICALZERO )?
				trueExpr_->eval(globalCoordonates) : falseExpr_->eval(globalCoordonates);
	}

	/** see superclass*/
	void eval(const std::vector<DVector>& globalCoordonates , DVector& resVector) const {
		DVector cond_res(resVector.size());
		DVector true_res(resVector.size());
		cond_->eval(globalCoordonates,cond_res);
		trueExpr_->eval(globalCoordonates,true_res);
		falseExpr_->eval(globalCoordonates,resVector);
	    for (unsigned int ii = 0 ; ii < resVector.size() ; ii++){
	    	resVector[ii] = (cond_res[ii] > FITOB_NUMERICALZERO )? true_res[ii] : resVector[ii];
	    }
	}

    inline string toString() const { return "COND(" +
    	cond_->toString() + ") ? ("+ trueExpr_->toString() + ") : ("+ falseExpr_->toString() + ")"; }
private:

    /** condition expression */
	ExpressionBasis* cond_;
    /** expression on the true branch*/
	ExpressionBasis* trueExpr_;
    /** expression on the false branch*/
	ExpressionBasis* falseExpr_;

};

}

#endif /* FITOBCONDEXPRESSION_HPP_ */

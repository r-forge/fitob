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
 * FitobIfElseAssigmentExpression.hpp
 *
 *  Created on: Jul 20, 2010
 *      Author: benk
 */

#ifndef FITOBIFELSEASSIGMENTEXPRESSION_HPP_
#define FITOBIFELSEASSIGMENTEXPRESSION_HPP_

#include "FitobExpressionBasis.hpp"

namespace fitob {

  using namespace std;


class IfElseAssigmentExpression : public fitob::ExpressionBasis {
public:

	/** */
	IfElseAssigmentExpression( const ExpressionBasis* cond ,
			                        const ExpressionBasis* trueB ,
			                        const ExpressionBasis* falseB );

	virtual ~IfElseAssigmentExpression() {;}

	 /** see superclass */
	inline bool isConstantExpression(const MeshContext& context) const {
		return false;
	}

	/** see superclass */
	inline double eval(const DVector& globalCoordonates) const {
		// evaluate the expression
		double val = cond_->eval(globalCoordonates);
		// if true then evaluate the true branch
		if (val > FITOB_NUMERICALZERO)
		   return true_->eval(globalCoordonates);
		else{
			val = globalCoordonates[1];
		   if (hasElseBanch_) val = false_->eval(globalCoordonates);
           // return either false branch value or the mesh value
		   return val;
		}
	}

	/** see superclass */
	void eval(const std::vector<DVector>& globalCoordonates , DVector& resVector) const {
		DVector cond_res(resVector.size());
		DVector false_res(resVector.size() , globalCoordonates[0][1]);
		cond_->eval(globalCoordonates,cond_res);
		true_->eval(globalCoordonates,resVector);
		if (hasElseBanch_) false_->eval(globalCoordonates,false_res);
	    for (unsigned int ii = 0 ; ii < resVector.size() ; ii++){
	    	if (cond_res[ii] > FITOB_NUMERICALZERO){
	    		//resVector[ii] = resVector[ii]
	    	} else {
	    	   resVector[ii] = false_res[ii];
	    	}
	    }
	}


   inline string toString() const {
	   return "IfElseAssigmentExpression( cond: );"; /* + cond_->toString() +
   	  + " TRUE: "+ true_->toString() + "  ELSE:"+ hasElseBanch_ + ")";*/ }

 private:

    /** */
    const ExpressionBasis* cond_;

	 /** */
    const ExpressionBasis* true_;

	 /** */
    const ExpressionBasis* false_;

    /** */
    bool hasElseBanch_;
};

}

#endif /* FITOBIFELSEASSIGMENTEXPRESSION_HPP_ */

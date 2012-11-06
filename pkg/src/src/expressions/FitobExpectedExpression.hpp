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

#ifndef FITOBEXPECTEDEXPRESSION_HPP_
#define FITOBEXPECTEDEXPRESSION_HPP_

#include "FitobExpressionBasis.hpp"

namespace fitob{

  using namespace std;

  /** The expected value of one expression. <br>
   *  In case of PDE evaluation this is just the coordinate itself. <br>
   *  In case of MC context this will be some different algorithm . <br>
   *  If we need to use the regression then if has to be the form EXPECT(V!) where export V;
   *  EXPECT(1+V), should also work, and 1+EXPECT(V!) as well*/
  class ExpectedExpression: public fitob::ExpressionBasis {
  public:

	 /** */
	  ExpectedExpression(ExpressionBasis* expr);

	 /** */
	 virtual ~ExpectedExpression();

	/** see superclass*/
	inline bool isConstantExpression(const MeshContext& context) const {
		return exp_->isConstantExpression(context);
	}

	/** see superclass*/
	double eval(const DVector& globalCoordonates) const;

	/** see superclass*/
	void eval(const std::vector<DVector>& globalCoordonates , DVector& resVector) const;

    inline string toString() const{
    	return "ExpectedExpression(" +
    	exp_->toString() + ") ";
    }

    /** */
    inline void setITriggeredRegression(bool b) const { ITriggeredRegression_ = b; }

    /** */
    inline bool getITriggeredRegression() const { return ITriggeredRegression_; }

  private:

	 /** */
	 ExpressionBasis* exp_;

	 /** */
	 mutable bool ITriggeredRegression_;
  };

}

#endif /* FITOBEXPECTEDEXPRESSION_HPP_ */

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

#ifndef FITOBDISCOUNTEXPRESSION_HPP_
#define FITOBDISCOUNTEXPRESSION_HPP_

#include "FitobExpressionBasis.hpp"
#include <math.h>

namespace fitob{

  using namespace std;

  /** Class to represent the discount factor of the context */
  class DiscountExpression: public fitob::ExpressionBasis {
  public:

	/** */
	DiscountExpression(ExpressionBasis* t1, ExpressionBasis* t2);

	/** */
	virtual ~DiscountExpression();

	/** see superclass*/
	bool isConstantExpression(const MeshContext& context) const;

	/** see superclass*/
	double eval(const DVector& globalCoordonates) const;

	/** see superclass*/
	void eval(const std::vector<DVector>& globalCoordonates , DVector& resVector) const;

    inline string toString() const {
    	return "DiscountExpression( " +
    	t1_->toString() + "," + t2_->toString() + " ) ";
    }

  private:

	 /** first tenior */
	 ExpressionBasis* t1_;

	 /** second tenior */
	 ExpressionBasis* t2_;
  };

}

#endif /* FITOBDISCOUNTEXPRESSION_HPP_ */

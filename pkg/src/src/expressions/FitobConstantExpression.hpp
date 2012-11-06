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
 * FitobConstantExpression.hpp
 *
 *  Created on: Feb 26, 2010
 *      Author: benk
 */

#ifndef FITOBCONSTANTEXPRESSION_HPP_
#define FITOBCONSTANTEXPRESSION_HPP_

#include "src/expressions/FitobExpressionBasis.hpp"
#include "src/utils/fitobdefs.hpp"
#include <boost/lexical_cast.hpp>

namespace fitob{

  using namespace std;

  class ConstantExpression : public ExpressionBasis{
  public:

	/** */
	ConstantExpression(const double value);

	/** */
	virtual ~ConstantExpression();

	/** see superclass */
    inline bool isConstantExpression(const MeshContext& context) const { return true;}

	/** see superclass */
    inline double eval(const DVector& globalCoordonates) const { return constValue_; }

	/** see superclass */
	void eval(const std::vector<DVector>& globalCoordonates , DVector& resVector) const {
	   for (unsigned int ii = 0 ; ii < resVector.size() ; ii++) resVector[ii] = constValue_;
	}

    /** */
    inline string toString() const {
    	return "ConstantExpression :" +
    	boost::lexical_cast<std::string>(constValue_);
    }

    /** returns the constant value */
    inline double getConstValue() const { return constValue_; }

  private:

	/** constant value*/
	const double constValue_;

  };
}
#endif /* FITOBCONSTANTEXPRESSION_HPP_ */

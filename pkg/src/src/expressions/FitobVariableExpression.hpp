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
 * FitobVariableExpression.hpp
 *
 *  Created on: Feb 27, 2010
 *      Author: benk
 */

#ifndef FITOBVARIABLEEXPRESSION_HPP_
#define FITOBVARIABLEEXPRESSION_HPP_

#include <boost/lexical_cast.hpp>

#include "src/expressions/FitobExpressionBasis.hpp"
#include "src/evalcontext/FitobMeshContext.hpp"
#include "src/variables/FitobVariable.hpp"
#include "src/utils/fitobdefs.hpp"

namespace fitob{

  using namespace std;

 class VariableExpression: public ExpressionBasis {
 public:

	VariableExpression(const Variable* variable, double signVar, bool MClookback = false);

	virtual ~VariableExpression();

	/** see superclass*/
    inline bool isConstantExpression(const MeshContext& context) const {
    	//std::cout << "VariableExpression::isConstantExpression: context.getDom().getAxisLevel(variable_->getGlobalIndex(): "
    	//		<< context.getDom().getAxisLevel(variable_->getGlobalIndex()) << " bool:"
    	//		<< (bool)(context.getDom().getAxisLevel(variable_->getGlobalIndex()) == 0) << std::endl;
    	return (bool)(context.getDom().getAxisLevel(variable_->getGlobalIndex()) == 0);
    }

	/** see superclass*/
    inline double eval(const DVector& globalCoordonates) const {
    	return sign_ * globalCoordonates[variable_->getGlobalIndex()];}

	/** see superclass*/
	void eval(const std::vector<DVector>& globalCoordonates , DVector& resVector) const {
	   for (unsigned int ii = 0 ; ii < resVector.size() ; ii++){
		   resVector[ii] = (sign_ * globalCoordonates[ii][(variable_->getGlobalIndex())]);
	   }
	}

    /** */
    inline string toString() const{
    	return "VariableExpression: " +
    	boost::lexical_cast<std::string>(sign_) + " * " +
    	variable_->getVariableName();
    }

	/** return the variable */
	inline const Variable* getVariable() const { return variable_; }

	/** In case of MC we are interested if the variable is from the previous time step or not <br>
	 *  This is important to know if we have to take the values from the previous time step*/
	inline bool isVariableLookingBack() const { return MClookback_; }

 private:

	/** The variable where this is the end */
	const Variable* variable_;

	/** the sign of the constant */
	double sign_;

	/** if this variable in MC context has to look back e.g. V = S! */
	bool MClookback_;
 };
}

#endif /* FITOBVARIABLEEXPRESSION_HPP_ */

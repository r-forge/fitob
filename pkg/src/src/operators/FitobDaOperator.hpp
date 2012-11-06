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
 * FitobDaOperator.hpp
 *
 *  Created on: Feb 28, 2010
 *      Author: benk
 */

#ifndef FITOBDAOPERATOR_HPP_
#define FITOBDAOPERATOR_HPP_

#include "FitobOperatorBasis.hpp"
#include "src/utils/fitobdefs.hpp"
#include "src/variables/FitobVariable.hpp"
#include "src/expressions/FitobExpressionBasis.hpp"
#include "src/utils/FitobXMLConfiguration.hpp"

#include <boost/lexical_cast.hpp>

namespace fitob{

  using namespace std;

  class FitobCalculator;

  /** The Da operation , shift for axis variables or for constants */
  class DaOperator: public fitob::OperatorBasis {

  public:

	/** Ctor */
	DaOperator( const boost::shared_ptr<XMLConfiguration>& xmlconfiguration ,
			    Variable* assignedVariable, ExpressionBasis* expression);

	virtual ~DaOperator();

	/** function to do forward estimation */
	virtual void forwardEstimation( boost::ptr_vector<MeshContext>& contextStack ,
			                        FitobCalculator* calc ,
			                        int& stackIndex ,
			                        double& timeStamp );

	/** function for diffusion axis size (enlargement) estimation */
	virtual void forwardEstimation_DiffusionEnlarement(
			boost::ptr_vector<MeshContext>& contextStack ,
			                          FitobCalculator* calc ,
			                          DVector& gradValues,
			                          int globalVariableIndex,
				                      int& stackIndex ,
			                          double timeStamp ) const;

	/** function to do backward calculation */
	virtual void backwardCalculation( boost::ptr_vector<MeshContext>& contextStack ,
			                          FitobCalculator* calc ,
				                      int& stackIndex ,
			                          double& timeStamp );

	/** function for simulate the scenarios forward */
	virtual void forwardMCSimulation( boost::ptr_vector<MeshContext>& contextStack ,
									  FitobCalculator* calc , MCMachine* MC ,
									  int& stackIndex , double& timeStamp ) ;

	/** function for evaluate the scenarios forward */
	virtual void forwardMCEvaluation( boost::ptr_vector<MeshContext>& contextStack ,
									  FitobCalculator* calc , MCMachine* MC ,
									  int& stackIndex , double& timeStamp ) ;

	/** function for evaluate the scenarios backward */
	virtual void backwardMCEvaluation( boost::ptr_vector<MeshContext>& contextStack ,
									  FitobCalculator* calc , MCMachine* MC ,
									  int& stackIndex , double& timeStamp ) ;


	const Variable* getVariable() const {  return assignedVariable_; }

	const ExpressionBasis* getExpression() const {  return expression_; }

	string toString() { return operationName_ + " Var:" +
		assignedVariable_->getVariableName() + " GlbI:" + boost::lexical_cast<std::string>(assignedVariable_->getGlobalIndex()) +
		" Exp:"+expression_->toString();}

  private:

	/** */
	void doShift_forward(DVector& axisVector , const Domain& actDom , double actTime ) const;

	/** recursive function to determine the maximal size */
	void calcForward_rec(int index, const Domain& actDom ,
			double& minDaA ,
			double& maxDaA ,
			DVector& evalCoord,
			IVector& minP,
			IVector& maxP ) const ;

	/** Pointer to the variable which should be shifted*/
	Variable*  assignedVariable_;

	/** expression which shifts one variable */
	ExpressionBasis* expression_;

	/** If we have a DA operation with a constant expression, that is just a simple
	 * value assigment to a constant import variable, in that case the evaluation is much simpler*/
	bool expressionIsConstant_;

	/** flag to indicate weather we always estimate the new domain (in forward estimation)
	 * */
	bool alwaysIncreasingDa_;
  };
}

#endif /* FITOBDAOPERATOR_HPP_ */

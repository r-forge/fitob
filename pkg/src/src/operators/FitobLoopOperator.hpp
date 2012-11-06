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
 * FitobLoopOperator.hpp
 *
 *  Created on: Feb 28, 2010
 *      Author: benk
 */

#ifndef FITOBLOOPOPERATOR_HPP_
#define FITOBLOOPOPERATOR_HPP_

#include "FitobOperatorBasis.hpp"
#include "src/expressions/FitobExpressionBasis.hpp"

namespace fitob{

  using namespace std;

  /** The operator to repeat one operation sequence */
  class LoopOperator: public fitob::OperatorBasis {
  public:

	LoopOperator(ExpressionBasis* loopCondition,
			  OperatorBasis* loopBody);

	virtual ~LoopOperator();

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

	string toString() { return operationName_ + " LoopCond:" +
		loopCondition_->toString() + " loopBody:" + loopBody_->toString();}

private:

	/** */
	int nrLoop_;

	/** expression in the loop which will be casted to int*/
	ExpressionBasis* loopCondition_;

	/** body of the loop */
	OperatorBasis* loopBody_;

  };
}

#endif /* FITOBLOOPOPERATOR_HPP_ */

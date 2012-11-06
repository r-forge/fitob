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
 * FitobIfOperator.hpp
 *
 *  Created on: Feb 28, 2010
 *      Author: benk
 */

#ifndef FITOBIFOPERATOR_HPP_
#define FITOBIFOPERATOR_HPP_

#include "FitobOperatorBasis.hpp"
#include "src/utils/fitobdefs.hpp"
#include "src/expressions/FitobExpressionBasis.hpp"

namespace fitob{

  using namespace std;

  /** The operator IF and IF ... ELSE <br>
   * The false branch is optional and the presence of that branch is marked by an internal flag <br>
   * IF operator can basically be of two different types, depending on the nature of the condition expression <br>
   * - The first case is when the condition expression is constant in this case the IF operator is just a branching
   * operator. In the forward estimation the branching decision can be stored and by the backward calculation recalled.
   * This is a rather simpler case <br> <br>
   * - The second case is when the condition expression is not constant (this means that in the expression there are
   * one or more import variables which form the axis of the actual mesh (this we can find out already in the forward
   * estimation phase) or the export variable is involved in the condition expression).
   * In this case a branching would cause the "division" of the mesh, this is rather complicated and not allowed
   * at the current version. <br>
   * Therefore for not constant expressions we only allow to have one MeshAssigmentOperator operator either branches of the
   * IF operator, so in this case the IF operator will be a condition MeshAssigmentOperator operator, which alone does not exist
   * If the body of the TRUE or FALSE branch is not only a single MeshAssigmentOperator */
  class IfOperator: public fitob::OperatorBasis {
  public:

	/** empty Ctor*/
	IfOperator();

	virtual ~IfOperator();

	/** initialize the object */
	void init(ExpressionBasis* condition,
			  OperatorBasis* trueBrach){
		hasElseBranch_ = false;
		condition_ = condition;
		trueBranch_ = trueBrach;
	};

	/** initialize the object */
	void init(ExpressionBasis* condition,
			  OperatorBasis* trueBrach,
			  OperatorBasis* falseBrach){
		hasElseBranch_ = true;
		condition_ = condition;
		trueBranch_ = trueBrach;
		falseBranch_ = falseBrach;
	};

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

	/** apply constrains on the grid */
	void applyConstraints(DVector& globalCoords) const;

	string toString() {
		if (hasElseBranch_){
			return operationName_ + " Exp:" + condition_->toString() +
			+ " TrueBranch:" +trueBranch_->toString() +
			  " FalseBranch:"+falseBranch_->toString();
		} else {
			return operationName_ + " Exp:" + condition_->toString() +
			+ " TrueBranch:" +trueBranch_->toString();
		}
    }

  private:

	/** TRUE if the else operator has false branch, FALSE otherwise */
	bool hasElseBranch_;

	/** In the forward estimation in case of constant expression in the condition
	 * we should store the value of the evaluation of the condition */
	bool constConditionValue_;

    /** TRUE if the condition expression is constant, FALSE otherwise <br>
     * the condition is not constant if the condition expression contains either
     * axis variables (import variables which are not constant) or the export variable <br>*/
	bool constCondition_;

	/** expression with the condition */
	ExpressionBasis* condition_;

	/** the operation(s) on the true branch */
	OperatorBasis* trueBranch_;

	/** the operation(s) on the false branch, if there exists a false branch */
	OperatorBasis* falseBranch_;
  };
}

#endif /* FITOBIFOPERATOR_HPP_ */

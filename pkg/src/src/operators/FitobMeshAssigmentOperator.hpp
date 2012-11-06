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
 * FitobMeshAssigmentOperator.hpp
 *
 *  Created on: Feb 28, 2010
 *      Author: benk
 */

#ifndef FITOBMESHASSIGMENTOPERATOR_HPP_
#define FITOBMESHASSIGMENTOPERATOR_HPP_

#include "FitobOperatorBasis.hpp"
#include "src/expressions/FitobExpressionBasis.hpp"
#include "src/variables/FitobVariable.hpp"

namespace fitob{

  using namespace std;

  /** This is the operation to assign values to the mesh (point by point)*/
  class MeshAssigmentOperator: public fitob::OperatorBasis {
  public:

	MeshAssigmentOperator(const Variable* expVariable, ExpressionBasis* assignExpression);

	virtual ~MeshAssigmentOperator();

	/** function to do forward estimation */
	virtual void forwardEstimation( boost::ptr_vector<MeshContext>& contextStack ,
			                        FitobCalculator* calc ,
			                        int& stackIndex ,
			                        double& timeStamp );

	/** function to do backward calculation */
	virtual void backwardCalculation( boost::ptr_vector<MeshContext>& contextStack ,
			                          FitobCalculator* calc ,
				                      int& stackIndex ,
			                          double& timeStamp );

	/** function to apply the constraint on the mesh */
	void applyConstraints(DVector& globalCoords) const;

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

	/** returns const pointer on the expression which will be used for the mesh assignment <br>
	 * this is needed for thr IF ELSE construct to build a special expression with condition and assignment expression */
	const ExpressionBasis* returnExpression() const { return assignExpression_; }

	/** returns the export variable which is the target of thie operations*/
	inline const Variable* getExportvariable() const { return expVariable_; }

	string toString() {
		return operationName_ + " ExpVar " +" Exp:" + assignExpression_->toString();}

  private:

	/** the export variable which should be assigned */
	const Variable* expVariable_;

	/** the assignment expression */
	ExpressionBasis* assignExpression_;

};
}

#endif /* FITOBMESHASSIGMENTOPERATOR_HPP_ */

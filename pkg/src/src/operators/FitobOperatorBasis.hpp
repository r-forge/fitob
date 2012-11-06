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
 * FitobOperatorBasis.hpp
 *
 *  Created on: Feb 25, 2010
 *      Author: benk
 */

#ifndef FITOBOPERATORBASIS_HPP_
#define FITOBOPERATORBASIS_HPP_

#include "src/utils/fitobdefs.hpp"
#include "src/evalcontext/FitobMeshContext.hpp"

namespace fitob{

  // forward declaration of these objects
  class OperatorSequence;
  class FitobCalculator;
  class MCMachine;


/** */
  typedef enum{ DaOp   = 0 ,
                IfOp   = 1 ,
                LoopOp = 2 ,
                AssigmentOp = 3 ,
                MeshAssigmentOp = 7 ,
                OperatorSeq = 4 ,
                ThOperator = 5 ,
                PlOperator = 8 ,
                NoOp = 6 } OperatorType;

  using namespace std;

  /** This is the superclass of any type of operator */
  class OperatorBasis : boost::noncopyable , public fitob::VerbClass {
  public:

	/** Ctor */
	OperatorBasis();

	/** Ctor */
	OperatorBasis(const string& name, const OperatorType& opType);

	virtual ~OperatorBasis();

	/** return the name of the operator */
	string& getOpName() { return operationName_;}

	/** return the name of the operator */
	virtual string toString() { return operationName_;}

	/** function to do forward estimation */
	virtual void forwardEstimation( boost::ptr_vector<MeshContext>& contextStack ,
			                        FitobCalculator* calc ,
			                        int& stackIndex ,
			                        double& timeStamp ){
		// throw an exception if this is not overwritten
		FITOB_ERROR_EXIT("OperatorBasis::forwardEstimation , must be overwritten! ");
     }

	/** function for diffusion axis size (enlargement) estimation */
	virtual void forwardEstimation_DiffusionEnlarement(
			boost::ptr_vector<MeshContext>& contextStack ,
			                          FitobCalculator* calc ,
			                          DVector& gradValues,
			                          int globalVariableIndex,
				                      int& stackIndex ,
			                          double timeStamp ) const {
		// this method should be implemented only by the Da operator, for the case
		// when the Da operator shifts one diffusion axis
	}

	/** function to do backward calculation */
	virtual void backwardCalculation( boost::ptr_vector<MeshContext>& contextStack ,
			                          FitobCalculator* calc ,
				                      int& stackIndex ,
			                          double& timeStamp ){
		// throw an exception if this is not overwritten
		FITOB_ERROR_EXIT("OperatorBasis::backwardCalculation , must be overwritten! ");
	}

	/** function to apply the constraints */
	virtual void applyConstraints(DVector& globalCoords) const{
        // this will not be overwritten only by some operators, if not and is called then exception will be thrown
		FITOB_ERROR_EXIT("OperatorBasis::applyConstraints , must be overwritten! , PROBABLY UNALLOWED OPERATOR IN CONSTRAIN !!!");
	}

	/** function for simulate the scenarios forward */
	virtual void forwardMCSimulation( boost::ptr_vector<MeshContext>& contextStack ,
									  FitobCalculator* calc , MCMachine* MC ,
									  int& stackIndex , double& timeStamp ) {
		// throw an exception if this is not overwritten
		FITOB_ERROR_EXIT("OperatorBasis::forwardMCSimulation , must be overwritten! ");
     }

	/** function for evaluate the scenarios forward */
	virtual void forwardMCEvaluation( boost::ptr_vector<MeshContext>& contextStack ,
									  FitobCalculator* calc , MCMachine* MC ,
									  int& stackIndex , double& timeStamp ) {
		// throw an exception if this is not overwritten
		FITOB_ERROR_EXIT("OperatorBasis::forwardMCEvaluation , must be overwritten! ");
     }

	/** function for evaluate the scenarios backward */
	virtual void backwardMCEvaluation( boost::ptr_vector<MeshContext>& contextStack ,
									  FitobCalculator* calc , MCMachine* MC ,
									  int& stackIndex , double& timeStamp ) {
		// throw an exception if this is not overwritten
		FITOB_ERROR_EXIT("OperatorBasis::backwardMCEvaluation , must be overwritten! ");
     }

	/** get the time of operator execution start */
	double operatorTime() const {return operatorTime_;}

	/** set the time of the start of operator execution */
	void setOperatorTime(double v) { operatorTime_ = v; }

	/** returns type of the operator */
	const OperatorType& operatorType() const { return operatorType_; }

  protected:

	/** name of the operation */
	string operationName_;

	/** type of the operator*/
	OperatorType operatorType_;

	/** time stamp when the execution of this operator starts */
	double operatorTime_;

 };
}
#endif /* FITOBOPERATORBASIS_HPP_ */

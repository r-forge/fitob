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
 * FitobOperatorSequence.hpp
 *
 *  Created on: Feb 28, 2010
 *      Author: benk
 */

#ifndef FITOBOPERATORSEQUENCE_HPP_
#define FITOBOPERATORSEQUENCE_HPP_

#include "FitobOperatorBasis.hpp"
#include "src/utils/fitobdefs.hpp"

namespace fitob{

  using namespace std;

  /** This class is one bunch of operators */
  class OperatorSequence: public fitob::OperatorBasis {

  public:

	/** Ctor */
	OperatorSequence();

	virtual ~OperatorSequence();

	/** add one operator to the sequence*/
	void addOperatorToSequence(OperatorBasis* op){
		operatorSequqnce_.resize(nrOperations_+1);
		operatorSequqnce_[nrOperations_] = op;
		nrOperations_++;
	}

	/** return the operator */
	OperatorBasis* returnOperator(const int i){
		return operatorSequqnce_[i];
	}

	/** return the number of operator */
	int nrOperators(){
		return nrOperations_;
	}

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

	/** function to impose the constraints */
	void applyConstraints(DVector& globalCoords) const;

	string toString() {
		string outstr = "\nOPSEQ:";
		for (int i = 0 ; i < nrOperations_ ; i++){
			//std::cout << "call toString()";
			//std::cout << operatorSequqnce_[i]->toString() <<std::endl;
			outstr = outstr + operatorSequqnce_[i]->toString() + "\n";
		}
		return outstr + "END-OPSEQ:\n";}

  private:

	/** the operator stored in array */
	std::vector<OperatorBasis*> operatorSequqnce_;

	/** Number of operators */
	int nrOperations_;

  };
}

#endif /* FITOBOPERATORSEQUENCE_HPP_ */

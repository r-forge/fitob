/*
 * FitobOperatorSequence.cpp
 *
 *  Created on: Feb 28, 2010
 *      Author: benk
 */

#include "FitobOperatorSequence.hpp"
#include "FitobOperatorBasis.hpp"

using namespace std;
using namespace fitob;

OperatorSequence::OperatorSequence() :
 OperatorBasis("OperatorSequence",  (OperatorSeq) )
 , operatorSequqnce_() , nrOperations_(0)
{
	setVerb(0);
}

OperatorSequence::~OperatorSequence() {
}


void OperatorSequence::forwardEstimation(
		boost::ptr_vector<MeshContext>& contextStack ,
		                        FitobCalculator* calc ,
		                        int& stackIndex ,
		                        double& timeStamp ){

  // here we just call each operator in forward direction
  for (int ind = 0 ; ind < nrOperations_ ; ind++){
	  operatorSequqnce_[ind]->forwardEstimation( contextStack , calc , stackIndex , timeStamp );
	  FITOB_OUT_LEVEL3(verb(),"OperatorSequence::forwardEstimation:" << stackIndex <<
			  " contextStack.size():" << contextStack.size() << " ind:"<<ind);
  }
}

void OperatorSequence::forwardEstimation_DiffusionEnlarement(
		                          boost::ptr_vector<MeshContext>& contextStack ,
		                          FitobCalculator* calc ,
		                          DVector& gradValues,
		                          int globalVariableIndex,
			                      int& stackIndex ,
		                          double timeStamp ) const {
	// if the index already below the index then nothing to do
	// this is the exit condition for instructions below "timeStamp"
	if ( stackIndex >= (int)contextStack.size() ) return; // EXIT IF NEEDED

    // here do traversal in forward direction
	  for (int ind = 0 ; ind < nrOperations_ ; ind++){
		  operatorSequqnce_[ind]->forwardEstimation_DiffusionEnlarement( contextStack , calc , gradValues ,
			               globalVariableIndex , stackIndex , timeStamp);
		  FITOB_OUT_LEVEL3(verb(),"OperatorSequence::forwardEstimation_DiffusionEnlarement:" << stackIndex
				  << " contextStack.size():" << contextStack.size() << " ind:"<<ind);
	  }
}

void OperatorSequence::backwardCalculation(
		boost::ptr_vector<MeshContext>& contextStack ,
		                          FitobCalculator* calc ,
			                      int& stackIndex ,
		                          double& timeStamp ){
  // Operators must be traversed in reverse order !!!!
  // RESPECT REVERSE ORDER !!!
	  for (int ind = nrOperations_ - 1 ; ind >= 0 ; ind-- ){
		  operatorSequqnce_[ind]->backwardCalculation( contextStack , calc , stackIndex , timeStamp );
	  }
}


void OperatorSequence::applyConstraints(DVector& globalCoords) const{
	// just call the constrain methods for all operations in the sequence FORWARD ORDER
	  for (int ind = 0 ; ind < nrOperations_ ; ind++){
		  operatorSequqnce_[ind]->applyConstraints(globalCoords);
	  }
}


void OperatorSequence::forwardMCSimulation( boost::ptr_vector<MeshContext>& contextStack ,
								  FitobCalculator* calc , MCMachine* MC ,
								  int& stackIndex , double& timeStamp ) {
	  // just call the constrain methods for all operations in the sequence FORWARD ORDER
	  for (int ind = 0 ; ind < nrOperations_ ; ind++){
		  operatorSequqnce_[ind]->forwardMCSimulation(  contextStack , calc , MC , stackIndex , timeStamp );
	  }
}


void OperatorSequence::forwardMCEvaluation( boost::ptr_vector<MeshContext>& contextStack ,
								  FitobCalculator* calc , MCMachine* MC ,
								  int& stackIndex , double& timeStamp ) {
	  // just call the constrain methods for all operations in the sequence FORWARD ORDER
	  for (int ind = 0 ; ind < nrOperations_ ; ind++){
		  operatorSequqnce_[ind]->forwardMCEvaluation(  contextStack , calc , MC , stackIndex , timeStamp );
	  }
}


void OperatorSequence::backwardMCEvaluation( boost::ptr_vector<MeshContext>& contextStack ,
								  FitobCalculator* calc , MCMachine* MC ,
								  int& stackIndex , double& timeStamp ) {
	  // Operators must be traversed in reverse order !!!!
	  // RESPECT REVERSE ORDER !!!
	  for (int ind = nrOperations_ - 1 ; ind >= 0 ; ind-- ){
		  operatorSequqnce_[ind]->backwardMCEvaluation( contextStack , calc , MC , stackIndex , timeStamp );
	  }
}

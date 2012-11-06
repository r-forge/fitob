/*
 * FitobIfOperator.cpp
 *
 *  Created on: Feb 28, 2010
 *      Author: benk
 */

#include "FitobIfOperator.hpp"
#include "src/operators/FitobMeshAssigmentOperator.hpp"
#include "src/expressions/FitobIfElseAssigmentExpression.hpp"
#include "src/scripteval/FitobCalculator.hpp"
#include "src/montecarlo/FitobMCMachine.hpp"

using namespace std;
using namespace fitob;

IfOperator::IfOperator() : OperatorBasis("IfOperator" , (IfOp) )
 , constConditionValue_(true) , condition_() , trueBranch_() , falseBranch_() {

	setVerb(6);
}

IfOperator::~IfOperator() {
}

void IfOperator::forwardEstimation(
		boost::ptr_vector<MeshContext>& contextStack ,
		                        FitobCalculator* calc ,
		                        int& stackIndex ,
		                        double& timeStamp ){

	// see if the expression in the condition constant is
	// we use the context on the top of the stack
	MeshContext& actualContext = contextStack[stackIndex];

	constCondition_ = condition_->isConstantExpression(actualContext);
	//FITOB_OUT_LEVEL3( verb() , " IfOperator::applyConstraints CONST: " << constCondition_ << " actualContext:" << actualContext.getDom().toString());

	if (constCondition_)
	{  // expression is constant so we have a simple IF branch operator
		constConditionValue_ = (condition_->eval(actualContext.minGlobCoord()) > FITOB_NUMERICALZERO ) ? true : false ;
	}
	else
	{ // expression is not constant we most have an Mesh assigment operator in both bodies
		// first test if the IF operator has else branch
		if (hasElseBranch_)
		{
		   MeshAssigmentOperator* trueBranch = dynamic_cast<MeshAssigmentOperator*>(trueBranch_);
		   FITOB_ERROR_TEST(trueBranch != 0," could not caste TRUE branch into MeshAssigmentOperator");

		   MeshAssigmentOperator* falseBranch = dynamic_cast<MeshAssigmentOperator*>(falseBranch_);
		   FITOB_ERROR_TEST(falseBranch != 0 ," could not caste FALSE branch into MeshAssigmentOperator");

		   // nothing to do here more
		}
		else // IF operator has no false (else) branch
		{
		   MeshAssigmentOperator* trueBranch = dynamic_cast<MeshAssigmentOperator*>(trueBranch_);
		   FITOB_ERROR_TEST(trueBranch == 0 ," could not caste TRUE branch into MeshAssigmentOperator");

		   // nothing to do here more
		}
	}
}

/** function for diffusion axis size (enlargement) estimation */
void IfOperator::forwardEstimation_DiffusionEnlarement(
		boost::ptr_vector<MeshContext>& contextStack ,
		                          FitobCalculator* calc ,
		                          DVector& gradValues,
		                          int globalVariableIndex,
			                      int& stackIndex ,
		                          double timeStamp ) const {

	// if the index already below the index then nothing to do
	// this is the exit condition for instructions below "timeStamp"
	if ( stackIndex >= (int)contextStack.size() ) return; // EXIT IF NEEDED

	if (constCondition_)
	{  // expression is constant so we have a simple IF branch operator
		if (hasElseBranch_)
		{
		   if (constConditionValue_){ // in case of const condition we should have evaluated the condition
			   trueBranch_->forwardEstimation_DiffusionEnlarement(
					 contextStack , calc , gradValues , globalVariableIndex , stackIndex , timeStamp );
		   }else{
			   falseBranch_->forwardEstimation_DiffusionEnlarement(
					 contextStack , calc , gradValues , globalVariableIndex , stackIndex , timeStamp );
		   }
		}
		else
		   if (constConditionValue_)
		   {
			   trueBranch_->forwardEstimation_DiffusionEnlarement(
					 contextStack , calc , gradValues , globalVariableIndex , stackIndex , timeStamp );
		   }
	}
	// in other cases we do not have to do anything
}

void IfOperator::backwardCalculation(
		boost::ptr_vector<MeshContext>& contextStack ,
		                          FitobCalculator* calc ,
			                      int& stackIndex ,
		                          double& timeStamp ){
	if (constCondition_)
	{  // expression is constant so we have a simple IF branch operator
	   // according to the precalculated constant condition calculation we take the branch
		if (hasElseBranch_)
		{
			 if (constConditionValue_){
				 trueBranch_->backwardCalculation( contextStack , calc , stackIndex , timeStamp);
			 } else {
				 falseBranch_->backwardCalculation( contextStack , calc , stackIndex , timeStamp);
			 }
		}
		else
		{
			 if (constConditionValue_){
				 trueBranch_->backwardCalculation( contextStack , calc , stackIndex , timeStamp);
			 }
		}
	}
	else
	{// expression is not constant we most have an Mesh assigment operator in both bodies
		// first test if the IF operator has else branch
		if (hasElseBranch_)
		{
		   MeshAssigmentOperator* trueBranch = dynamic_cast<MeshAssigmentOperator*>(trueBranch_);
		   FITOB_ERROR_TEST(trueBranch != 0 , " could not caste TRUE branch into MeshAssigmentOperator");

		   MeshAssigmentOperator* falseBranch = dynamic_cast<MeshAssigmentOperator*>(falseBranch_);
		   FITOB_ERROR_TEST(falseBranch != 0 , " could not caste FALSE branch into MeshAssigmentOperator");

		   // ---- make mesh assignment ---
		   IfElseAssigmentExpression ifElseExp(condition_,trueBranch->returnExpression() , falseBranch->returnExpression() );
		   // if we get here then we just have a simple mesh assignment , no if condition will be considered here
		   contextStack[stackIndex].getMesh()->setValues(
				   &(ifElseExp) , contextStack[stackIndex].minGlobCoord());
		   // do eventually plotting
		   calc->getPlotter()->plot( &(contextStack[stackIndex]) , calc->getScriptModel()->getModelName() );
		}
		else // IF operator has no false (else) branch
		{
		   MeshAssigmentOperator* trueBranch = dynamic_cast<MeshAssigmentOperator*>(trueBranch_);
		   FITOB_ERROR_TEST(trueBranch != 0 , " could not caste TRUE branch into MeshAssigmentOperator");

		   // ---- make mesh assignment ---
		   IfElseAssigmentExpression ifExp(condition_,trueBranch->returnExpression() , NULL );
		   // if we get here then we just have a simple mesh assignment , no if condition will be considered here
		   contextStack[stackIndex].getMesh()->setValues(
				   &(ifExp) , contextStack[stackIndex].minGlobCoord());
		   // do eventually plotting
		   calc->getPlotter()->plot( &(contextStack[stackIndex]) , calc->getScriptModel()->getModelName() );
		}
	}
}


void IfOperator::applyConstraints(DVector& globalCoords) const {
/*	// code for constant expression
    if (constCondition_)
	{  // expression is constant so we have a simple IF branch operator
	   // according to the precalculated constant condition calculation we take the branch
		FITOB_OUT_LEVEL3( verb() , " IfOperator::applyConstraints CONST ");
		if (hasElseBranch_){
			 if (constConditionValue_){  trueBranch_->applyConstraints( globalCoords );}
			 else {  falseBranch_->applyConstraints( globalCoords ); }
		} else{
			 if (constConditionValue_){ trueBranch_->applyConstraints( globalCoords );}
		}
	}
	else */

    // ---- Now we assume that in the constraint block all IF operators are not constant ----
	// expression is not constant we most have an Mesh assigment operator in both bodies

        // evaluate the expression to decide which way to go
		double expression_val = condition_->eval(globalCoords);
        //FITOB_OUT_LEVEL3( verb() , " IfOperator::applyConstraints NOT CONST " << expression_val);
		if (hasElseBranch_)
		{
			 if (expression_val > FITOB_NUMERICALZERO){
				 trueBranch_->applyConstraints( globalCoords );
			 } else {
				 falseBranch_->applyConstraints( globalCoords );
			 }
		}
		else // IF operator has no false (else) branch
		{
			 if (expression_val > FITOB_NUMERICALZERO){
				 trueBranch_->applyConstraints( globalCoords );
			 }
		}

}



void IfOperator::forwardMCSimulation( boost::ptr_vector<MeshContext>& contextStack ,
								  FitobCalculator* calc , MCMachine* MC ,
								  int& stackIndex , double& timeStamp ) {

	if (constCondition_)
	{  // expression is constant so we have a simple IF branch operator
	   // according to the precalculated constant condition calculation we take the branch
		if (hasElseBranch_)
		{
			 if (constConditionValue_){
				 trueBranch_->forwardMCSimulation( contextStack , calc , MC , stackIndex , timeStamp);
			 } else {
				 falseBranch_->forwardMCSimulation( contextStack , calc , MC , stackIndex , timeStamp);
			 }
		}
		else
		{
			 if (constConditionValue_){
				 trueBranch_->forwardMCSimulation( contextStack , calc , MC , stackIndex , timeStamp);
			 }
		}
	}
	else
	{// expression is not constant we most have an Mesh assigment operator in both bodies
		// first test if the IF operator has else branch
		if (hasElseBranch_)
		{
		   MeshAssigmentOperator* trueBranch = dynamic_cast<MeshAssigmentOperator*>(trueBranch_);
		   FITOB_ERROR_TEST(trueBranch != 0," could not caste TRUE branch into MeshAssigmentOperator");

		   MeshAssigmentOperator* falseBranch = dynamic_cast<MeshAssigmentOperator*>(falseBranch_);
		   FITOB_ERROR_TEST(falseBranch != 0 , " could not caste FALSE branch into MeshAssigmentOperator");

		   // ---- make mesh assignment ---
		   IfElseAssigmentExpression ifElseExp(condition_,trueBranch->returnExpression() , falseBranch->returnExpression() );

		   FITOB_ERROR_TEST( trueBranch->getExportvariable()->getGlobalIndex() == falseBranch->getExportvariable()->getGlobalIndex() ,
				   " IfOperator::forwardMCSimulation , export variables must be the same : " << trueBranch->getExportvariable()->toString()
				   << " and : " << falseBranch->getExportvariable()->toString() );

		   // set the export values in the MC step
		   (MC->getMCStep(MC->nrMCSteps()-1)).applyExpression( &(MC->getMCStep(MC->nrMCSteps()-1)) ,
				   trueBranch->getExportvariable() , &(ifElseExp) , calc );
		}
		else // IF operator has no false (else) branch
		{
		   MeshAssigmentOperator* trueBranch = dynamic_cast<MeshAssigmentOperator*>(trueBranch_);
		   if (trueBranch == 0){
			   FITOB_ERROR_EXIT(" could not caste TRUE branch into MeshAssigmentOperator");
		   }

		   // ---- make mesh assignment ---
		   IfElseAssigmentExpression ifExp(condition_,trueBranch->returnExpression() , NULL );

		   // set the export values in the MC step
		   (MC->getMCStep(MC->nrMCSteps()-1)).applyExpression( &(MC->getMCStep(MC->nrMCSteps()-1)) ,
				   trueBranch->getExportvariable() , &(ifExp) , calc);
		}
	}
}


void IfOperator::forwardMCEvaluation( boost::ptr_vector<MeshContext>& contextStack ,
								  FitobCalculator* calc , MCMachine* MC ,
								  int& stackIndex , double& timeStamp ) {
	if (constCondition_)
	{  // expression is constant so we have a simple IF branch operator
	   // according to the precalculated constant condition calculation we take the branch
		if (hasElseBranch_)
		{
			 if (constConditionValue_){
				 trueBranch_->forwardMCEvaluation( contextStack , calc , MC , stackIndex , timeStamp);
			 } else {
				 falseBranch_->forwardMCEvaluation( contextStack , calc , MC , stackIndex , timeStamp);
			 }
		}
		else
		{
			 if (constConditionValue_){
				 trueBranch_->forwardMCEvaluation( contextStack , calc , MC , stackIndex , timeStamp);
			 }
		}
	}
	else
	{// expression is not constant we most have an Mesh assigment operator in both bodies
		// first test if the IF operator has else branch
		if (hasElseBranch_)
		{
		   MeshAssigmentOperator* trueBranch = dynamic_cast<MeshAssigmentOperator*>(trueBranch_);
		   FITOB_ERROR_TEST( trueBranch != 0 , " could not caste TRUE branch into MeshAssigmentOperator");

		   MeshAssigmentOperator* falseBranch = dynamic_cast<MeshAssigmentOperator*>(falseBranch_);
		   FITOB_ERROR_TEST(falseBranch != 0 , " could not caste FALSE branch into MeshAssigmentOperator");

		   // ---- make mesh assignment ---
		   IfElseAssigmentExpression ifElseExp(condition_,trueBranch->returnExpression() , falseBranch->returnExpression() );

		   FITOB_ERROR_TEST( trueBranch->getExportvariable()->getGlobalIndex() == falseBranch->getExportvariable()->getGlobalIndex() ,
				   " IfOperator::forwardMCSimulation , export variables must be the same : " << trueBranch->getExportvariable()->toString()
				   << " and : " << falseBranch->getExportvariable()->toString() );

		   // set the export values in the MC step
		   (MC->getMCStep(MC->nrMCSteps()-1)).applyExpression( &(MC->getMCStep(MC->nrMCSteps()-1)) ,
				   trueBranch->getExportvariable() , &(ifElseExp) , calc );
		}
		else // IF operator has no false (else) branch
		{
		   MeshAssigmentOperator* trueBranch = dynamic_cast<MeshAssigmentOperator*>(trueBranch_);
		   FITOB_ERROR_TEST(trueBranch != 0 , " could not caste TRUE branch into MeshAssigmentOperator");

		   // ---- make mesh assignment ---
		   IfElseAssigmentExpression ifExp(condition_,trueBranch->returnExpression() , NULL );

		   // set the export values in the MC step
		   (MC->getMCStep(MC->nrMCSteps()-1)).applyExpression( &(MC->getMCStep(MC->nrMCSteps()-1)) ,
				   trueBranch->getExportvariable() , &(ifExp) , calc );
		}
	}
}


void IfOperator::backwardMCEvaluation( boost::ptr_vector<MeshContext>& contextStack ,
								  FitobCalculator* calc , MCMachine* MC ,
								  int& stackIndex , double& timeStamp ) {
	if (constCondition_)
	{  // expression is constant so we have a simple IF branch operator
	   // according to the precalculated constant condition calculation we take the branch
		if (hasElseBranch_)
		{
			 if (constConditionValue_){
				 trueBranch_->backwardMCEvaluation( contextStack , calc , MC , stackIndex , timeStamp);
			 } else {
				 falseBranch_->backwardMCEvaluation( contextStack , calc , MC , stackIndex , timeStamp);
			 }
		}
		else
		{
			 if (constConditionValue_){
				 trueBranch_->backwardMCEvaluation( contextStack , calc , MC , stackIndex , timeStamp);
			 }
		}
	}
	else
	{// expression is not constant we most have an Mesh assigment operator in both bodies
		// first test if the IF operator has else branch
		if (hasElseBranch_)
		{
		   MeshAssigmentOperator* trueBranch = dynamic_cast<MeshAssigmentOperator*>(trueBranch_);
		   FITOB_ERROR_TEST( trueBranch != 0 , " could not caste TRUE branch into MeshAssigmentOperator");

		   MeshAssigmentOperator* falseBranch = dynamic_cast<MeshAssigmentOperator*>(falseBranch_);
		   FITOB_ERROR_TEST(falseBranch != 0 , " could not caste FALSE branch into MeshAssigmentOperator");

		   // ---- make mesh assignment ---
		   IfElseAssigmentExpression ifElseExp(condition_,trueBranch->returnExpression() , falseBranch->returnExpression() );

		   FITOB_ERROR_TEST( trueBranch->getExportvariable()->getGlobalIndex() == falseBranch->getExportvariable()->getGlobalIndex() ,
				   " IfOperator::forwardMCSimulation , export variables must be the same : " << trueBranch->getExportvariable()->toString()
				   << " and : " << falseBranch->getExportvariable()->toString() );

		   // set the export values in the MC step
		   (MC->getMCStep(MC->nrMCSteps()-1)).applyExpression( &(MC->getMCStep(MC->nrMCSteps()-1)) ,
				   trueBranch->getExportvariable() , &(ifElseExp) , calc );
		}
		else // IF operator has no false (else) branch
		{
		   MeshAssigmentOperator* trueBranch = dynamic_cast<MeshAssigmentOperator*>(trueBranch_);
		   FITOB_ERROR_TEST(trueBranch != 0 , " could not caste TRUE branch into MeshAssigmentOperator");

		   // ---- make mesh assignment ---
		   IfElseAssigmentExpression ifExp(condition_,trueBranch->returnExpression() , NULL );

		   // set the export values in the MC step
		   (MC->getMCStep(MC->nrMCSteps()-1)).applyExpression( &(MC->getMCStep(MC->nrMCSteps()-1)) ,
				   trueBranch->getExportvariable() , &(ifExp) , calc );
		}
	}
}

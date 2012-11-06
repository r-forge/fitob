/*
 * FitobDaOperator.cpp
 *
 *  Created on: Feb 28, 2010
 *      Author: benk
 */

#include "FitobDaOperator.hpp"
#include "src/evalcontext/FitobDaOperatorEvaluable.hpp"
#include "src/scripteval/FitobCalculator.hpp"
#include "src/montecarlo/FitobMCMachine.hpp"

using namespace std;
using namespace fitob;

DaOperator::DaOperator( const boost::shared_ptr<XMLConfiguration>& xmlconfiguration ,
		Variable* assignedVariable, ExpressionBasis* expression) :
 OperatorBasis("DaOperator" , (DaOp) ) , assignedVariable_(assignedVariable)
 , expression_(expression) , expressionIsConstant_(false) , alwaysIncreasingDa_(false) {
	setVerb(0);

	if (xmlconfiguration->getStringConfiguration("thetaconfigurations.daoperator.INCREASING_DA.<xmlattr>.value") == "true"){
	  alwaysIncreasingDa_ = true;
	}
}

DaOperator::~DaOperator() {
}


void DaOperator::forwardEstimation(
		boost::ptr_vector<MeshContext>& contextStack ,
		                        FitobCalculator* calc ,
		                        int& stackIndex ,
		                        double& timeStamp ){

	// we can assume that the shifted variable is important
	FITOB_OUT_LEVEL3(verb(),"DaOperator::forwardEstimation stackIndex:" << stackIndex << " contextStack.size():" << contextStack.size()
			<< "  globalVariableIndex:" << assignedVariable_->getGlobalIndex() );
	MeshContext& actualContext = contextStack[stackIndex];

	// see if we have a constant expression
	expressionIsConstant_ = expression_->isConstantExpression(actualContext);
	FITOB_OUT_LEVEL1(verb(),"DaOperator::forwardEstimation expression_:" << expression_->toString() );

	// create the cloned domain
	const Domain& actDom = actualContext.getDom();
	Domain newDomain(actualContext.getDom());

	DVector shiftAxis = actDom.getGradedAxis(assignedVariable_->getGlobalIndex());
	//do the shift, of the shiftAxis, using actDom
	doShift_forward(shiftAxis , actDom , timeStamp);

	// sort shiftAxis vector
	std::sort( shiftAxis.begin(),shiftAxis.end() );

	// update the domain
	newDomain.setAxis( shiftAxis , assignedVariable_->getGlobalIndex() );

	// add the new context to the stack
	//contextStack.resize(contextStack.size()+1);
	contextStack.push_back(new MeshContext(newDomain,timeStamp));
	FITOB_OUT_LEVEL3(verb(),"DaOperator::forwardEstimation stackIndex = " << contextStack.size()-1);
	FITOB_OUT_LEVEL3(verb(),"DaOperator::forwardEstimation New Domain:" << newDomain.toString() );
	stackIndex = contextStack.size()-1;
}

void DaOperator::forwardEstimation_DiffusionEnlarement(
		boost::ptr_vector<MeshContext>& contextStack ,
		                          FitobCalculator* calc ,
		                          DVector& gradValues,
		                          int globalVariableIndex,
			                      int& stackIndex ,
		                          double timeStamp ) const {

	// if the index already below the index then nothing to do
	// this is the exit condition for instructions below "timeStamp"
	FITOB_OUT_LEVEL3(verb(),"DaOperator::forwardEstimation_DiffusionEnlarement stackIndex:" << stackIndex
			<< " contextStack.size():" << contextStack.size() );
	if ( stackIndex >= (int)contextStack.size() ) return; // EXIT IF NEEDED

	// this method should be implemented only by the Da operator, for the case
	// when the Da operator shifts one diffusion axis

	FITOB_OUT_LEVEL3(verb(),"DaOperator::forwardEstimation_DiffusionEnlarement: test var:" << assignedVariable_->getGlobalIndex()
			<< "  globalVariableIndex:" << globalVariableIndex );

	if (assignedVariable_->getGlobalIndex() == globalVariableIndex){
        // do the shift using stackIndex
		FITOB_OUT_LEVEL3(verb(),"DaOperator::forwardEstimation_DiffusionEnlarement: ");
		// we must clone the actual domain and replace the impacted axis, so that the actual values will be used
		// this is important for the cases e.g. S = S*0.9;
		Domain tmpDom = contextStack[stackIndex].getDom();
		tmpDom.setAxis( gradValues , globalVariableIndex );
		doShift_forward( gradValues , tmpDom , contextStack[stackIndex].getTime() );
	}
	// increase the stack index
	stackIndex = stackIndex + 1;
}


void DaOperator::doShift_forward(DVector& axisVector , const Domain& actDom , double actTime ) const
{
   // --------------------- DA OPERATOR FORWARD ESTIMATION -------------

	int nrImportVariable = actDom.nrImportVariables();

	DVector tmpEval(actDom.nrGlobalVariables());
	tmpEval[0] = actTime;
	// it is important that if the expression has export variable than this will be set to ZERO
	tmpEval[1] = 0.0;
	for (int jj = 0 ; jj < nrImportVariable ; jj++ )
		tmpEval[jj+2] = actDom.getGradedAxisMin(jj+2);

	if ( (axisVector.size() > 1) && (expressionIsConstant_ == true) )
	{   // in this case one not constant axis should be constant ... we do not allow this
		FITOB_ERROR_EXIT("DaOperator::doShift_forward not constant axis can not become constant! : " << axisVector.size());
	}

	if (expressionIsConstant_)
	{
		if (axisVector.size() != 1) {  axisVector.resize(1);  }
		// set the only value
		axisVector[0] = expression_->eval(tmpEval);
	}
	else
	{
		// determine the order of the shift
		IVector startVect(nrImportVariable, 0 );
		IVector increment(nrImportVariable, 1 );

		if (alwaysIncreasingDa_){
			// here we only set the increments to zero if the axis is constant
		    for (int aI = 0 ; aI < nrImportVariable ; aI++)
		    {
			   FITOB_OUT_LEVEL3(verb(),"DaOperator (increasing DA) maxPoints: " << fabs(actDom.getGradedAxisMax(aI+2) - actDom.getGradedAxisMin(aI+2)));
			   if (fabs(actDom.getGradedAxisMax(aI+2) - actDom.getGradedAxisMin(aI+2)) < FITOB_NUMERICALZERO )
			   {
				  startVect[aI] = 0; increment[aI] = 0; // const variable
			   }
		    }
		}
		else
		{
		    double minDaA, maxDaA; // these will be the final intervals for the shifted axis
		    DVector evalCoord(nrImportVariable+2); // coordinates to evaluate expression
		    IVector minP(nrImportVariable);
		    IVector maxP(nrImportVariable);
		    // ----
		    FITOB_OUT_LEVEL3(verb(),"DaOperator::doShift_forward , calling calcForward_rec , nrImportVariable:" << nrImportVariable);
		    FITOB_OUT_LEVEL3(verb(),"DaOperator::doShift_forward domain:" << actDom.toString() );
		    evalCoord[0] = actTime; // time
		    evalCoord[1] = 0.0; // export variable
		    calcForward_rec( nrImportVariable - 1 , actDom , minDaA , maxDaA , evalCoord, minP, maxP);

		    // minP contains the min Points combination , maxP is has the combinations (MIN,MAX,..) for the max
		    for (int aI = 0 ; aI < nrImportVariable ; aI++)
		    {
			   FITOB_OUT_LEVEL3(verb(),"DaOperator maxPoints: " << minP[aI] << "," << maxP[aI] << " ="
					<< fabs(actDom.getGradedAxisMax(aI+2) - actDom.getGradedAxisMin(aI+2)));
			   if (fabs(actDom.getGradedAxisMax(aI+2) - actDom.getGradedAxisMin(aI+2)) < FITOB_NUMERICALZERO )
			   {
				  startVect[aI] = 0; increment[aI] = 0; // const variable
			   }
			   else
			   {  // this is a real axis ( non constant import variable )
			      increment[aI] = maxP[aI] - minP[aI];
			      startVect[aI] = (maxP[aI] < minP[aI]) ? actDom.getGradedAxis(aI+2).size()-1 : 0 ;
			   }
			   FITOB_OUT_LEVEL3(verb(),"DaOperator inc. start. : " << startVect[aI] << "," << increment[aI] );
		    }
            // -------- end of determining directions
		}

		// we make the forward estimation with maximal resolution, because it is cheap
		int nrGradPerAxis = actDom.getMaxGradAxisSize();
		axisVector.resize(nrGradPerAxis);
		//FITOB_ERROR_TEST((bool)(nrGradPerAxis == (int)axisVector.size()),
		//		"DaOperator::doShift_forward axis size do not match nrGradPerAxis:" << nrGradPerAxis << " axisVector.size():"<<axisVector.size())

		// Shift the graded axis as we stored the increments and the start vectors
		for (int gI = 0 ; gI < nrGradPerAxis ; gI++){
			for (int aI = 0 ; aI < nrImportVariable ; aI++ ){
				tmpEval[aI+2] = actDom.getGradedAxis(aI+2)[startVect[aI]];
				FITOB_OUT_LEVEL3(verb(),"aI:" << aI << " V:" << tmpEval[aI+2] <<
						" startVect[aI]" << startVect[aI] << " increment[aI]:"<<increment[aI]);
			}
			// eval the Da operator
			axisVector[gI] = expression_->eval(tmpEval);
			FITOB_OUT_LEVEL3(verb(),"DaOperator::doShift_forward res:" << axisVector[gI] << " Expression:" << expression_->toString());
            // increment the index vectors with the calculated increments
			for (int aI = 0 ; aI < nrImportVariable ; aI++ ){
				// todo: increment this only if the axis is not constant
				startVect[aI] += increment[aI]; // we store into the startVect the actual position
			}
		}
	}
}


/** recursive function to determine the maximal size */
void DaOperator::calcForward_rec(int index, const Domain& actDom ,
		double& minDaA ,
		double& maxDaA ,
		DVector& evalCoord,
		IVector& minP,
		IVector& maxP) const {

	FITOB_OUT_LEVEL3(verb(), "DaOperator::calcForward_rec, index:" << index );
	if (index < 0){
		maxDaA = expression_->eval(evalCoord);
		minDaA = maxDaA;
		FITOB_OUT_LEVEL3(verb(), "DaOperator::calcForward_rec, res:" << minDaA << " Expression:" << expression_->toString());
		if (verb() > 4){
			for (unsigned int mm = 0 ; mm < evalCoord.size() ; mm++)
				std::cout << evalCoord[mm] << ",";
			std::cout << std::endl;
		}
		return;
	}

	// clone these arrays
	DVector evalNextLevel = evalCoord;
	IVector minPNextLevel1 = minP;
	IVector maxPNextLevel1 = maxP;
	IVector minPNextLevel2 = minP;
	IVector maxPNextLevel2 = maxP;

	double minVal1, maxVal1;
	double minVal2, maxVal2;
	// call recursively for min
	evalNextLevel[index+2] = actDom.getGradedAxisMin(index+2);
	FITOB_OUT_LEVEL3(verb(), "DaOperator::calcForward_rec, call REC MIN:" << (index-1) );
	calcForward_rec( index - 1 , actDom , minVal1, maxVal1 ,
			evalNextLevel , minPNextLevel1 , maxPNextLevel1 );
	// call recursively for max
	evalNextLevel[index+2] = actDom.getGradedAxisMax(index+2);
	FITOB_OUT_LEVEL3(verb(), "DaOperator::calcForward_rec, call REC MAX:" << (index-1) );
	calcForward_rec( index - 1 , actDom , minVal2, maxVal2 ,
			evalNextLevel , minPNextLevel2 , maxPNextLevel2 );

	// choose the minimum point
	FITOB_OUT_LEVEL3(verb(), "DaOperator::calcForward_rec, TestMIN:" << index << " minVal1: " << minVal1 << " minVal2: " << minVal2);
    if (minVal1 < minVal2){
    	for(int m=0 ; m < index ; m++) minP[m]=minPNextLevel1[m];
    	minP[index] = 0; //min
    	minDaA = minVal1;
    	FITOB_OUT_LEVEL3(verb(), "DaOperator::calcForward_rec, TestMIN Val1");
    } else {
    	for(int m=0 ; m < index ; m++) minP[m]=minPNextLevel2[m];
    	minP[index] = 1; //max
    	minDaA = minVal2;
    	FITOB_OUT_LEVEL3(verb(), "DaOperator::calcForward_rec, TestMIN Val2");
    }
	FITOB_OUT_LEVEL3(verb(), "DaOperator::calcForward_rec, TestMIN Done: "<< (index-1));

    // choose the maximum point
	FITOB_OUT_LEVEL3(verb(), "DaOperator::calcForward_rec, TestMAX:" << index << " maxVal1: " << maxVal1 << " maxVal2: " << maxVal2);
    if (maxVal1 < maxVal2){
    	for(int m=0 ; m < index ; m++) maxP[m]=maxPNextLevel2[m];
    	maxP[index] = 1; //max
    	maxDaA = maxVal2;
    	FITOB_OUT_LEVEL3(verb(), "DaOperator::calcForward_rec, TestMAX Val2");
    } else {
    	for(int m=0 ; m < index ; m++) maxP[m]=maxPNextLevel1[m];
    	maxP[index] = 0; //min
    	maxDaA = maxVal1;
    	FITOB_OUT_LEVEL3(verb(), "DaOperator::calcForward_rec, TestMAX Val1");
    }
    FITOB_OUT_LEVEL3(verb(), "DaOperator::calcForward_rec, TestMAX Done: "<< (index-1));
}


void DaOperator::backwardCalculation(
		boost::ptr_vector<MeshContext>& contextStack ,
		                          FitobCalculator* calc ,
			                      int& stackIndex ,
		                          double& timeStamp ){

		//       create new mesh and set the values on the new mesh
		//       based on Da operator and the old mesh
		//       - create context first with the old mesh (from FitobCalculator) and old Domain (stackIndex)
		//       - create new mesh with the new mesh and with new Domain (stackIndex - 1)
		//       - set the values on the new mesh
		//       - set the new mesh to be the actual mesh in FitobCalculator

    // create "evaluable" object
	DaOperatorEvaluable evalObj(&(contextStack[stackIndex]) , this);

	if (expressionIsConstant_)
	{
		// if expression is constant then copy the old grid
		boost::shared_ptr<MeshBase>& newGrid = contextStack[stackIndex].getMeshPtr();
		contextStack[stackIndex-1].setMesh(newGrid);
	}
	else
	{
		//  create new Grid only when there is really an axis shif and not a constant shift
		boost::shared_ptr<MeshBase> newGrid;
		newGrid = calc->getGridFactory()->createMesh( &(contextStack[stackIndex-1].getDom()) , calc );
		// set the grid to be the grid of this gridcontext
		contextStack[stackIndex-1].setMesh(newGrid);
		// set the values of the new grid using the values of the old grid
		FITOB_OUT_LEVEL1(verb()," DaOperator::backwardCalculation set new values on the mesh");
		newGrid.get()->setValues( &(evalObj) , contextStack[stackIndex].getDom().getAverage() );

		//apply constraints
		DVector globCoord = contextStack[stackIndex].getDom().getAverage();
		contextStack[stackIndex].getMesh()->applyConstraints( &(calc->getScriptModel().get()->constrainOpSeq()) , globCoord );
	}

	// do eventually plotting
	FITOB_OUT_LEVEL1(verb()," DaOperator::backwardCalculation plotting");
	calc->getPlotter()->plot( &(contextStack[stackIndex-1]) , calc->getScriptModel()->getModelName() );
    // decrement the stack index
	stackIndex = stackIndex - 1;
}



void DaOperator::forwardMCSimulation( boost::ptr_vector<MeshContext>& contextStack ,
								  FitobCalculator* calc , MCMachine* MC ,
								  int& stackIndex , double& timeStamp ) {

	FITOB_OUT_LEVEL1(verb()," DaOperator::forwardMCSimulation START timeStamp:" << timeStamp << " , stackIndex:" << stackIndex );
	if (expressionIsConstant_)
	{
		// nothing to do here, a constant expression is already in the domain
	}
	else
	{
		// create the new step, this is not necessary if the expression is constant
		// we do not need to create a new MCStep since we can operate directly on the actual MC step
		// todo: this might be necessary only when the number of variables are changing
		MC->addMCStep( contextStack[stackIndex].getDom() , contextStack[stackIndex].getTime() , stackIndex );

		// apply the Da operator indifferent if the target is a diffusion axis or not
		(MC->getMCStep(MC->nrMCSteps()-1)).applyExpression( &(MC->getMCStep(MC->nrMCSteps()-1)) ,
				assignedVariable_ , expression_ , calc);
	}
	stackIndex = stackIndex + 1;
	FITOB_OUT_LEVEL1(verb()," DaOperator::forwardMCSimulation END timeStamp:" << timeStamp << " , stackIndex:" << stackIndex );
}


void DaOperator::forwardMCEvaluation( boost::ptr_vector<MeshContext>& contextStack ,
								  FitobCalculator* calc , MCMachine* MC ,
								  int& stackIndex , double& timeStamp ) {
	// nothing to do here
	if (expressionIsConstant_)
	{
		// nothing to do here, a constant expression is already in the domain
	}
	else
	{
		MC->getMCStep(stackIndex+1).setExportValues( &(MC->getMCStep(stackIndex)) );
		stackIndex = stackIndex + 1;
	}
}


void DaOperator::backwardMCEvaluation( boost::ptr_vector<MeshContext>& contextStack ,
								  FitobCalculator* calc , MCMachine* MC ,
								  int& stackIndex , double& timeStamp ) {
	// nothing to do here
	// nothing to do here
	if (expressionIsConstant_)
	{
		// nothing to do here, a constant expression is already in the domain
	}
	else
	{
		MC->getMCStep(stackIndex-1).setExportValues( &(MC->getMCStep(stackIndex)) );
		stackIndex = stackIndex - 1;
	}
}

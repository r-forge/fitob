/*
 * FitobScriptModel.cpp
 *
 *  Created on: Feb 28, 2010
 *      Author: benk
 */

#include "FitobScriptModel.hpp"

#include  <boost/algorithm/string/trim.hpp>
#include  <boost/algorithm/string/predicate.hpp>
#include  <boost/lexical_cast.hpp>

#include "src/operators/FitobDaOperator.hpp"
#include "src/operators/FitobIfOperator.hpp"
#include "src/operators/FitobLoopOperator.hpp"
#include "src/operators/FitobMeshAssigmentOperator.hpp"
#include "src/operators/FitobThetaOperator.hpp"
#include "src/operators/FitobNoOperator.hpp"
#include "src/operators/FitobPlotOperator.hpp"

#include "src/diffusionmodel/FitobFactorModelBase.hpp"
#include "src/diffusionmodel/FitobFactorModelGeneral.hpp"

#include "src/expressions/FitobAddExpression.hpp"
#include "src/expressions/FitobAndExpression.hpp"
#include "src/expressions/FitobConstantExpression.hpp"
#include "src/expressions/FitobCondExpression.hpp"
#include "src/expressions/FitobDivExpression.hpp"
#include "src/expressions/FitobEqExpression.hpp"
#include "src/expressions/FitobExpExpression.hpp"
#include "src/expressions/FitobGteExpression.hpp"
#include "src/expressions/FitobGtExpression.hpp"
#include "src/expressions/FitobLogExpression.hpp"
#include "src/expressions/FitobDiscountExpression.hpp"
#include "src/expressions/FitobExpectedExpression.hpp"
#include "src/expressions/FitobLteExpression.hpp"
#include "src/expressions/FitobLtExpression.hpp"
#include "src/expressions/FitobMaxExpression.hpp"
#include "src/expressions/FitobMinExpression.hpp"
#include "src/expressions/FitobMultExpression.hpp"
#include "src/expressions/FitobNeqExpression.hpp"
#include "src/expressions/FitobOrExpression.hpp"
#include "src/expressions/FitobPowExpression.hpp"
#include "src/expressions/FitobSqrtExpression.hpp"
#include "src/expressions/FitobSubsExpression.hpp"
#include "src/expressions/FitobVariableExpression.hpp"

#define MyDebugLevel 0

using namespace fitob;
using namespace std;

ScriptModel::ScriptModel(
		const boost::shared_ptr<XMLConfiguration> xmlconfiguration ,
		tree_parse_info<> scriptTree) :
operatorBody_() , constraintOpSeq_() ,
xmlconfiguration_(xmlconfiguration) , hasConstraints_(false), factorIndex_(0), interestRateModel_(0) {
	nrExportVariables_ = 0;
	importvariables_.resize(0);
	exportVariables_.resize(0);
	scriptDefinedModels_.resize(0);
	setVerb(3);
	FITOB_OUT_LEVEL3( verb() , "Entering ScriptModel Ctor with tree")
	iter_t const modelTree = scriptTree.trees.begin();
	parseTree(modelTree);
}

ScriptModel::ScriptModel(const boost::shared_ptr<XMLConfiguration> xmlconfiguration ) :
operatorBody_() , constraintOpSeq_() ,
xmlconfiguration_(xmlconfiguration) , hasConstraints_(false) , factorIndex_(0), interestRateModel_(0) {
	nrExportVariables_ = 0;
	importvariables_.resize(0);
	exportVariables_.resize(0);
	scriptDefinedModels_.resize(0);
	setVerb(3);
	FITOB_OUT_LEVEL3( verb() , "Entering ScriptModel Ctor without tree")
}

ScriptModel::~ScriptModel() {
	// TODO: write Dtor
	// todo: delete the objects
}

void ScriptModel::parseScript(tree_parse_info<> scriptTree){
	iter_t const modelTree = scriptTree.trees.begin();
	FITOB_OUT_LEVEL3( verb() , "Entering ScriptModel parseTree ")
	parseTree(modelTree);
}

void ScriptModel::parseTree(iter_t const& i){

  // this node has the declarations
  FITOB_OUT_LEVEL3( verb() , "Get ModelBody ");
  iter_t const modelBody = i;
  iter_t const modelName = i->children.begin();

  // the first item must be the script name
  FITOB_OUT_LEVEL3( verb() , " get model Name");
  modelName_ = string(modelName->value.begin(), modelName->value.end());
  FITOB_OUT_LEVEL3( verb() , "Model name: " << modelName_);
  FITOB_OUT_LEVEL3( verb() , "Model Body size: " << modelBody->children.size() );

  // create the time variable
  timeVariable_ = Variable("Time", Time , 0 );

  // the index of each risk factor among the risk factors
  int factorIndex = 0;

  // the 1 .. N-1 nodes here must be declaration of import/export variables
  for (unsigned int ii = 1 ; ii < (modelBody->children.size()-1) ; ii++){
	  iter_t const variableDecl = modelBody->children.begin()+ii;

      if (variableDecl->value.id() == ThetaParser::DeclarationID){
         if ( boost::algorithm::starts_with(
        		  string(variableDecl->value.begin(), variableDecl->value.end()) ,
        		  string("import") )
          ){ // we have an import variable
        	 iter_t const variableName = variableDecl->children.begin();

        	 // THE FIRST TWO VARIABLES
        	 const std::string nameImp = string(variableName->value.begin(), variableName->value.end());
        	 int globalIndex = (importvariables_.size() + nrExportVariables_ + 1);

        	 FITOB_OUT_LEVEL3( verb() , "Import: ");
        	 FITOB_OUT_LEVEL3( verb() , nameImp);

        	 // add variable, in this case we might do it manually
        	 globalIndex = addVariable(nameImp);
         }
         else if ( boost::algorithm::starts_with(
        		  string(variableDecl->value.begin(), variableDecl->value.end()) ,
        		  string("export") )
          ){ // we have an export variable
        	 iter_t const variableName = variableDecl->children.begin();

        	 // THE FIRST TWO VARIABLES
        	 const std::string nameImp = string(variableName->value.begin(), variableName->value.end());
        	 const int globalIndex = nrExportVariables_+1;
        	 nrExportVariables_ = nrExportVariables_ + 1;

        	 FITOB_OUT_LEVEL3( verb() , "Export: ");
        	 FITOB_OUT_LEVEL3( verb() , nameImp);
        	 const VariableType VarType = Export;
        	 exportVariables_.resize(exportVariables_.size()+1);
        	 // add the extra export variable
        	 exportVariables_[exportVariables_.size()-1] = (new Variable(nameImp , VarType , globalIndex ));
         }
         else {
        	 // throw some error
         	 FITOB_ERROR_MSG(" No import or export type found !")
         }
      }
      else
      if (variableDecl->value.id() == ThetaParser::SDEDeclarationID){
     	    FITOB_OUT_LEVEL3( verb() , "Enter SDEDeclaration: ");

  		    iter_t const assignedVariable = variableDecl->children.begin();

  		    std::string nameVar
  		        = deleteSpaces(string(assignedVariable->value.begin(), assignedVariable->value.end()));

  		    int globalIndex = addVariable(nameVar);
     	    FITOB_OUT_LEVEL3( verb() , "SDEDeclaration: varName: " << nameVar << " Glb.Index:" << globalIndex);

     	    // parse the first two expressions they are needed anyhow
  		    iter_t const modelExp = variableDecl->children.begin()+1;
  		    FITOB_OUT_LEVEL3( verb() , "  SDEDeclaration: parsing expr 1");
  		    ExpressionBasis* exp1 = parseExpression(modelExp->children.begin());
  		    FITOB_OUT_LEVEL3( verb() , "  SDEDeclaration: parsing expr 2");
  		    ExpressionBasis* exp2 = parseExpression(modelExp->children.begin()+1);
  		    // a full model declaration
  		    if (modelExp->value.id() == ThetaParser::modefullID){
                FITOB_OUT_LEVEL3( verb() , "  SDEDeclaration: parsing expr 3");
                ExpressionBasis* exp3 = parseExpression(modelExp->children.begin()+2);
                // create the model
                scriptDefinedModels_.resize(factorIndex_+1);
                scriptDefinedModels_[factorIndex_] = (
                		new FactorModelGeneral( importvariables_[globalIndex-(nrExportVariables_+1)] , factorIndex_ ,exp1, exp2, exp3, 0 ) );
  		    } else
  		    // a full model declaration
  		    if (modelExp->value.id() == ThetaParser::modesimplID){
                // create the model
                scriptDefinedModels_.resize(factorIndex_+1);
                scriptDefinedModels_[factorIndex_] = (
                		new FactorModelGeneral( importvariables_[globalIndex-(nrExportVariables_+1)] , factorIndex_ ,exp1, exp1, exp2, 0 ) );
  		    }
  		    else
  		    // an interest rate model declaration
  		    if (modelExp->value.id() == ThetaParser::modeInterestID) {
                FITOB_OUT_LEVEL3( verb() , "  SDEDeclaration: parsing expr 3");
                ExpressionBasis* exp3 = parseExpression(modelExp->children.begin()+2);
                // create the model
                scriptDefinedModels_.resize(factorIndex_+1);
                scriptDefinedModels_[factorIndex_] = (
                		new FactorModelGeneral( importvariables_[globalIndex-(nrExportVariables_+1)] , factorIndex_ ,exp1, exp1, exp2, exp3 ) );
                interestRateModel_ =  scriptDefinedModels_[factorIndex_];
  		    } else {
  		    	// todo: throw error
  		    	FITOB_ERROR_MSG(" No valid SDE declaration !")
  		    }
  		    factorIndex_ = factorIndex_ + 1;
  		    FITOB_OUT_LEVEL3( verb() , "SDEDeclaration END factorIndex = " << factorIndex_);
      }
      else{
    	  // throw some error
      	  FITOB_ERROR_MSG(" No declaration found !")
      }
  }

  // since there might be several import variables declaration before export variable
  // we run over all import variables and set the global index again, if they do not have the global index es needed
  for (unsigned int impVar = 0 ; impVar < importvariables_.size() ; impVar++){
	  FITOB_OUT_LEVEL3( verb() , "Change from:" << importvariables_[impVar]->getGlobalIndex() << "  to:" << 1+nrExportVariables_+impVar);
	  FITOB_OUT_LEVEL3( verb() , "Change Variable name:" << importvariables_[impVar]->toString());
	  importvariables_[impVar]->setGlobalIndex(1+nrExportVariables_+impVar);
  }

  // model operators
  iter_t const OperatorsBody = modelBody->children.begin()+(modelBody->children.size()-1);

  FITOB_OUT_LEVEL3(verb(),"Root Operator Nr: " << OperatorsBody->children.size() );
  for (unsigned int ii = 0 ; ii < OperatorsBody->children.size() ; ii++){
	  // add one operator
	  FITOB_OUT_LEVEL3( verb() , "-------------Add one operator-----------" );
	  operatorBody_.addOperatorToSequence(
			  parseOperator((OperatorsBody->children.begin()+ii) ));

  }

  //FITOB_OUT_LEVEL3(verb(), constraintOpSeq_.toString());
  FITOB_OUT_LEVEL3(verb(), "BODY: \n" <<operatorBody_.toString());
}

/** parse one operator */
OperatorBasis* ScriptModel::parseOperator(iter_t const& i){

// --------- GO PARSING THE SCRIPT ----------------
	if (i->value.id() == ThetaParser::BodyID){

   	    FITOB_OUT_LEVEL3( verb() , "Enter BodyID: ");
		// parse the operator sequence
		OperatorSequence* opSequqnce = new OperatorSequence();
		for (unsigned int ii = 0 ; ii < i->children.size() ; ii++){
			  iter_t const OperatorsBody = i->children.begin()+ii;
			  // add all the operators to the sequence
			  opSequqnce->addOperatorToSequence(parseOperator(OperatorsBody));
		}
   	    FITOB_OUT_LEVEL3( verb() , "Return OperatorSequence: ");
		return opSequqnce;

	} else if (i->value.id() == ThetaParser::AssigmentOpID){
   	    FITOB_OUT_LEVEL3( verb() , "Enter AssigmentOpID: ");

		iter_t const assignedVariable = i->children.begin();

		std::string nameVar
		  = deleteSpaces(string(assignedVariable->value.begin(), assignedVariable->value.end()));

		int globalIndex = addVariable(nameVar);
   	    FITOB_OUT_LEVEL3( verb() , "AssigmentOpID: varName: " << nameVar << " Glb.Index:" << globalIndex);
		iter_t const assignedExp = i->children.begin()+1;
   	    FITOB_OUT_LEVEL3( verb() , "  AssigmentOpID: assig expression parsing");
		ExpressionBasis* asExp = parseExpression(assignedExp);
   	    FITOB_OUT_LEVEL3( verb() , "  AssigmentOpID: done assig expression parsing" );
        if ( (globalIndex >= 1) && (globalIndex < nrExportVariables_+1) ) {  //  means that the mesh value are set (export variable)
        	MeshAssigmentOperator* meshAssOp =
        			new MeshAssigmentOperator( exportVariables_[globalIndex-1] , asExp );
       	    FITOB_OUT_LEVEL3( verb() , "Return MeshAssigmentOperator: " << asExp->toString() );
        	return meshAssOp;
        } else{ // we have a Da operation, since this assigmen targets an import variable
       	    FITOB_OUT_LEVEL3( verb() , "  AssigmentOpID: Create Da opeator: " << globalIndex);
        	DaOperator* daOp =
        			new DaOperator( xmlconfiguration_ , importvariables_[globalIndex-(nrExportVariables_+1)], asExp ); // Important to subtract 2
       	    FITOB_OUT_LEVEL3( verb() , "Return DaOperator: ");
        	return daOp;
        }

	} else if (i->value.id() == ThetaParser::ThetaOpID){

   	    FITOB_OUT_LEVEL3( verb() , "Enter ThetaOpID: ");
		iter_t const ThetaExpression = i->children.begin();
		ExpressionBasis* thetaExp = parseExpression(ThetaExpression);
		ThetaOperator* theta = new ThetaOperator(thetaExp);
   	    FITOB_OUT_LEVEL3( verb() , "Return ThetaOperator: ");
		return theta;

	} else if (i->value.id() == ThetaParser::CommentOpID){

   	    FITOB_OUT_LEVEL3( verb() , "Enter CommentOpID: ");
        // do nothing
		NoOperator* noop = new NoOperator();
   	    FITOB_OUT_LEVEL3( verb() , "Return NoOperator: ");
		return noop;

	} else if (i->value.id() == ThetaParser::LoopOpID){

   	    FITOB_OUT_LEVEL3( verb() , "Enter LoopOpID: ");
		iter_t const LoopExpression = i->children.begin();
		string varName = string(LoopExpression->value.begin(), LoopExpression->value.end());
		varName = deleteSpaces(varName);
		FITOB_OUT_LEVEL3( verb() , "Loop Var Name: " << varName);
		if (boost::algorithm::equals( varName , string("inf"))){
			// here we have to deal with constraints
			FITOB_OUT_LEVEL3( verb() , " parse loop operators body ");
			OperatorBasis* opBasis = parseOperator(i->children.begin()+1);
            // add one contraint operator, which always must hold
	   	    FITOB_OUT_LEVEL3( verb() , "Add Op. to ContstrainOperatorSequence: ");
			constraintOpSeq_.addOperatorToSequence(opBasis);
			hasConstraints_ = true;
	   	    FITOB_OUT_LEVEL3( verb() , "Return NoOperator (because of loop (inf) ): ");
			NoOperator* noop = new NoOperator();
			return noop;
		} else{
			// this is a simple loop operator
			FITOB_OUT_LEVEL3( verb() , " parse loop expression ");
			ExpressionBasis* loopExp = parseExpression(LoopExpression);
            // get the one or more operators
			FITOB_OUT_LEVEL3( verb() , " parse loop operators body ");
			OperatorBasis* opBasis = parseOperator(i->children.begin()+1);
            // create and return
			LoopOperator* loopOp = new LoopOperator(loopExp , opBasis);
	   	    FITOB_OUT_LEVEL3( verb() , "Return LoopOperator: ");
			return loopOp;
		}

	} else if (i->value.id() == ThetaParser::IfElseOpID){

   	    FITOB_OUT_LEVEL3( verb() , "Enter IfElseOpID: ");
		iter_t const LoopExpression = i->children.begin();
		ExpressionBasis* IfExp = parseExpression(i->children.begin());
        // get the one or more operators
		FITOB_OUT_LEVEL3( verb() , "   IfElseOpID:parse true branch: ");
		OperatorBasis* opTrue = parseOperator(i->children.begin()+1);
		FITOB_OUT_LEVEL3( verb() , "   IfElseOpID:parse false branch: ");
		OperatorBasis* opFalse = parseOperator(i->children.begin()+2);
		FITOB_OUT_LEVEL3( verb() , "   IfElseOpID:done parsing true and false branches ");
        // create and return
		IfOperator* ifElseOp = new IfOperator();
		ifElseOp->init( IfExp , opTrue , opFalse);
   	    FITOB_OUT_LEVEL3( verb() , "Return IfOperator (with else): ");
		return ifElseOp;

	} else if (i->value.id() == ThetaParser::IfOpID){

   	    FITOB_OUT_LEVEL3( verb() , "Enter IfOpID: ");
		ExpressionBasis* IfExp = parseExpression(i->children.begin());
        // get the one or more operators
		OperatorBasis* opTrue = parseOperator(i->children.begin()+1);
        // create and return
		IfOperator* ifOp = new IfOperator();
		ifOp->init( IfExp , opTrue );
   	    FITOB_OUT_LEVEL3( verb() , "Return IfOperator: ");
		return ifOp;

	} else if (i->value.id() == ThetaParser::PlOpID){

   	    FITOB_OUT_LEVEL3( verb() , "Enter PlOpID: ");
		ExpressionBasis* plotExp = parseExpression(i->children.begin());
        // create and return
		PlotOperator* plOp = new PlotOperator( plotExp );
   	    FITOB_OUT_LEVEL3( verb() , "Return PlotOperator: ");
		return plOp;

	} else {
		// throw some error
    	FITOB_ERROR_MSG(" No operation found! ")
		NoOperator* noop = new NoOperator();
		return noop;
	}
	FITOB_ERROR_MSG(" No match to operation found! ")
	NoOperator* noop = new NoOperator();
	return noop;
}

/** parse one expression */
ExpressionBasis* ScriptModel::parseExpression(iter_t const& i){

	if (i->value.id() == ThetaParser::constantID){

		 string value = string(i->value.begin(), i->value.end());
		 value = deleteSpaces(value);
		 const double d = boost::lexical_cast<double>(value);
	   	 FITOB_OUT_LEVEL3( verb() ,"Return ConstantExpression: " << value );
	   	 ConstantExpression* cstExp = new ConstantExpression(d);
	   	 FITOB_OUT_LEVEL3( verb() ,"Return ConstantExpression: " << cstExp->toString() );
         return (cstExp);

	} else if (i->value.id() == ThetaParser::factorID){

		iter_t const BaseNode = i->children.begin();
		iter_t const ExpNode = i->children.begin()+1;
		ExpressionBasis* baseExp = parseExpression(BaseNode);
		ExpressionBasis* expExp  = parseExpression(ExpNode);
	   	FITOB_OUT_LEVEL3( verb() ,"Return PowExpression: ");
		return (new PowExpression( baseExp , expExp ));

    } else if (i->value.id() == ThetaParser::termID){

		 string value = string(i->value.begin(), i->value.end());
		 value = deleteSpaces(value);
		 iter_t const BaseNode = i->children.begin();
		 iter_t const ExpNode = i->children.begin()+1;
		 ExpressionBasis* firstExp = parseExpression(BaseNode);
		 ExpressionBasis* secondExp  = parseExpression(ExpNode);
		 // decide which expression to create
         if (boost::algorithm::equals( value , string("*"))){
     	   	 FITOB_OUT_LEVEL3( verb() , "Return MultExpression: ");
             return (new MultExpression(firstExp,secondExp));
         }else{
     	   	 FITOB_OUT_LEVEL3( verb() , "Return DivExpression: ");
             return (new DivExpression(firstExp,secondExp));
         }

    } else if (i->value.id() == ThetaParser::expressionID){

		 string value = string(i->value.begin(), i->value.end());
		 value = deleteSpaces(value);
		 iter_t const firstNode = i->children.begin();
		 iter_t const secondNode = i->children.begin()+1;
 	   	 FITOB_OUT_LEVEL3(verb(), " expressionID  parse first expression");
		 ExpressionBasis* firstExp = parseExpression(firstNode);
 	   	 FITOB_OUT_LEVEL3(verb(), " expressionID  parse second expression");
		 ExpressionBasis* secondExp  = parseExpression(secondNode);
		 // decide which expression to create
        if (boost::algorithm::equals( value , string("+"))){
    	   	FITOB_OUT_LEVEL3( verb() , "Return AddExpression: ");
            return (new AddExpression(firstExp,secondExp));
        }else  if (boost::algorithm::equals( value , string("-"))){
    	   	FITOB_OUT_LEVEL3( verb() , "Return SubsExpression: ");
            return (new SubsExpression(firstExp,secondExp));
        }else  if (boost::algorithm::equals( value , string("<="))){ //&lt;=
    	   	FITOB_OUT_LEVEL3( verb() , "Return LteExpression: ");
            return (new LteExpression(firstExp,secondExp));
        }else  if (boost::algorithm::equals( value , string(">="))){ //&gt;=
    	   	FITOB_OUT_LEVEL3( verb() , "Return GteExpression: ");
            return (new GteExpression(firstExp,secondExp));
        }else  if (boost::algorithm::equals( value , string("=="))){
    	   	FITOB_OUT_LEVEL3( verb() , "Return EqExpression: ");
            return (new EqExpression(firstExp,secondExp));
        }else  if (boost::algorithm::equals( value , string("~="))){
    	   	FITOB_OUT_LEVEL3( verb() , "Return NeqExpression: ");
            return (new NeqExpression(firstExp,secondExp));
        }else  if (boost::algorithm::equals( value , string("&&"))){ //&amp;&amp;
    	   	FITOB_OUT_LEVEL3( verb() , "Return AndExpression: ");
            return (new AndExpression(firstExp,secondExp));
        }else  if (boost::algorithm::equals( value , string("||"))){
    	   	FITOB_OUT_LEVEL3( verb() , "Return OrExpression: ");
            return (new OrExpression(firstExp,secondExp));
        }else  if (boost::algorithm::equals( value , string(">"))){ //"&gt;"
    	   	FITOB_OUT_LEVEL3( verb() , "Return GtExpression: ");
            return (new GtExpression(firstExp,secondExp));
        }else  if (boost::algorithm::equals( value , string("<"))){ //"&lt;"
    	   	FITOB_OUT_LEVEL3( verb() , "Return LtExpression: ");
            return (new LtExpression(firstExp,secondExp));
        }

    } else if (i->value.id() == ThetaParser::variableID){

		 string varName = string(i->value.begin(), i->value.end());
		 bool MClookback = false;
		 varName = deleteSpaces(varName);
		 double signVar = 1.0;
		 if (varName[0] == '-'){
		     // delete the "-" sign
			 varName[0] = ' ';
			 varName = deleteSpaces(varName);
             signVar = -1.0;
         }
		 // delete the "!" if there are any at the end of the variable name
		 int endChar = varName.size()-1;
		 if (varName[endChar] == '!'){
			 MClookback = true;
			 varName[endChar] = ' ';
			 varName = deleteSpaces(varName);
		 }
		 // once we left only with the variable name then add
         int globalIndex = addVariable(varName);
 	   	FITOB_OUT_LEVEL3(verb(),"Return VariableExpression: " << globalIndex << "sign:"<<signVar);
         if (globalIndex == 0) return (new VariableExpression( &timeVariable_ , signVar ));
         if ((globalIndex >= 1) && (globalIndex < 1+nrExportVariables_))
        	 return (new VariableExpression( exportVariables_[globalIndex-1] , signVar , MClookback));
         // the else branch
         return (new VariableExpression( importvariables_[globalIndex-(nrExportVariables_+1)] , signVar , MClookback )); // Important to subtract 2

    } else if (i->value.id() == ThetaParser::endtermID){

    	FITOB_ERROR_MSG(" No end term should be present ! ")
        return parseExpression(i->children.begin());

    } else if (i->value.id() == ThetaParser::minFctID){

		 iter_t const firstNode = i->children.begin();
		 iter_t const secondNode = i->children.begin()+1;
		 ExpressionBasis* firstExp = parseExpression(firstNode);
		 ExpressionBasis* secondExp  = parseExpression(secondNode);
 	   	 FITOB_OUT_LEVEL3(verb(),"Return minFctID: ");
		 return (new MinExpression(firstExp , secondExp ));

    } else if (i->value.id() == ThetaParser::maxFctID){

		 iter_t const firstNode = i->children.begin();
		 iter_t const secondNode = i->children.begin()+1;
		 ExpressionBasis* firstExp = parseExpression(firstNode);
		 ExpressionBasis* secondExp  = parseExpression(secondNode);
 	   	 FITOB_OUT_LEVEL3(verb(),"Return maxFctID: ");
		 return (new MaxExpression(firstExp , secondExp ));

    } else if (i->value.id() == ThetaParser::sqrtFctID){

		 iter_t const firstNode = i->children.begin();
		 ExpressionBasis* firstExp = parseExpression(firstNode);
 	   	 FITOB_OUT_LEVEL3(verb(),"Return sqrtFctID: ");
		 return (new SqrtExpression(firstExp));

    } else if (i->value.id() == ThetaParser::expFctID){

		 iter_t const firstNode = i->children.begin();
		 ExpressionBasis* firstExp = parseExpression(firstNode);
 	   	 FITOB_OUT_LEVEL3(verb(),"Return expFctID: ");
		 return (new ExpExpression(firstExp));

    } else if (i->value.id() == ThetaParser::logFctID){

		 iter_t const firstNode = i->children.begin();
		 ExpressionBasis* firstExp = parseExpression(firstNode);
 	   	 FITOB_OUT_LEVEL3(verb(),"Return logFctID: ");
		 return (new LogExpression(firstExp));

    } else if (i->value.id() == ThetaParser::expectedFctID){

		 iter_t const firstNode = i->children.begin();
		 ExpressionBasis* firstExp = parseExpression(firstNode);
 	   	 FITOB_OUT_LEVEL3(verb(),"Return expectedFctID: ");
		 return (new ExpectedExpression(firstExp));

    } else if (i->value.id() == ThetaParser::discountFctID){

		 iter_t const firstNode = i->children.begin();
		 iter_t const secondNode = i->children.begin()+1;
		 ExpressionBasis* firstExp = parseExpression(firstNode);
		 ExpressionBasis* secondExp  = parseExpression(secondNode);
 	   	 FITOB_OUT_LEVEL3(verb(),"Return discountFctID: ");
		 return (new DiscountExpression(firstExp,secondExp));

    } else if (i->value.id() == ThetaParser::condFctID){

 		 iter_t const firstNode = i->children.begin();
 		 iter_t const secondNode = i->children.begin()+1;
 		 iter_t const thirdNode = i->children.begin()+2;
 		 ExpressionBasis* firstExp = parseExpression(firstNode);
 		 ExpressionBasis* secondExp  = parseExpression(secondNode);
 		 ExpressionBasis* thirdExp  = parseExpression(thirdNode);
  	   	 FITOB_OUT_LEVEL3(verb(),"Return condFctID: ");
 		 return (new CondExpression(firstExp,secondExp,thirdExp));

    }else{

		// throw some error
    	FITOB_ERROR_MSG(" No expression found! ")
	}
	FITOB_ERROR_MSG(" No Match to Expression found!!! ")
    ConstantExpression* cexp = new ConstantExpression(0.0);
    return cexp;
}

// ------------------------ get
const Variable* ScriptModel::getVariable(int globalIndex){
	if (globalIndex == 0) return &timeVariable_;
	if ((globalIndex>0) && (globalIndex < nrExportVariables_+1))  return exportVariables_[globalIndex-1];
	else  return importvariables_[globalIndex-(nrExportVariables_+1)];

}

// ------------------------- Add one variable ------------------ //
int ScriptModel::addVariable(const std::string& variableName){

	std::string nameImp = string(variableName);
	nameImp = deleteSpaces(nameImp);
	//FITOB_OUT_LEVEL3(verb(),"ScriptModel::addVariable VarName:" << nameImp << " variableName:" << variableName << "nrExportVariables:"<<nrExportVariables_);
	// test for the Time variable and for export variable
	if (boost::algorithm::equals( nameImp , string("Time")))
		return 0;

	for(unsigned int i = 0 ; i < exportVariables_.size() ; i++){
		if (boost::algorithm::equals( variableName , exportVariables_[i]->getVariableName()))
			return i+1; // return the global index
	}

	for(unsigned int i = 0 ; i < importvariables_.size() ; i++){
		if (boost::algorithm::equals( nameImp , importvariables_[i]->getVariableName()))
			return i+nrExportVariables_+1; // return the global index = (import var index + nrExportVariables_)
	}

   	const int globalIndex = (importvariables_.size() + nrExportVariables_+1);
   	const VariableType VarType = Import;
   	importvariables_.resize(importvariables_.size()+1);
   	importvariables_[importvariables_.size()-1] = (new Variable(nameImp , VarType , globalIndex ));
   	return globalIndex;
}


const Variable* ScriptModel::addVariableVar(const std::string& variableName) {
	std::string nameImp = string(variableName);
	nameImp = deleteSpaces(nameImp);
	//FITOB_OUT_LEVEL3(verb(),"ScriptModel::addVariable VarName:" << nameImp << " variableName:" << variableName << "nrExportVariables:"<<nrExportVariables_);
	// test for the Time variable and for export variable
	if (boost::algorithm::equals( nameImp , string("Time")))
		return &timeVariable_;

	for(unsigned int i = 0 ; i < exportVariables_.size() ; i++){
		if (boost::algorithm::equals( variableName , exportVariables_[i]->getVariableName()))
			return exportVariables_[i]; // return the global index
	}

	for(unsigned int i = 0 ; i < importvariables_.size() ; i++){
		if (boost::algorithm::equals( nameImp , importvariables_[i]->getVariableName()))
			return importvariables_[i];
	}

   	const int globalIndex = (importvariables_.size() + nrExportVariables_+1);
   	const VariableType VarType = Import;
   	importvariables_.resize(importvariables_.size()+1);
   	importvariables_[importvariables_.size()-1] = (new Variable(nameImp , VarType , globalIndex ));
   	return importvariables_[importvariables_.size()-1];
}

// ------------------- return global index of one variable ---------------------
int ScriptModel::getGlobalVariableIndex(const std::string& variableName) const {

	// test for the Time variable and for export variable
	if (boost::algorithm::equals( variableName , string("Time")))
		return 0;

	for(unsigned int i = 0 ; i < exportVariables_.size() ; i++){
		if (boost::algorithm::equals( variableName , exportVariables_[i]->getVariableName()))
			return i+1; // return the global index
	}

	for(unsigned int i = 0 ; i < importvariables_.size() ; i++){
		if (boost::algorithm::equals( variableName , importvariables_[i]->getVariableName()))
			return i+nrExportVariables_+1; // return the global index
	}
	return -1; //Wall
}

// --------------------- return the import variable index -----------------
int ScriptModel::getImportVariableIndex(const std::string& importName) const {

	for(unsigned int i = 0 ; i < importvariables_.size() ; i++){
		if (boost::algorithm::equals( importName , importvariables_[i]->getVariableName()))
			return i; // variable found return the import index
	}
	return -1; //Wall
}

// --------------------- return the export variable index -----------------
int ScriptModel::getExportVariableIndex(const std::string& importName) const {

	for(unsigned int i = 0 ; i < exportVariables_.size() ; i++){
		if (boost::algorithm::equals( importName , exportVariables_[i]->getVariableName()))
			return i; // variable found return the import index
	}
	return -1; //Wall
}

// ---- clean string from empty spaces -----
std::string ScriptModel::deleteSpaces(string input){
	// delete empty spaces from the string
	boost::algorithm::trim_right(input);
	boost::algorithm::trim_left(input);
	return input;
}

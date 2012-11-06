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
 * FitobScriptModel.hpp
 *
 *  Created on: Feb 28, 2010
 *      Author: benk
 */

#ifndef FITOBSCRIPTMODEL_HPP_
#define FITOBSCRIPTMODEL_HPP_

#include <boost/spirit/include/classic_core.hpp>
#include <boost/spirit/include/classic_ast.hpp>
#include <boost/spirit/include/classic_tree_to_xml.hpp>
#include "src/utils/fitobdefs.hpp"
#include "src/parser/FitobScriptParser.hpp"
#include "src/parser/fitob_script_grammar.hpp"

#include "src/variables/FitobVariable.hpp"
#include "src/operators/FitobOperatorBasis.hpp"
#include "src/operators/FitobOperatorSequence.hpp"
#include "src/expressions/FitobExpressionBasis.hpp"
#include "src/utils/FitobXMLConfiguration.hpp"

#include <iostream>
#include <stack>
#include <functional>
#include <string>
#include <cassert>
#include <fstream>

/** ------------- */

using namespace BOOST_SPIRIT_CLASSIC_NS;

typedef char const*         iterator_t;
typedef tree_match<iterator_t> parse_tree_match_t;
typedef parse_tree_match_t::tree_iterator iter_t;

namespace fitob {

    //forward declarations
	class FactorModelBase;
	class FactorModelGeneral;


using namespace std;

  /** Class containing all the information , model resulted from the script <br>
   * The informations are mainly the declarations and the operator sequence <br>
   *
   * This object stores the variables and its name. There are different sets of variables. <br>
   *  -First set of variable is the global variables, where all variables are found. <br>
   *  0th index is the time variable, the second is the export variable , and the rest are the
   *  import or the internal variables (now internal and import variables are all considered same) <br> <br>
   *  -Second set of variables are the so called "local" variables (or coordinates), this set is the collection of
   *  variables at one given time which are axis of the mesh (are not constant) <br> <br>
   *  -Third set is the set of import variables, (or domain variables), these define the domain. <br>*/
  class ScriptModel : VerbClass{
  public:

	/** The Ctor whith argument , so that the script will be right away parsed*/
	ScriptModel( const boost::shared_ptr<XMLConfiguration> xmlconfiguration ,
			tree_parse_info<> scriptTree );

	/** The Ctor whith no argument */
	ScriptModel( const boost::shared_ptr<XMLConfiguration> xmlconfiguration );

	virtual ~ScriptModel();

	/** */
	void parseScript(tree_parse_info<> scriptTree);

	/** This function checks if the variable exists <br>
	 * if no, then adds as IMPORT variable and returns the global index
	 * if yes, then just returns the global index of the variable* <br>
	 * spaces are deleted from the input string  <br>
	 * returns the global index ! */
	const Variable* addVariableVar(const std::string& variableName);

	/** return the variable of one global index */
	const Variable* getVariable(int globalIndex);

	/** returns the global variable index */
	int getGlobalVariableIndex(const std::string& variableName) const ;

	/** returns the total number of global variables (incl. export and time)*/
	inline int getNrGlobalVariables() const { return importvariables_.size() + nrExportVariables_ + 1 ; }

	/** returns the import variable index (the import index not the global index)*/
	int getImportVariableIndex(const std::string& importName) const ;

	/** returns the import variable name */
	const string& getImportVariableName(int importIndex) const { return importvariables_[importIndex]->getVariableName(); }

	/** returns the Export variable index (the import index not the global index)*/
	int getExportVariableIndex(const std::string& importName) const ;

	/** returns the export variable name */
	const string& getExportVariableName(int exportIndex) const { return exportVariables_[exportIndex]->getVariableName(); }

	/** returns the total number of import variables */
	inline int getNrImportVariables() const { return importvariables_.size(); }

	/** returns the total number of export variables */
	inline int getNrExportVariables() const { return nrExportVariables_; }

	/** returns the starting operation sequence of the script */
	OperatorSequence& bodyOpSeq() {return operatorBody_;}

	/** returns the constraining operation sequence of the script */
	OperatorSequence& constrainOpSeq() {return constraintOpSeq_;}

	/** returns the name of the model*/
	const string& getModelName() const { return modelName_; }

	/** The toString method */
	string toString() { return (constraintOpSeq_.toString() + operatorBody_.toString()); }

	/** number of models that are defined in the script */
	const int getNrDefinedModel() const { return factorIndex_; };

	/** return the model defined in the script */
	const FactorModelBase* getDefinedModel(int index) const { return scriptDefinedModels_[index]; }

	/** signals if there has been defined an interest rate model or not */
	bool hasInterestRateModel() const { return (interestRateModel_ != 0); };

	/** returns the interest rate model defined in the script*/
	const FactorModelBase* getInterestRateModel() const { return interestRateModel_; }

  private:

	/** This function checks if the variable exists <br>
	 * if no, then adds as IMPORT variable and returns the global index
	 * if yes, then just returns the global index of the variable* <br>
	 * spaces are deleted from the input string  <br>
	 * returns the global index ! */
	int addVariable(const std::string& variableName);


	/** This function parses the result from the script <br>
	 * It has to build the declarations <br>
	 * the operator sequence
	 * the expressions*/
	void parseTree(iter_t const& i);

	/** Parse one operator (recursively) */
	OperatorBasis* parseOperator(iter_t const& i);

	/** This function only parses the expression which is needed */
	ExpressionBasis* parseExpression(iter_t const& i);

	/** This is a function which deletes the spaces from names */
	std::string deleteSpaces(string input);

	/** The time variable (global Index = 0)*/
	Variable timeVariable_;

	/** The only single export variable (global Index = 1 .. N) , in case of PDE is for now N=1 */
	std::vector<Variable*> exportVariables_;

	/** count the nr of export variables*/
	int nrExportVariables_;

	/** The list of the import (and) internal variables*/
	std::vector<Variable*> importvariables_;

    /** The operators of the script */
	OperatorSequence operatorBody_;

    /** The operator sequences which enforce eventual constraints for the problem */
	OperatorSequence constraintOpSeq_;

	/** The name of the model */
	std::string modelName_;

    /** XMl configuration which should be global */
    const boost::shared_ptr<XMLConfiguration> xmlconfiguration_;

	/** Is true is there is constraint for the problem
	 * otherwise it must be false */
	bool hasConstraints_;

	/** the variable that includes */
	std::vector<FactorModelBase*> scriptDefinedModels_;

	/** number of defined models in the script*/
	int factorIndex_;

	/** the interest rate model */
	FactorModelBase*  interestRateModel_;

  };

}

#endif /* FITOBSCRIPTMODEL_HPP_ */

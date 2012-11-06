/*
 * FitobScriptParser.cpp
 *
 *  Created on: Feb 26, 2010
 *      Author: benk
 */

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/split.hpp>

#include "FitobScriptParser.hpp"

#define BOOST_SPIRIT_DUMP_PARSETREE_AS_XML

#if defined(BOOST_SPIRIT_DUMP_PARSETREE_AS_XML)
#include <map>
#endif
#include <iostream>
#include <stack>
#include <functional>
#include <string>
#include <cassert>
#include <fstream>
#include <sstream>

using namespace fitob;
using namespace std;
using namespace BOOST_SPIRIT_CLASSIC_NS;

ScriptParser::ScriptParser(const string& scripFile) {
        // the parser structure
	    ThetaParser        parser;
        // open the script file
	    ifstream           in(scripFile.data());

	    stringstream stringstr_out (stringstream::out);

	    // test if we could open the script file
	    if(!in){
	    	FITOB_ERROR_MSG("Cannot open file:" << scripFile);
	    }

	    // read in the content of the script file
	    char                str_tmp[25000];
	    string              str;
	    typedef vector< string > split_vector_type;
	    split_vector_type   SplitVec;
	    while(in){
	      in.getline(str_tmp, 25000);      // Delimiter defaults to newline
	      str = str + str_tmp;
	      // delete the comments (everything behind % sign)
	      boost::algorithm::split( SplitVec , str , boost::algorithm::is_any_of("%") );
	      str = SplitVec[0];
	    }
	    in.close();
        // print out the content of the file when necesarry
	    FITOB_OUT_LEVEL1(1, str << "EOF" );
        // parse the script file
	    info_ = ast_parse(str.c_str(), parser);
        // test if the parsing has succeeded
        if (info_.full)
        {
#if defined(BOOST_SPIRIT_DUMP_PARSETREE_AS_XML)
            // dump parse tree as XML
            std::map<parser_id, std::string> rule_names;
            rule_names[ThetaParser::constantID] = "constant";
            rule_names[ThetaParser::factorID] = "factor";
            rule_names[ThetaParser::termID] = "term";
            rule_names[ThetaParser::expressionID] = "expression";
            rule_names[ThetaParser::variableID] = "variable";
            rule_names[ThetaParser::endtermID] = "endterm";
            rule_names[ThetaParser::minFctID] = "minFct";
            rule_names[ThetaParser::maxFctID] = "maxFct";
            rule_names[ThetaParser::sqrtFctID] = "sqrtFct";
            rule_names[ThetaParser::expFctID] = "expFct";
            rule_names[ThetaParser::logFctID] = "logFct";
            rule_names[ThetaParser::expectedFctID] = "expectedFct";
            rule_names[ThetaParser::discountFctID] = "discountFct";
            rule_names[ThetaParser::condFctID] = "condFct";
            rule_names[ThetaParser::BodyID] = "Body";
            rule_names[ThetaParser::OpID] = "Op";
            rule_names[ThetaParser::ThetaOpID] = "ThetaOp";
            rule_names[ThetaParser::AssigmentOpID] = "AssigmentOp";
            rule_names[ThetaParser::CommentOpID] = "CommentOp";
            rule_names[ThetaParser::LoopOpID] = "LoopOp";
            rule_names[ThetaParser::IfOpID] = "IfOp";
            rule_names[ThetaParser::IfElseOpID] = "IfElseOp";
            rule_names[ThetaParser::PlOpID] = "PlOp";
            rule_names[ThetaParser::ScriptID] = "Script";
            rule_names[ThetaParser::DeclarationID] = "Declaration";
            rule_names[ThetaParser::SDEDeclarationID] = "SDEDeclaration";
            rule_names[ThetaParser::SDEModelDefID] = "SDEModelDef";
            rule_names[ThetaParser::modesimplID] = "modesimpl";
            rule_names[ThetaParser::modefullID] = "modefull";
            rule_names[ThetaParser::modeInterestID] = "modeInterest";

            tree_to_xml(stringstr_out, info_.trees, str.c_str(), rule_names);
#endif
            FITOB_OUT_LEVEL1( 1, "parsing succeeded" );
	     } else {
	        FITOB_ERROR_MSG("parsing failed");
	     }

}

ScriptParser::~ScriptParser() {
  // nothing to do
}

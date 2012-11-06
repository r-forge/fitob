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
 * FitobScriptParser.hpp
 *
 *  Created on: Feb 26, 2010
 *      Author: benk
 */

#ifndef FITOBSCRIPTPARSER_HPP_
#define FITOBSCRIPTPARSER_HPP_

#include <boost/spirit/include/classic_core.hpp>
#include <boost/spirit/include/classic_ast.hpp>
#include <boost/spirit/include/classic_tree_to_xml.hpp>
#include "src/parser/fitob_script_grammar.hpp"
#include "src/utils/fitobdefs.hpp"

namespace fitob{

using namespace std;

/** This class reads the script file and with boost spirit it returns into
 * a tree structure which then needs to be parsed to build the final set of instruction */
class ScriptParser {

public:

	/** Ctor with filename */
	ScriptParser(const string& scripFile);

	/** */
	virtual ~ScriptParser();

	/** returns the parsed information */
	const tree_parse_info<>& getParseInfo() const { return info_; }

private:

	/** the parsed (XML) tree parsed with BOOST SPIRIT */
	mutable tree_parse_info<> info_;

};
}

#endif /* FITOBSCRIPTPARSER_HPP_ */

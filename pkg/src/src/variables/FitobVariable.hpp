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
 * FitobVariable.hpp
 *
 *  Created on: Feb 26, 2010
 *      Author: benk
 */

#ifndef FITOBVARIABLE_HPP_
#define FITOBVARIABLE_HPP_

#include "src/utils/fitobdefs.hpp"

typedef enum{ Internal = 0 ,
              Export   = 1 ,
              Import   = 2 ,
              Time     = 3 } VariableType;

namespace fitob{

  using namespace std;

  /** Class to represent all the (types) variables in the script <br>
   * The convention is the following regarding how to store the global variables <br>
   * At the first place is always the Time <br>
   * Second is the export variable <br>
   * The next positions are the import and internal variables  */
  class Variable : public fitob::VerbClass {

  public:

	/** Ctor*/
	Variable();

	/** Ctor for the variable object */
    Variable(const string& name, const VariableType type, const int globalIndex);

	virtual ~Variable();

	/** return the type of the variable */
	inline const VariableType getVariableType() const {return type_;}

	/** return the name of the variable */
	inline const string& getVariableName() const {	return name_; }

	/** test if the Variable is of the required type */
	inline const bool isTypeOF(VariableType type) const {return (type == type_);}

	/** returns the global Index of this variable */
	inline const int getGlobalIndex() const {return globalIndex_;}

	/** sets the global Index for this variable */
	inline void setGlobalIndex(int gI) { globalIndex_ = gI;}

	/** return a string descriptor */
	inline string toString() const { return name_;}

	/** operator to test the equality of two variables */
	inline bool operator==(const Variable& rhs) const {
		return ((rhs.name_ == this->name_) && (rhs.type_ == this->type_));
	};

  private:

	/** name of the variable*/
	mutable string name_;

	/** type of the variable*/
	mutable VariableType type_;

	/** the global index of the variable */
	mutable int globalIndex_;
  };
}

#endif /* FITOBVARIABLE_HPP_ */

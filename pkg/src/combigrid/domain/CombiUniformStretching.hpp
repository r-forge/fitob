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
 * CombiUniformStretching.hpp
 *
 *  Created on: Apr 4, 2011
 *      Author: benk
 */

#ifndef COMBIUNIFORMSTRETCHING_HPP_
#define COMBIUNIFORMSTRETCHING_HPP_

#include "combigrid/domain/AbstractStretchingMaker.hpp"

namespace combigrid {

/** uniform stretching used only for testing purposes */
class UniformStretching : public AbstractStretchingMaker{
public:

	UniformStretching():AbstractStretchingMaker() {;}

	virtual ~UniformStretching() {;}

	void get1DStretching(
			int level , double min, double max,
			std::vector<double>& stretching) const;

};

}

#endif /* COMBIUNIFORMSTRETCHING_HPP_ */

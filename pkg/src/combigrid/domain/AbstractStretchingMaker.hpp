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
 * AbstractStretchingMaker.hpp
 *
 *  Created on: Apr 4, 2011
 *      Author: benk
 */

#ifndef ABSTRACTSTRETCHINGMAKER_HPP_
#define ABSTRACTSTRETCHINGMAKER_HPP_

#include "combigrid/utils/combigrid_ultils.hpp"

namespace combigrid {
/** class to create stretching in 1D*/
class AbstractStretchingMaker{
public:
	/**
	 * @param level [IN] level of the array
	 * @param min [IN] minimum value of the domain
	 * @param max [IN] maximum value of the domain
	 * @param stretching [OUT] 2^level + 1 elements */
	virtual void get1DStretching(
			int level , double min, double max,
			std::vector<double>& stretching) const = 0;

};
}

#endif /* ABSTRACTSTRETCHINGMAKER_HPP_ */

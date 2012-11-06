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
 * CombiAtanSpecialStretching.hpp
 *
 *  Created on: Apr 5, 2011
 *      Author: benk
 */

#ifndef COMBIATANSPECIALSTRETCHING_HPP_
#define COMBIATANSPECIALSTRETCHING_HPP_

#include "combigrid/domain/AbstractStretchingMaker.hpp"

namespace combigrid {

/** Stretching formula with the following matlab formula: <br>
 * 	L = 6; <br>
 * 	N = -1:(1/2^L):1; <br>
 * 	Fakt = 5.3; <br>
 * 	V = tan(((pi/2) - 1/Fakt)*N); <br>
 * 	V1 = 3*atan( (N).^5+0.1*(N) ); <br>
 * 	V = V ./ max(V); <br>
 * 	V1 = V1 ./ max(V1); <br>
*/
class AtanSpecialStretching : public AbstractStretchingMaker {
public:

	AtanSpecialStretching():AbstractStretchingMaker(){;}

	virtual ~AtanSpecialStretching(){;}

	void get1DStretching(
			int level , double min, double max,
			std::vector<double>& stretching) const;
};

}

#endif /* COMBIATANSPECIALSTRETCHING_HPP_ */

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
 * CombiTanStretching.hpp
 *
 *  Created on: Apr 5, 2011
 *      Author: benk
 */

#ifndef COMBITANSTRETCHING_HPP_
#define COMBITANSTRETCHING_HPP_

#include "combigrid/domain/AbstractStretchingMaker.hpp"

namespace combigrid {

/** The stretching function of this class is: <br>
 *octave:28> intFact = 1/4; <br>
  octave:29> x=(-pi/2+intFact):2^-L:(pi/2-intFact); <br>
  octave:30> y = tan(x); <br>
 * */
class TanStretching : public AbstractStretchingMaker{
public:
	/** Ctor
	 * @param intFact must be smaller than one*/
	TanStretching(double intFact = 1.0/7.0 ):AbstractStretchingMaker() , intFact_(intFact){
		if (intFact_ > 1.5) intFact_ = 1.0/1.5;
		if (intFact_ < 0.01) intFact_ = 1.0/10.0;
	}

	virtual ~TanStretching(){;}

	void get1DStretching(
			int level , double min, double max,
			std::vector<double>& stretching) const;

private:

	/** internal factor for the formula */
	double intFact_;

};

}

#endif /* COMBITANSTRETCHING_HPP_ */

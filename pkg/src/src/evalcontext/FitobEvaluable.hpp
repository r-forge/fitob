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
 * FitobEvaluable.hpp
 *
 *  Created on: Jul 3, 2010
 *      Author: benk
 */

#ifndef FITOBEVALUABLE_HPP_
#define FITOBEVALUABLE_HPP_

namespace fitob{

  using namespace std;

  /** The base class for all grids (mesh) */
  class Evaluable {
  public:

	/** classical evaluation function where the input parameters are the coordinates <br>
	 * (has dummy implementation , returns 0.0)
	 * @param globalCoordonates , global coordinates of the point which should be evaluated*/
	virtual double eval(const DVector& globalCoordonates) const = 0;

	/** more performance oriented evaluation , to avoid to much function calls
	 * @param globalCoordonates [in] vector of input coordinates
	 * @param resVector [out] result vector */
	virtual void eval(const std::vector<DVector>& globalCoordonates , DVector& resVector) const = 0;

  };
}

#endif /* FITOBEVALUABLE_HPP_ */

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
 * CombiLinearBasisFunction.hpp
 *
 *  Created on: Feb 21, 2011
 *      Author: benk
 */

#ifndef COMBILINEARBASISFUNCTION_HPP_
#define COMBILINEARBASISFUNCTION_HPP_


#include "combigrid/utils/combigrid_ultils.hpp"
#include "combigrid/basisfunction/CombiBasisFunctionBasis.hpp"

namespace combigrid{

   /** Linear basis function */
   class LinearBasisFunction : public combigrid::BasisFunctionBasis {
   public:

	   /** empty Ctror */
	   LinearBasisFunction() {;}

	   /** first method which returns the contribution of the first point in the 1D cell
	    * @param coord  1D coordonate idealy should be [0,1] but for extrapolation could be different [-1,2]*/
	   virtual double functionEval1(double coord) const {
           return (1.0 - coord);
	   }

	   /** second method which returns the contribution of the second point in the 1D cell
	    * @param coord  1D coordonate idealy should be [0,1] but for extrapolation could be different [-1,2]*/
	   virtual double functionEval2(double coord) const {
		   return (coord);
	   }

	   /** return the default basis function*/
	   static const BasisFunctionBasis* getDefaultBasis() { return defaultBasis_;}

   private:

	   /** default basis function */
	   static const BasisFunctionBasis* defaultBasis_;

   };
}

#endif /* COMBILINEARBASISFUNCTION_HPP_ */

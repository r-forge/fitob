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
 * FitobMultigridFG.hpp
 *
 *  Created on: Jul 12, 2010
 *      Author: benk
 */

#ifndef FITOBMULTIGRIDFG_HPP_
#define FITOBMULTIGRIDFG_HPP_

#include "src/utils/fitobdefs.hpp"
#include "src/mesh/FitobFullGrid.hpp"
#include "src/pdesolver/FitobMultigridFGBase.hpp"

namespace fitob {

   using namespace std;

   /** Special class for multigrid full grid <br> */
   class MultigridFG : public MultigridFGBase {

   public:

	  /** Ctor for the highest level contructor */
	  MultigridFG(const FullGridBase* fullgrid);

	  /** Creates a MGFG for the lower level , if flag is true*/
      MultigridFG(const MultigridFG *mgfg , bool coarse);

      /** empty Dtor */
	  virtual ~MultigridFG() {;}

	  /** if inner points > 10 , then no further refinement is possible */
	  virtual bool caBeRefined() const { return (this->innerPoints() > 10); };

	  /** apply restrictions on the mesh */
	  virtual void applyConstraints(const OperatorSequence* constraintOpSeq_ ,  const DVector& globalCoords );

	  /** plot the grid for debugging purpose */
	  virtual void plotMAT_grad(const string& filename , DVector& vect , int count ) const ;

	  /** return the type of the multi-grid */
	  virtual MultigridFGType getMultigrigFGType() const  { return MG_FG_WITH_BOUNDARY; }
};

}

#endif /* FITOBMULTIGRIDFG_HPP_ */

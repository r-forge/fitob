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
 * FitobMultigridFG_WB.hpp
 *
 *  Created on: Jan 3, 2011
 *      Author: benk
 */

#ifndef FITOBMULTIGRIDFG_WB_HPP_
#define FITOBMULTIGRIDFG_WB_HPP_

#include "FitobMultigridFGBase.hpp"

namespace fitob {

/** Multigrid full grid for the multigrid method <br>.
 * Similar to the fullgrids */
class MultigridFG_WB: public fitob::MultigridFGBase {
public:
	  /** Ctor for the highest level contructor */
	  MultigridFG_WB(const FullGridBase* fullgrid);

	  /** Creates a MGFG for the lower level , if flag is true*/
	  MultigridFG_WB(const MultigridFG_WB *mgfg , bool coarse);

      /** empty Dtor */
	  virtual ~MultigridFG_WB() {;}

	  /** */
	  virtual bool caBeRefined() const { return (maxlevel_ > 2); };

	  /** apply restrictions on the mesh */
	  virtual void applyConstraints(const OperatorSequence* constraintOpSeq_ ,  const DVector& globalCoords );

	  /** plot the grid for debugging purpose */
	  virtual void plotMAT_grad(const string& filename , DVector& vect , int count ) const ;

	  /** return the type of the multigrid */
	  virtual MultigridFGType getMultigrigFGType() const  { return MG_FG_WITHOUT_BOUNDARY; }

private:

	  /** the maximal level*/
	  int maxlevel_;
};

}

#endif /* FITOBMULTIGRIDFG_WB_HPP_ */

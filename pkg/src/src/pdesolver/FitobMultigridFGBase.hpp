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
 * FitobMultigridFGBase.hpp
 *
 *  Created on: Jan 3, 2011
 *      Author: benk
 */

#ifndef FITOBMULTIGRIDFGBASE_HPP_
#define FITOBMULTIGRIDFGBASE_HPP_

#include "src/utils/fitobdefs.hpp"
#include "src/mesh/FitobFullGridBase.hpp"

namespace fitob {

    /** types of the full grid for multigrid */
    typedef enum{ MG_FG_WITH_BOUNDARY   = 0 ,
	              MG_FG_WITHOUT_BOUNDARY   = 1 } MultigridFGType;

   using namespace std;

   /** Base class for the full grid which is used in the multi grid */
   class MultigridFGBase : public VerbClass {

   public:

	  /** Ctor for the highest level contructor */
	  MultigridFGBase(const FullGridBase* fullgrid);

	  /** Creates a MGFG for the lower level , if flag is true*/
	  MultigridFGBase(const MultigridFGBase *mgfg , bool coarse);

      /** empty Dtor */
	  virtual ~MultigridFGBase() {;}

	  /** writes the solution back to the full grid*/
	  inline void writeBackSolution(FullGridBase* fullgrid){
			// write back
			FITOB_ERROR_TEST( (int)u_.size() == fullgrid->totalPoints() , " MultigridFG::writeBackSolution , sizes must match "
					<< u_.size() << " - " << fullgrid->totalPoints() );
			// assuming identical mesh we write back the solution
			for (int i=0 ; i < fullgrid->totalPoints() ; i++){
				*(fullgrid->valP(i)) = u_[i];
			}
	  }

	  /** dimension of this grid */
	  inline int dim() const { return dim_; }

	  /** level in the required axis*/
	  inline int level(int i) const { return levels_[i]; }

	  /** number of points in one axis*/
	  inline int nrPoints(int i) const { return nrPoints_[i]; }

	  /** function for the smoother  to access the unknowns directly */
	  inline DVector& u() { return u_; }

	  inline const DVector& u() const { return u_; }

	  /** fast access to the right hand side*/
	  inline DVector& rhs() {return rhs_; }

	  inline const DVector& rhs() const { return rhs_; }

	  /** fast access to the residuum vector on the mesh*/
	  inline DVector& res() { return residuum_; }

	  inline const DVector& res() const { return residuum_; }

	  /** fast access to the diagonal entry of the mesh, this is to save FLOP by calculating it. <br>
	   * These coefficients should be calculated before the iteration, by the smoother algorithm. */
	  inline DVector& d() { return d_; }

	  inline const DVector& d() const { return d_; }

	  /** the scaling of the selected axis */
	  inline const DVector& axisScaling(int axis) const { return axisScaling_[axis]; }

	  inline const Domain* domain() const { return domain_; }

	  inline const IVector& offset() const { return offset_; }

	  /** the number of inner points*/
	  int innerPoints() const { return innerPoints_; }

	  /** apply restrictions on the mesh */
	  virtual void applyConstraints(const OperatorSequence* constraintOpSeq_ ,  const DVector& globalCoords ) = 0;

	  /** plot the grid for debugging purpose */
	  virtual void plotMAT_grad(const string& filename , DVector& vect , int count ) const = 0;

	  /** for the multi-gird to know how deep the V-cycle , which grid can be further refined*/
	  virtual bool caBeRefined() const = 0;

	  /** return the type of the multigrid */
	  virtual MultigridFGType getMultigrigFGType() const = 0;

	  /** create a multigrid-FG out af a full grid*/
	  static MultigridFGBase* createMultigridFG(FullGridBase* fg);

	  /** create a multigrid-FG out af a multigrid-FG and makes it corser (refine strategy for the multigrid)*/
	  static MultigridFGBase* createMultigridFG(MultigridFGBase* mgfg , bool corse);

   protected:

	  /** dimension */
	  int dim_;

	  /** */
	  int innerPoints_;

	  /** Level of the Mesh per axis*/
	  IVector levels_;

      /** The number of points in each direction */
	  IVector nrPoints_;

	  /** Integer values which helps to calculate the stencil */
	  IVector offset_;

	  /** one Double vector for the axis scaling purposes */
	  std::vector<DVector> axisScaling_;

	  /** Value of the right hand side of the equation */
	  DVector rhs_;

	  /** Values of the unknown */
	  DVector u_;

	  /** Values of matrix A_{i,i} values , to save FLOP in one iteration, should be calculated before the iteration */
	  DVector d_;

	  /** Residuum value */
	  DVector residuum_;

	  /** Domain corresponding to the mesh, IMPORTANT: the scaling should not be used here */
	  const Domain* domain_;
};

}

#endif /* FITOBMULTIGRIDFGBASE_HPP_ */

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
 * FitobMeshBase.hpp
 *
 *  Created on: Feb 27, 2010
 *      Author: benk
 */

#ifndef FITOBMESHBASE_HPP_
#define FITOBMESHBASE_HPP_

#include "src/utils/fitobdefs.hpp"
#include "src/evalcontext/FitobDomain.hpp"
#include "src/evalcontext/FitobEvaluable.hpp"

namespace fitob{

  /** types and characteristics of the grid */
  typedef enum{ GRID_WITH_BOUNDARY   = 0 ,
				GRID_WITHOUT_BOUNDARY   = 1 } GridType;

  // forward declaration of this
  class OperatorSequence;

  using namespace std;

  /** The base class for all grids (mesh) */
  class MeshBase : public VerbClass{
  public:

	/** Ctor of the basis class
	 * @param name, name of the mesh type*/
	MeshBase(const Domain* dom , const string& name);

	virtual ~MeshBase() {;}

	/** classical evaluation function where the input parameters are the coordinates <br>
	 * (has dummy implementation , returns 0.0)
	 * @param globalCoords , global coordinates of the point which should be evaluated*/
	virtual double eval(const DVector& globalCoords) const = 0;

	/** more performance oriented evaluation , to avoid to much function calls <br>
	 * both vectors must be presized !! <br>
	 * @param globalCoordonates [in] vector of input coordinates
	 * @param resVector [out] result vector */
	// TODO: implement this in all subclasses
	virtual void eval(const std::vector<DVector>& globalCoordonates , DVector& resVector) const = 0;

	/** This method is used to set the values on the grid
	 * @param func [in] the object which can be evaluable at global coordinates (can be function or mesh in the background)
	 * @param globalCoords [in] global coordinates, time part is needed, since Domain does not contains the time, which
	 * might be needed for function evaluation
	 *   */
	virtual void setValues(const Evaluable* func , const DVector& globalCoords) = 0;

	/** function to apply the constraints to the grid values. <br>
	 * @param constraintOpSeq_ [in] the sequence of operators (which must have special form), representing the contraints
	 * @param globalCoords [in] global coordinates (only needed for the time, if this is used in the constraints)*/
	virtual void applyConstraints(const OperatorSequence* constraintOpSeq_ ,
			                      const DVector& globalCoords ) = 0;

	/** Return the domain belonging to the*/
	const Domain* domain() const {return meshDomain_;}

	/** toString method for debugging */
	virtual const string& toString() const {return name_;}

	/** this grid has boundary points */
	virtual GridType getGridType() const { return gridType_; }

	/** grid specific plotting
	 * @param filename filename with path but without extension, the extension has to be added
	 * by the gird specific plotting */
	virtual void gridSpecificPlot(const string& filename ) const {
		//FITOB_OUT(" WARNING: gridSpecificPlot is not overwritten");
	}

  protected :

	/** */
	void setGridType(GridType gT ) { gridType_ = gT;}

	GridType gridType_;

  private:

	/** the name of the mesh*/
	string name_;

	/** each mesh must have a domain*/
    const Domain* meshDomain_;

  };
}

#endif /* FITOBMESHBASE_HPP_ */

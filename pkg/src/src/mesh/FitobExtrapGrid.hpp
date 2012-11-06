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
 * FitobCombiGrid.hpp
 *
 *  Created on: March 23, 2011
 *      Author: benk
 */

#ifndef FITOBEXTRAPGRID_HPP_
#define FITOBEXTRAPGRID_HPP_

#include "src/utils/fitobdefs.hpp"
#include "src/mesh/FitobSparseGrid.hpp"
#include "src/mesh/FitobMeshBase.hpp"
#include "src/mesh/FitobFullGridBase.hpp"

// forward declaration
namespace combigrid{
  class AbstractCombiGrid;
  class CombiSchemeBasis;
}

namespace fitob {

   using namespace std;

/** Class for extrapolation grid, which does not have the combination method. <br>
 *  Basically there is no projection to the sparse grid space , hierarchization not possible
 * , only the extrapolation technique is working <br>
 * This means that the set of full grids are keept separately */
class ExtrapGrid : public fitob::SparseGrid{
public:

	ExtrapGrid(const Domain* dom , int sparseGridType , double diagonalCutOffLevel ,
			  bool use_Opticom ,
			  const DVector& adaptiveTruncation = DVector(0) ,
			  GridType gridType = GRID_WITH_BOUNDARY );

	virtual ~ExtrapGrid();

	/** see upper class */
	double eval(const DVector& globalCoords) const;

	/** see upper class */
	void eval(const std::vector<DVector>& globalCoordonates , DVector& resVector) const;

	/** see upper class */
	void setValues(const Evaluable* func , const DVector& globalCoords);

	/** see upper class */
	void applyConstraints(const OperatorSequence* constraintOpSeq_ ,
			                      const DVector& globalCoords );

	/** dehierarchize the sparse grid (needed before the combination technique)*/
	void deHierarchize() { /* nothing to do here, because it is pure combi grid*/ ;}

	/** herarchize the sparse grid (after value setting or combi technique)*/
	void hierarchize() { /* nothing to do here, because it is pure combi grid*/ ;}

// ------------ methods needed for the combination technique ------------------

	/** decompose the sparse grid into many full grids (combination technique) <br>
	 * In this case NO OP*/
	void deCompose();

	/** combine the full grids into the sparse grid (combination technique) <br>
	 * In this case NO OP*/
	void reCompose(double &errorIndicator, bool change_coef);

	/** return the nr of Full Grids*/
	int nrFullGrid() const { return nrFullGrids_; }

    FullGridBase& getFG(int index) { return fullGridContainer_[index]; }

    /** return the underlying combigrid */
    const combigrid::AbstractCombiGrid* getAbstractCombiGrid() const { return combinationGrid_; }

    combigrid::AbstractCombiGrid* getAbstractCombiGrid() { return combinationGrid_; }

private:

	/** combi grid */
	combigrid::AbstractCombiGrid* combinationGrid_;

	/** combi scheme for the combi grid */
	combigrid::CombiSchemeBasis* combiScheme_;

	/** nr of full grids for the combi technique */
	int nrFullGrids_;

	/** The full grid container for the combination technique */
	boost::ptr_vector<FullGridBase> fullGridContainer_;

	/** The full grid container for the combination technique */
	std::vector< IVector > fullGridLevels_;

	/** maximal global level, this will be used when there is no dimension adaptivity */
	int global_Max_level_;

	/** type of the sparse grid */
	int sparseGridType_;

	/** this parameter is only used by one special sparse grid (cut off diagonal)*/
	int diagonal_cut_of_level_;
};

}

#endif /* FITOBEXTRAPGRID_HPP_ */

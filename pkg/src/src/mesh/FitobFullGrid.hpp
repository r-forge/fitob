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
 * FitobFullGrid.hpp
 *
 *  Created on: Jul 5, 2010
 *      Author: benk
 */

#ifndef FITOBFULLGRID_HPP_
#define FITOBFULLGRID_HPP_

#include "src/mesh/FitobMeshBase.hpp"
#include "src/mesh/FitobFullGridBase.hpp"

namespace fitob {

 using namespace std;

 /** full mesh (grid) class, this is a very essential class for the whole project <br>
  * This full grid has boundary points. <br> */
 class FullGrid : public FullGridBase {

 public:

	/** standard Ctor, creates mesh as it is defined in the domain */
	FullGrid(const Domain* dom);

	/** special Ctor, where specialLevels might be different than those in the "dom" object <br>
	 * This constructor will probably be used in the combination technique */
	FullGrid(const Domain* dom , IVector& specialLevels);

	/** special Ctor, where specialLevels might be different than those in the "dom" object <br>
	 * This constructor will probably be used in the combination technique
	 * @param dom [in] domain
	 * @param specialLevels [in] the level in each dimension (in case of combi, anisotrop)
	 * @param values vector with the grid values (from the SGpp combi technique)*/
	FullGrid(const Domain* dom , IVector& specialLevels , DVector* values );

	/** Dtor which deletes the dynamically allocated vector*/
	virtual ~FullGrid();

	/** see super class for more docu */
	double eval(const DVector& globalCoords) const;

	/** see super class for more docu */
	void eval(const std::vector<DVector>& globalCoordonates , DVector& resVector) const;

	/** see super class for more docu */
	void setValues(const Evaluable* func , const DVector& globalCoords);

	/** see super class for more docu */
	void applyConstraints(const OperatorSequence* constraintOpSeq_ ,
			                      const DVector& globalCoords );

	/** return information in a string about the full grid*/
	string toStringFG() const;

	/** plot the grid in one file */
	void plotGrid(const string& filename , int count = 0);

	/** this grid has boundary points */
	virtual GridType getGridType() const { return GRID_WITH_BOUNDARY; }

 private:

	/** used by Ctor*/
	void init();

 };

}

#endif /* FITOBFULLGRID_HPP_ */

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
 * FitobSparseGrid.hpp
 *
 *  Created on: Aug 27, 2010
 *      Author: benk
 */

#ifndef FITOBSPARSEGRID_HPP_
#define FITOBSPARSEGRID_HPP_

#include "src/mesh/FitobMeshBase.hpp"
#include "src/mesh/FitobFullGridBase.hpp"

namespace fitob{

  using namespace std;

  /** The base class for all sparse grids <br>
   * Here we define additional methods to the MeshBase class methods */
  class SparseGrid : public MeshBase{
  public:

	/** Ctor of the basis class
	 * @param name, name of the mesh type*/
	SparseGrid(const Domain* dom , const string& name , bool use_Opticom):
		MeshBase(dom,name) , use_Opticom_(use_Opticom)
	{
	}

	virtual ~SparseGrid() {;}

	/** dehierarchize the sparse grid (needed before the combination technique)*/
	virtual void deHierarchize() = 0;

	/** herarchize the sparse grid (after value setting or combi technique)*/
	virtual void hierarchize() = 0;

	/** decompose the sparse grid into many full grids (combination technique) */
	virtual void deCompose() = 0 ;

	/** combine the full grids into the sparse grid (combination technique)
	 * @param errorIndicator [OUT] error indicator for the combination
	 * @param change_coef [IN] only for OPTICOM, if true then applied if the setting in XML allow*/
	virtual void reCompose(double &errorIndicator, bool change_coef) = 0 ;

	/** return the nr of Full Grids*/
	virtual int nrFullGrid() const = 0 ;

	/** returns the fullgrid for the combi technique */
	virtual FullGridBase& getFG(int index) = 0 ;

	/** this method shows if the full grids in the sparse grid can be combined. <br>
	 * Only for the real combination grid should this method return true */
	virtual bool canCombineFullGrids() const { return false; }

  protected:

	/** weather to use opticom for the combination technique */
	bool use_Opticom_;

  };
}

#endif /* FITOBSPARSEGRID_HPP_ */

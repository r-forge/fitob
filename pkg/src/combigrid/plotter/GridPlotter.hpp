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
 * GridPlotter.hpp
 *
 *  Created on: Jun 1, 2011
 *      Author: benk
 */

#ifndef GRIDPLOTTER_HPP_
#define GRIDPLOTTER_HPP_

#include "combigrid/combigrid/AbstractCombiGrid.hpp"

namespace combigrid {

class Evaluable{

public:
	Evaluable(const FullGridD* fg):fg_(fg) , cg_(0) {;}
	Evaluable(const AbstractCombiGrid* cg):fg_(0) , cg_(cg){;}

	inline double eval(std::vector<double>& coords ) const {
		if (fg_ == 0){
			return cg_->eval(coords);
		}else{
			return fg_->eval(coords);
		}
	}
private:
	const FullGridD* fg_;
	const AbstractCombiGrid* cg_;
};

/** plott one grid in a matlab file */
class GridPlotter {
public:

	/** empty Ctor*/
	GridPlotter(){;}

	/** empty Dtor*/
	virtual ~GridPlotter(){;}

	/** plot one Full grid */
	static void plotFullGrid(const std::string& filePath , const FullGridD* fg ,
			std::vector<double>& globalCoord_in , int resolution = 0);

	/** plot one combination grid */
	static void plotCombiGrid(const std::string& filePath , const AbstractCombiGrid* cg ,
			std::vector<double>& globalCoord_in , int resolution = 0);

private:

	static void plotObject(int dim ,
			const std::string& filePath ,
			const Evaluable* obj ,
			const GridDomain* domain ,
			std::vector<double>& globalCoord_in ,
			int resolution);

};

}

#endif /* GRIDPLOTTER_HPP_ */

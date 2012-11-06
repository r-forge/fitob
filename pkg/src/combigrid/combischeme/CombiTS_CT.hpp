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
 * CombiTS_CT.hpp
 *
 *  Created on: Feb 23, 2011
 *      Author: benk
 */

#ifndef COMBITS_CT_HPP_
#define COMBITS_CT_HPP_

#include "combigrid/utils/combigrid_ultils.hpp"
#include "CombiSchemeBasis.hpp"

namespace combigrid {

/** class of the two scale combination scheme (square root CT) <br>*/
class TS_CT : public CombiSchemeBasis{

public:

	/** Ctor
	 * @param dim dimension of the scheme
	 * @param level global level */
	TS_CT( int dim , int level );

	/** Ctor
	 * @param dim
	 * @param levels the level vector for the dimension adaptive case */
	TS_CT( int dim , const std::vector<int>& levels );

	/** Ctor for cases when in specific dimensions no combi should be done
	 * @param dim
	 * @param levels the level vector for the dimension adaptive case
	 * @param makeCombiInDimension */
	TS_CT( int dim , const std::vector<int>& levels ,
	 const std::vector<bool>& makeCombiInDimension );


	/** Ctor for manual steared TS scheme where the user specifies the higher and the lower levels
	 * @param minlevels the min levels
	 * @param maxlevels the max levels  */
	TS_CT( const std::vector<int>& minlevels ,
			 const std::vector<int>& maxlevels  );

private:

};
}

#endif /* COMBITS_CT_HPP_ */

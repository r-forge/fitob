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
 * CombiCombiSchemeBasis.hpp
 *
 *  Created on: Feb 23, 2011
 *      Author: benk
 */

#ifndef COMBICOMBISCHEMEBASIS_HPP_
#define COMBICOMBISCHEMEBASIS_HPP_

#include "combigrid/utils/combigrid_ultils.hpp"
//#include "combigrid/combigridkernel/CombiGridKernel.hpp"

namespace combigrid {

/** base class for any combi scheme. From this class should all the scheme classes derived e.g. S-CT, TS-CT, ... <br>
 * */
class CombiSchemeBasis {

public:

	/** Empty constructor */
	CombiSchemeBasis() {
		dim_ = 0;
		levels_vector_.resize(0);
		levels_.resize(dim_, 0);
		cofficients_.resize(0);
	}

	/** Ctor */
	CombiSchemeBasis(int dim, int level) {
		dim_ = dim;
		levels_vector_.resize(0);
		levels_.resize(dim_, level);
		cofficients_.resize(0);
	}

	/** Ctor */
	CombiSchemeBasis(int dim, const std::vector<int>& levels) {
		dim_ = dim;
		levels_vector_.resize(0);
		levels_ = levels;
		cofficients_.resize(0);
	}

	/** return the dimension of the combi scheme */
	inline int getDim() const {
		return dim_;
	}

	/** number of subsapces */
	inline int getNrSapces() const {
		return levels_vector_.size();
	}

	/** returns the level vector for one subspace */
	inline const std::vector<int>& getLevel(int i) const {
		return levels_vector_[i];
	}

	inline const std::vector<std::vector<int> >& getLevels() const {
		return levels_vector_;
	}

	/** the levels of the full grid which we extrapolate with the combi scheme*/
	inline const std::vector<int>& getMaxLevel() const {
		return levels_;
	}

	/** returns the coefficient for one subspace */
	inline double getCoef(int i) const {
		return cofficients_[i];
	}

	/** returns the coefficient for one subspace */
	inline std::vector<double> getCoef() const {
		return cofficients_;
	}

	/** Method adding a new full grid.
	 * Coefficients and levels are updated and the indices of the changed
	 * levels are returned.
	 */
	std::vector<int> updateScheme(std::vector<std::vector<int> > levelsNew,
			std::vector<double> coef);

	void setCoef(std::vector<double> newCoef);

	void setCoef(int i, double newCoef);

protected:

	/** function which removes the duplicate spaces */
	void removeDuplicates();

	/** the dimension of the scheme */
	int dim_;

	/** the level vector for each space */
	std::vector<std::vector<int> > levels_vector_;

	/** the coefficients for the spaces */
	std::vector<double> cofficients_;
	/** the levels of the full grid which we*/
	std::vector<int> levels_;

};

}

#endif /* COMBICOMBISCHEMEBASIS_HPP_ */

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
 * FitobFullGridBase.hpp
 *
 *  Created on: Jan 1, 2011
 *      Author: benk
 */

#ifndef FITOBFULLGRIDBASE_HPP_
#define FITOBFULLGRIDBASE_HPP_

#include "FitobMeshBase.hpp"

namespace fitob {

 using namespace std;

 /** full mesh (grid) class, this is a very essential class for the whole project <br>
  * This is one base abstract class. There are two types of full grids with and without boundary points. */
 class FullGridBase : public MeshBase {

 public:

	/** standard Ctor, creates mesh as it is defined in the domain */
	FullGridBase(const Domain* dom , const string& name) : MeshBase( dom , name ) ,
	       actualLevels_(), values_(0) , values_created_(true) {;}

	/** Dtor which deletes the dynamically allocated vector*/
	virtual ~FullGridBase() {
		if (values_created_){
			// if we created the vector, only then delete the vector
			(*values_).resize(0); delete values_;
		}
	}

	/** calculates from vector index one linear index */
	inline int getLinearIndex(const IVector& localIndex) const{
		  int ret = 0;
		  // calculate the total offset
		  for (int i = 0 ; i < nrAxis_ ; i++)	ret += axisOffset_[i]*localIndex[i];
		  return ret;
	}

	/** calculates from linear index one vector index */
	inline void getVectorIndex(int ind , IVector& localIndex) const{
		  // calculate the vector offsets
		  for (int i = 0 ; i < nrAxis_ ; i++)	localIndex[i] = (ind % axisOffset_[i]);
	}

	/** this returns the offset of one axis */
	inline int offs(int axisI) const { return axisOffset_[axisI]; }

	/** return the value at one linear index */
	inline double val(int i) const
	{
		return (*values_)[i];
	}

	/** set value at one linear position */
	inline void setVal( int i , double v)
	{
		(*values_)[i] = v;
	}

	/** return pointer to one part of the vector */
	inline double* valP(int i)
	{
		return &((*values_)[i]);
	}

	/** return the unknown vector */
	inline const DVector& unknVect() const { return (*values_); }

	/** the unknown vector (for Combi grid to change the values)*/
	inline DVector& unknVect() { return (*values_); }

	/** sets the unknown vector */
	inline void setUnknVect(DVector* newVect) {
		values_ = newVect; values_created_ = false;
		FITOB_ERROR_TEST( ( (nrTotalPoints_ == (int)(*values_).size()) || ((*values_).size() == 0) ) ,
			    		"If Input vector was given size must match , values.size(): " << (*values_).size() << " nrTotalPoints_" << nrTotalPoints_ );
	}

	/** returns the scaling of one axis
	 * @param i [in] axis index */
	inline const DVector& getScalingAxis(int i) const {
		return axisGrading_[i];
	}

	/** returns the dimension of the grid */
	inline int dim() const { return nrAxis_; }

	/** returns the level of one axis */
	inline int getAxisLevel(int i) const { return actualLevels_[i]; }

    /** nr points per axis i*/
	inline int nrPointsPerAxis(int i) const { return nrPointsPerAxis_[i];}

	/** total point number*/
	inline int totalPoints() const { return nrTotalPoints_; }

	/** return information in a string about the full grid*/
	virtual string toStringFG() const = 0;

	/** plot the grid in one file */
	virtual void plotGrid(const string& filename , int count = 0) = 0;

 protected:

	/** levels of the axis */
	IVector actualLevels_;

    /** the double vector which stores the values , we use direct pointers in
     * order to have optimal performance */
	//double* values_;
	DVector* values_;

	/** this flag shows if Fitob created the vector or it was previously created */
	bool values_created_;

	/** number of points per axis inclusive */
	IVector nrPointsPerAxis_;

	/** offsets of each axis */
	IVector axisOffset_;

	/** grading of each axis */
	std::vector<DVector> axisGrading_;

	/** the length of values_ */
	int nrTotalPoints_;

	/** nr of real axis */
	int nrAxis_;

 };

}

#endif /* FITOBFULLGRIDBASE_HPP_ */

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
 * FitobDomain.h
 *
 *  Created on: Feb 26, 2010
 *      Author: benk
 */

#ifndef FITOBDOMAIN_H_
#define FITOBDOMAIN_H_

#include "src/utils/fitobdefs.hpp"

namespace fitob{

  using namespace std;

  /** Domain stores the informations which are needed besides the values on the Mesh <br>:
   * Mapping from global to local coordinates and from local to global <br>
   * Also contains the values of the constant variables , and the axis of the mesh variables
   * The scaling of the axis should be stored here, so that it can be passed to the mesh (Dirk's code) <br>
   * <br>
   * Domain stores per axis at the maximum level the resolution of the grading of each axis
   * (each axis at the same maximal level) <br>.
   * If one wants to access the grading at a lower level then we just need to multiply the index
   * with an offset, which is easy to calculate <br>
   * This function could be implemented here (index,level)->the graded point on the 1D axis <br>
   * <br>
   * <br>
   * IMPORTANT: this object should also monitor which axes are constant and which are not, and should be able
   * to map between local and global coordinates: <br>
   * - local coordinates are the coordinates which only the mesh see, the mesh is created only for those dimensions
   * which really exist (min and max values are not equal),
   * - global coordinates are where it does not matter which
   * */

  class Domain : public VerbClass {
  public:

	/** The Ctor directly with graded points
	 * @param axisLevel the level of each axis (in the forward estimation they MUST have the same level)
	 * @param axisGrading one Double vector per axis which stores the points positions
	 * (for forward estimation the vectors have the same length)
	 * @param maxRequiredLevel the maximal level */
	Domain( const IVector &axisLevel ,
			const std::vector<DVector> &axisGrading ,
			int maxRequiredLevel ,
			int nrExportVariables = 1);

	/** Ctor with min and max values per axis (in old fashion way) <br>
	 * That means that the division is linear
	 * @param axisLevel, the axis level per axis defined
	 * @param minAxisValues the minimal value per axis
	 * @param maxAxisValues the maximal value per axis
	 * @param nrExportVariables
	 * @param evalPoint */
	Domain( const IVector &axisLevel ,
			const DVector &minAxisValues ,
			const DVector &maxAxisValues ,
			int nrExportVariables = 1 ,
			const DVector* evalPoint = NULL );

	/** Copy Ctor */
	Domain( const Domain& dom);

	/** Copy Ctor */
	Domain( const Domain* dom);

	/** empty Dtor */
	virtual ~Domain() {;}

	/** returns const vector which contains the axis value */
	inline const DVector& getGradedAxis(int axisGlobalIndex) const {
		return axisGrading_[axisGlobalIndex-(nrExportVariable_+1)];
	}

	/** If the domain is scaled , then we can create a non-scaled domain*/
	Domain* getNonScaledDomain() const ;

	/** returns the average values for each axis */
	inline const DVector& getAverage() const { return axisAverageValue_; }

	/** Set new axis, should be used only after transformation <br>
	 * Automatically the internal data structures are initialized */
	void setAxis(DVector& newAxis , int axisGlobalIndex);

	/** returns the minimum value on one given axis, might be needed at the grid evaluation */
	inline double getGradedAxisMin(int axisGlobalIndex) const {
		return axisGrading_[axisGlobalIndex-(nrExportVariable_+1)][0];
	}

	/** returns the minimum value on one given axis, might be needed at the grid evaluation */
	inline double getGradedAxisMax(int axisGlobalIndex) const {
		return axisGrading_[axisGlobalIndex-(nrExportVariable_+1)][ axisGrading_[axisGlobalIndex-(nrExportVariable_+1)].size() - 1 ];
	}

	/** get the Level associated with one axis */
	int getAxisLevel(int axisGlobalIndex) const {
		return axisLevel_[axisGlobalIndex-(nrExportVariable_+1)];
	}

	/** return the level vetor*/
	inline const IVector& getLevelVector() const { return axisLevel_;}

	/** The maximal level among all axis */
	inline int getMaximalLevel()  const { return maxRequiredLevel_; }

	/** The maximal length of the grading vector */
	inline int getMaxGradAxisSize() const { return (int)maxGradAxisSize_; }

	/** Maps from local to global coordinates */
	inline void localToGlobal(const DVector &localCoord , DVector &globalCoord ) const {
		  // some testing if the sizes match
		  FITOB_ERROR_TEST( (globalCoord.size() == (axisAverageValue_.size())) , "Domain::localToGlobal "
				   << globalCoord.size() << " !=" << (axisAverageValue_.size()));
		  // first fill with the average values
		  for (unsigned int jj = nrExportVariable_+1; jj < axisAverageValue_.size() ; jj++) globalCoord[jj] = axisAverageValue_[jj];
		  for (unsigned int jj = 0; jj < local_to_global_IndexMapping_.size() ; jj++)
		  {
			  //FITOB_OUT_LEVEL3(4," Domain::localToGlobal jj:" << jj << " loc->glob:" << local_to_global_IndexMapping_[jj] );
			  //FITOB_OUT_LEVEL3(4," Domain::localToGlobal size:" << local_to_global_IndexMapping_.size() );
			  globalCoord[local_to_global_IndexMapping_[jj]+(nrExportVariable_+1)] = localCoord[jj];
		  }
	}

	/** Maps from global to local coordinates */
	inline void globalToLocal(const DVector &globalCoord , DVector &localCoord ) const {
		  // some testing if the sizes match
		  FITOB_ERROR_TEST( (localCoord.size() == local_to_global_IndexMapping_.size()) , "Domain::globalToLocal , localCoord.size() = " <<
				  localCoord.size() << " , local_to_global_IndexMapping_.size() = " << local_to_global_IndexMapping_.size() << " , DOMAIN=" << this->toString());
		  for (unsigned int jj = 0; jj < local_to_global_IndexMapping_.size() ; jj++)
		  {
			  localCoord[jj] = globalCoord[local_to_global_IndexMapping_[jj]+(nrExportVariable_+1)];
			  //FITOB_OUT_LEVEL3(verb()," Domain::globalToLocal jj:" << jj << " ");
		  }
	}


	/** Maps the axis index from local to global */
	int inline localToGlobalIndex(int localIndex) const {
		return local_to_global_IndexMapping_[localIndex]+(nrExportVariable_+1);
	}

	/** Maps the axis index from global to local */
	int inline globalToLocalIndex(int axisGlobalIndex) const {
		return global_to_local_IndexMapping_[axisGlobalIndex-(nrExportVariable_+1)];
	}

	/** returns */
	int inline globalToImport(int globalVarIndex) const { return globalVarIndex - (nrExportVariable_+1); }

	/** returns */
	int inline importToGlobal(int importVarIndex) const { return importVarIndex + (nrExportVariable_+1); }

	/** returns */
	int inline exportToGlobal(int exportVarIndex) const { return 1+exportVarIndex; }

	/** returns total number of import variables */
	int inline nrGlobalVariables() const { return nrImportVariable_ + (nrExportVariable_+1); }

	/** returns total number of import variables */
	int inline nrImportVariables() const { return nrImportVariable_; }

	/** returns total number of export variables */
	int inline nrExportVariables() const { return nrExportVariable_; }

	/** returns the number import variables which are not constant (and form axis of the mesh)*/
	int inline nrRealAxis() const { return nrAxis_; }

	/** in the forward calculation we use the maximum level, but
	 * in the backward we might have adaptive levels for dimensions*/
	void applyDimAdaptLevels(const IVector& adaptLevels);

	/** returns the evaluation point in the actual domain ,
	 * these points are in local coordinates
	 * @return the local coordinates */
	const DVector& getEvaluationPoint() const { return evalPoint_; }

	/** Prints out the Domain*/
	const string toString() const;

  private:

	/** Initializes the internal structure , should be called after each "setAxis" */
	void init();

	/** The maximal required level for any axis */
	int maxRequiredLevel_;

	/** The maximal length of the graded vector (in the forward estimation is the same for non constant axis) <br>
	 * */
	unsigned int maxGradAxisSize_;

	/** Total number of import variables on the domain (const and non-const variable)*/
	int nrImportVariable_;

	/** number of non constant import variables*/
	int nrAxis_;

	/** number of export variables, in case of PDE this is for now 1, but for MC can be different */
	int nrExportVariable_;

	/** Vector of vectors, storing the grading of each axis at the same resolution<br>
	 * First index is the axis index second index is the grading index [axisIndex][gradindex] <br>
	 * This object also stores the constant variables , in that case the Vector of that Variable has length 1*/
	std::vector<DVector> axisGrading_;

	/** The average value of one variable (const or not const)*/
	DVector axisAverageValue_;

	/** Level stored per each axis (global axis), level 0 for constant axis */
	IVector axisLevel_;

	/** this stores the mapping from local to global coordinate indexing <br>
	 * [0,1,2] => [1,3,4] */
	IVector local_to_global_IndexMapping_;

	/** this stores the mapping from global to local coordinate indexing <br>
	 * [0,1,2,3,4] => [0,0,1,1,2,3]  */
	IVector global_to_local_IndexMapping_;

	/** final evaluation point in the actual domain */
	DVector evalPoint_;

  };
}

#endif /* FITOBDOMAIN_H_ */

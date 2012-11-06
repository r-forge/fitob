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
 * FitobGridFactory.hpp
 *
 *  Created on: Jul 5, 2010
 *      Author: benk
 */

#ifndef FITOBGRIDFACTORY_HPP_
#define FITOBGRIDFACTORY_HPP_

#include "src/utils/fitobdefs.hpp"

#include "src/mesh/FitobMeshBase.hpp"
#include "src/utils/FitobXMLConfiguration.hpp"
#include "src/evalcontext/FitobDomain.hpp"

namespace fitob {

   // forward declaration
   class FitobCalculator;

   using namespace std;

   /** Class to generate grids based on the configuration */
   class GridFactory : public VerbClass {
   public:

	  /** Ctor which, parses the config files and stores the option*/
	  GridFactory(boost::shared_ptr<XMLConfiguration>& XMLConfig);

	  virtual ~GridFactory() {;}

	  /** creates one mesh , based on configurations */
	  boost::shared_ptr<MeshBase> createMesh(const Domain* gridDomain ,
			                                 const FitobCalculator* calc) const;

	  /** The FITOB calculator for convergence analysis purpose might reset the
	   * dimension adaptive truncation levels */
	  void setDimensionAdaptiveTuncationLevel(DVector &newVect){
		  dimensionAdaptiveTruncationLevels_ = newVect;
		  //if (verb() > 3){
			//  for (int i = 0 ; i < dimensionAdaptiveTruncationLevels_.size() ; i++){
			//	  FITOB_OUT_LEVEL3(verb(), "GridFactory::setDimensionAdaptiveTuncationLevel i:"<<i<<" tr_lev=" <<dimensionAdaptiveTruncationLevels_[i]);
			//  }
		  //}
	  }

   private:

	  /** grids might have common characteristics , this is stored in this variable (e.g. has boundary points or not)*/
	  GridType typeOfGrid_;

	  /** const pointer to the XML configuration */
	  const XMLConfiguration* XMLConfig_;

	  /** switch variable to store the configuration*/
	  int gridtype_;

	  /** this parameter is only uses at one special sparse grid type */
	  double special_diagonal_cut_off_level_;

	  /** weather to use the opticom for combination technique to calculate the grids coefficients */
	  bool use_OptiCom_;

	  /** the flag which shows, weather to use scaling of the grid*/
	  bool use_Scaling_;

	  /** In case of the non scaled meshes we need to create new Domains which
	   * later must be deleted, avoiding memory leak, these domains will be stored here */
	  mutable boost::ptr_vector<Domain> domainDepo_;

	 /** In the case of T-CT me might have dimension adaptive truncation*/
	 DVector dimensionAdaptiveTruncationLevels_;
};

}

#endif /* FITOBGRIDFACTORY_HPP_ */

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
 * FitobMultigirdSolver.hpp
 *
 *  Created on: Jul 6, 2010
 *      Author: benk
 */

#ifndef FITOBMULTIGIRDSOLVER_HPP_
#define FITOBMULTIGIRDSOLVER_HPP_

#include "src/pdesolver/FitobSolverBase.hpp"
#include "src/pdesolver/FitobProlongationRestriction.hpp"
#include "src/mesh/FitobFullGridBase.hpp"
#include "src/mesh/FitobFullGrid.hpp"
#include "src/mesh/FitobFullGrid_WB.hpp"

namespace fitob {

    class MultigridFGBase;
	class SmootherBase;
	class SparseGrid;

  using namespace std;

  /** class to solve the PDE on the full grid, either directly on full grid or decompose
   * the sparse grid in full grids and use the combination technique */
  class MultigirdSolver : public SolverBase {
  public:

	MultigirdSolver(const XMLConfiguration* config);

	virtual ~MultigirdSolver() {;}

	/** solves the PDE on the full mesh
	 * @param grid
	 * @param context
	 * @param models
	 * @param timeStep */
	void solvePDE(MeshBase *grid , const MeshContext* context ,
	              const ModelCollection* models ,
	              const FitobCalculator* calc ,
	              double timeStep);

  private :

	/** solves the PDE with Gauss-Seidel on the full grid with boundary points !
	 * @param fgrid
	 * @param context
	 * @param models
	 * @param timeStep */
	void solvePDE_with_GS(FullGridBase *fgrid , const MeshContext* context ,
	                      const ModelCollection* models ,
	                      const FitobCalculator* calc ,
	                      double timeStep);

	/** solves the PDE with geom. multigrid on the full grid with boundary points !
	 * @param fgrid
	 * @param context
	 * @param models
	 * @param timeStep
	 * @param nrFullGrid , in case of combi grid, the FG number for debugging purposes*/
	void solvePDE_with_MG(FullGridBase *fgrid , const MeshContext* context ,
	                      const ModelCollection* models ,
	                      const FitobCalculator* calc ,
	                      double timeStep ,
	                      int nrFullGrid = 0);

	/** Recursive function for the V cycle */
	void recursive_V_cycle(
            boost::ptr_vector<MultigridFGBase> &gridContain ,
			boost::ptr_vector<SmootherBase> &smoothContain ,
			int maxDepth ,
			int actualDepth ,
			int count ,
			DVector &globCoord);

#if defined(FITOB_MPI)

	/** method for the single thread that distributes the load among processors*/
	static void* MPI_distribute_Load(void *obj);

#endif

    /** nr of presmooth iterations */
	int presmooth_;

	/** nr of postsmooth iterations */
	int postsmooth_;

	/** the maximum time between combi*/
	double minCombiDT_;

	/** the maximum time between combi*/
	double maxCombiDT_;

	/** eps value to control the step size */
	double combiStepEps_;

	/** vector to store the global coordinates*/
	DVector globalCoordinates_backup;

	/** In parallel case to sort the processing order of the full grids in case of combi grids*/
	IVector sortedFullGrids_;

	/** for fix step this is the fix step, for time step control the starting value*/
	double minMicroDT_;

	/** valid only in time step control, this is the upper limit */
	double maxMicroDT_;

	/** the epsilon for the time step control*/
	double timeStepControll_Epsilon_;

	/** The tolerance for the solver for one micro time step*/
	double solverEpsilon_;

	/** use the predicted value of the explicit scheme as initial value for the implicit scheme*/
	bool usePredictor_;

	/** apply the time step control formula, to change the time step*/
	bool makeTimeStepControll_;

	/** in the time step control formula either we use the L2 norm or the Inf norm*/
	bool useInfNormForTimeStepControll_;

	/** We might use an under relaxed predictor (in case of instability)*/
	double predictorUnderrelaxCoef_;

};

}

#endif /* FITOBMULTIGIRDSOLVER_HPP_ */

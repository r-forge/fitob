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
 * FitobMCStep.hpp
 *
 *  Created on: Mar 14, 2011
 *      Author: benk
 */

#ifndef FITOBMCSTEP_HPP_
#define FITOBMCSTEP_HPP_

#include "src/utils/fitobdefs.hpp"
#include "src/evalcontext/FitobDomain.hpp"
#include "src/mesh/FitobMeshBase.hpp"
#include "src/variables/FitobVariable.hpp"
#include "src/expressions/FitobExpressionBasis.hpp"

using namespace std;

namespace fitob {

/** Class to contain the whole state of all scenarios at one time step <br>
 * In the vector we store first the export variable and then the local variables
 * for one path of monte-carlo */
class MCStep : public VerbClass{
public:

	/** Ctor */
	MCStep(const Domain& domain, int nrScearions, double time, int evalContextIndex );

	virtual ~MCStep() {;}

	/** return the number of scenario, since int has a limit of ~ 10^9 ,
	 * this should be enough as a limit for scenarios  */
	inline int getNrScenario() const { return nrScearions_; }

	/** return how many non constant variables are per scenario */
	inline int nrVariables() const { return nrVariables_;}

	/** returns the global variable for one variable in the MC step */
	inline int getGlobalIndex( int varIndex ) const { return varIndex_to_globalIndex_[varIndex]; }

	/** returns the variable variable for one global index  in the MC step */
	inline int getVariableIndex( int globalIndex ) const { return globalIndex_to_varIndex_[globalIndex]; }

	/** return the domain of the step */
	const Domain& getDom() const { return domain_;}

	/** returns the global coordinate vector*/
	const DVector& getGlobalCoords() const { return constGlobalCoordinates_; }

	/** sets one value in the simulation vector*/
	inline void setSimulationValue(int nrScen , int variableIndex , double val) {
		simulationVector_[nrScen*nrVariables_ + variableIndex] = val;
	}

	/** gets one value from the simulation vector*/
	inline double getSimulationValue(int nrScen , int variableIndex ) const {
		return simulationVector_[nrScen*nrVariables_ + variableIndex];
	}

	/** calculate the average value of all scenario */
	void calculateAverage();

	/** method to operate on each variables of one MC step
	 * @param sourceMC [in] the MC step, which is the source
	 * @param targetVar [in] variable which (must be local) will be reset
	 * @param expr [in] the expression which modifies the variable
	 * @param calc [in] the calculator, the main central object*/
	void applyExpression(const MCStep* sourceMC ,
			const Variable* targetVar , const ExpressionBasis* expr ,
			const FitobCalculator* calc);

	/** sets the values from the source Monte-Carlo step*/
	void setValues(const MCStep* sourceMC );

	/** sets only the export values (for evaluations)*/
	void setExportValues(const MCStep* sourceMC );

	/** returns the average values for all export and import variables */
	const DVector& getAverage() const { return constGlobalCoordinates_;}

	/** return the mesh of the MC Step */
	inline const MeshBase* getMesh() const { return mesh_.get(); }

	/** return the mesh of the MC Step */
	inline boost::shared_ptr<MeshBase>& getMeshPrt() { return mesh_; }

	/** return the index of the stack in the */
	inline int getEvalContextIndex() const { return evalContextIndex_; }

	/** return the mesh of the context*/
	inline MeshBase* getMesh() { return mesh_.get(); }

	/** Mesh can be set later , or can be get from the EvalContext object using the stored index */
	void setMesh(MeshBase* mesh) { hasValidMesh_ = true; mesh_ = boost::shared_ptr<MeshBase>(mesh); }

	/** Mesh can be set later*/
	void setMesh(boost::shared_ptr<MeshBase>& mesh) { hasValidMesh_ = true; mesh_ = mesh; }

	/** converts the MC step into a string to display */
	const string toString() const;

private:

	/** nr of scenarios in one simulation */
	const int nrScearions_;

	/** number of "variable" variables in one scenario (nr local variables + export variables)*/
	const int nrVariables_;

	/** vector to contain all the values for variables <br>
	 * the size is nrScearions_ X nrVariables_ , where nrVariables_ is the fastest moving index */
	DVector simulationVector_;

	/** the vector with the global coordinates <br>
	 * the first element is the time */
	DVector constGlobalCoordinates_;

	/** the actual domain for this step*/
	const Domain& domain_;

	/** bool variable to show if this step has a valid mesh */
	bool hasValidMesh_;

	/** the index in the stack for the EvalContext object */
	int evalContextIndex_;

    /** The mesh belonging to this step, if there will be any needed
     * todo: we might have several meshes for one step (we might have different regressions
     * for the same step) */
    boost::shared_ptr<MeshBase> mesh_;

    /** mapping between variable index to global Index */
    IVector varIndex_to_globalIndex_;

    /** mapping between global index and variabl index */
    IVector globalIndex_to_varIndex_;
};

}

#endif /* FITOBMCSTEP_HPP_ */

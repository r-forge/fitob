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
 * FitobGridContext.hpp
 *
 *  Created on: Feb 27, 2010
 *      Author: benk
 */

#ifndef FITOBGRIDCONTEXT_HPP_
#define FITOBGRIDCONTEXT_HPP_

#include "src/evalcontext/FitobDomain.hpp"
#include "src/mesh/FitobMeshBase.hpp"
#include "src/utils/fitobdefs.hpp"
#include "src/evalcontext/FitobEvaluable.hpp"

namespace fitob{

  using namespace std;

  /** This class stores all the information which one grid needs*/
  class MeshContext : public Evaluable {

  public:

	MeshContext(const Domain& domain, double actualtime);

	virtual ~MeshContext();

	/** Domain of the evaluation context*/
	inline const Domain& getDom() const { return meshDomain_; }

	/** Time stamp of the context*/
	inline double getTime() const { return timeStamp_; }

	/** return the mesh of the context*/
	inline const MeshBase* getMesh() const { return mesh_.get(); }

	/** return the mesh of the context, smart Ptr*/
	inline boost::shared_ptr<MeshBase>& getMeshPtr() { return mesh_; }

	/** return the mesh of the context*/
	inline MeshBase* getMesh() { return mesh_.get(); }

	/** Mesh can be set later*/
	void setMesh(MeshBase* mesh) { hasValidMesh_ = true; mesh_ = boost::shared_ptr<MeshBase>(mesh); }

	/** Mesh can be set later*/
	void setMesh(boost::shared_ptr<MeshBase>& mesh) { hasValidMesh_ = true; mesh_ = mesh; }

	/** The context exist without mesh, but so only after setting the mesh in the eval context
	 * will the mesh be valid, till then will not be initialized */
	inline const bool hasValidMesh() const { return hasValidMesh_; }

	/** Evaluation vector for constant expressions */
	inline const DVector& minGlobCoord() const { return minGlobalVariables_; }

	/** classical evaluation function where the input parameters are the coordinates <br>
	 * (has dummy implementation , returns 0.0)
	 * @param globalCoordonates , global coordinates of the point which should be evaluated*/
	virtual double eval(const DVector& globalCoordonates) const {
		if (hasValidMesh_) {
			return mesh_->eval(globalCoordonates);
		}
		else {
			FITOB_ERROR_EXIT(" MeshBase:eval 1 no valid mesh" );
			return 0.0;
		}
	}

	/** more performance oriented evaluation , to avoid to much function calls
	 * @param globalCoordonates [in] vector of input coordinates
	 * @param resVector [out] result vector */
	virtual void eval(const std::vector<DVector>& globalCoordonates , DVector& resVector) const {
		if (hasValidMesh_){
			mesh_->eval( globalCoordonates , resVector );
		} else {
			FITOB_ERROR_EXIT(" MeshBase:eval 2 no valid mesh" );
		}
	}

  private:

	/** The domain of the mesh */
    Domain meshDomain_;

    /** Vector stores the minimum values of the import variables, the actual time
     * and 0 at the export variable place. <br>
     * This vector is used by constant expression evaluation. <br>
     * 1st is time <br>
     * 2nd is the export variable <br>
     * rest the import variables in their order */
    DVector minGlobalVariables_;

    /** the actual time when this mesh */
    double timeStamp_;

    /** If the context has a valid mesh, after Ctor the mesh is not valid*/
    bool hasValidMesh_;

    /** The mesh belonging to this context */
    boost::shared_ptr<MeshBase> mesh_;

  };
}

#endif /* FITOBGRIDCONTEXT_HPP_ */

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
 * FitobConstantMesh.hpp
 *
 *  Created on: Apr 20, 2010
 *      Author: benk
 */

#ifndef FITOBCONSTANTMESH_HPP_
#define FITOBCONSTANTMESH_HPP_

#include "FitobMeshBase.hpp"

namespace fitob{

using namespace std;

/** Class for mesh in 0D which is one point (one double value)*/
class ConstantMesh : public fitob::MeshBase {

public:

	ConstantMesh(const Domain* dom , double value);

	/** see super class for more docu */
	double eval(const DVector& globalCoords) const { return constValue_;}

	/** see super class for more docu */
	void eval(const std::vector<DVector>& globalCoordonates , DVector& resVector) const
	{
		for (unsigned int ii = 0 ; ii < resVector.size() ; ii++)
			resVector[ii] = constValue_;
	}

	/** see super class for more docu */
	void setValues(const Evaluable* func , const DVector& globalCoords) {
		DVector tmpVect = globalCoords;
		tmpVect[1] = constValue_;
		constValue_ = func->eval(tmpVect);
	}

	/** see super class for more docu */
	void applyConstraints(const OperatorSequence* constraintOpSeq_ ,
	                      const DVector& globalCoords ) { /*todo */};

	virtual ~ConstantMesh() {;}

private:

	/** constant value of the mesh*/
	double constValue_;

};
}

#endif /* FITOBCONSTANTMESH_HPP_ */

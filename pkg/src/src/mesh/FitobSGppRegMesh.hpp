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
 * FitobSGppRegMesh.hpp
 *
 *  Created on: Mar 6, 2012
 *      Author: benk
 */

#ifndef FITOBSGPPREGMESH_HPP_
#define FITOBSGPPREGMESH_HPP_

#include "FitobMeshBase.hpp"

namespace sg{
   namespace datadriven{

       class LearnerBase;

       class LearnerBaseSP;
   }
}

namespace fitob {

/** The mesh class used only used for SGpp based regression solver (of Alex H.)*/
class SGppRegMesh : public MeshBase {

public:

	/** Constructor */
	SGppRegMesh(const Domain* dom
			    );

	virtual ~SGppRegMesh();

	virtual double eval(const DVector& globalCoords) const;

	virtual void eval(const std::vector<DVector>& globalCoordonates , DVector& resVector) const;

	virtual void setValues(const Evaluable* func , const DVector& globalCoords);

	virtual void applyConstraints(const fitob::OperatorSequence*, const DVector&);

private:

};

}

#endif /* FITOBSGPPREGMESH_HPP_ */

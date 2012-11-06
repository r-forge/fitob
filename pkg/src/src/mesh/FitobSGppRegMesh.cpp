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
 * FitobSGppRegMesh.cpp
 *
 *  Created on: Mar 6, 2012
 *      Author: benk
 */

#include "FitobSGppRegMesh.hpp"

using namespace fitob;

SGppRegMesh::SGppRegMesh( const Domain* dom
	    ): MeshBase(dom,"SGppRegMesh")
	    {
	//FITOB_OUT_LEVEL3( 6 , "SGppRegMesh::SGppRegMesh DOMAIN = " << domain()->toString() );
}

SGppRegMesh::~SGppRegMesh() {

}

double SGppRegMesh::eval(const DVector& globalCoords) const {
	double tmpD = 0.0;
	DVector localCoords(domain()->nrRealAxis());
	domain()->globalToLocal( globalCoords , localCoords);
	//FITOB_OUT_LEVEL3( 6 , "domain()->toString()=" << domain()->toString() );
	//FITOB_OUT_LEVEL3( 6 , "globalCoords[0]=" << globalCoords[0] << " , globalCoords[1]=" << globalCoords[1]);
	//FITOB_OUT_LEVEL3( 6 , "localCoords[0]=" << localCoords[0] );
	// scale the local coordinates to the unit square
	//FITOB_OUT_LEVEL3( 6 , "XOLD(end+1)=" << localCoords[0] << ";");
	for (int i = 0 ; i < domain()->nrRealAxis() ; i++ ) {
		double minV = domain()->getGradedAxisMin(  domain()->localToGlobalIndex(i)  );
		double maxV =domain()->getGradedAxisMax(  domain()->localToGlobalIndex(i)  );
		localCoords[i] = (localCoords[i]- minV) / (maxV - minV);
	}
	//FITOB_OUT_LEVEL3( 6 , "X(end+1)=" << localCoords[0] << "; Y(end+1)= " << tmpD << ";");
	return tmpD;
}

void SGppRegMesh::eval(const std::vector<DVector>& globalCoordonates , DVector& resVector) const {
	int nrPointsEval = (int)resVector.size();
	int nrGlobalCoords = domain()->nrGlobalVariables();
	DVector localCoords( domain()->nrRealAxis() );
	DVector globalCoords_tmp( nrGlobalCoords );

}

void SGppRegMesh::setValues(const Evaluable* func , const DVector& globalCoords){
	// empty implementation, since the values are not set here
}

void SGppRegMesh::applyConstraints(const fitob::OperatorSequence*, const DVector&){
	// empty implementation, since the values are not set here
}

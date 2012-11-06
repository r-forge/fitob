/*
 * FitobGridContext.cpp
 *
 *  Created on: Feb 27, 2010
 *      Author: benk
 */

#include "FitobMeshContext.hpp"

using namespace fitob;
using namespace std;

MeshContext::MeshContext(const Domain& domain, double actualtime) :
meshDomain_(domain) , timeStamp_(actualtime) , hasValidMesh_(false) , mesh_(){

	// create the vector for constant expression evaluation
	minGlobalVariables_.resize(domain.nrImportVariables()+2); // +2 because of time and the export variable
    // we fill in the vector with the minimum value of each axis
	minGlobalVariables_[0] = actualtime; // 1st is the TIME
	minGlobalVariables_[1] = 0; // 2nd is the EXPORT variable
	for (int ii=0 ; ii < domain.nrImportVariables() ; ii++)
		minGlobalVariables_[ii+2] = domain.getGradedAxisMin(ii+2);

}

MeshContext::~MeshContext() {
	// nothing to do here
}

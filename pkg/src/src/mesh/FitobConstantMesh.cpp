/*
 * FitobConstantMesh.cpp
 *
 *  Created on: Apr 20, 2010
 *      Author: benk
 */

#include "FitobConstantMesh.hpp"

using namespace std;
using namespace fitob;

ConstantMesh::ConstantMesh(const Domain* dom , double value): MeshBase(dom,"Constant Mesh"), constValue_(value) {
    // test if the recieved domain is realy a constant domain
	FITOB_ERROR_TEST( dom->nrRealAxis() < 1, " ConstantMesh, Domain is not constant!")

}

/*
 * FitobMeshBase.cpp
 *
 *  Created on: Feb 27, 2010
 *      Author: benk
 */

#include "FitobMeshBase.hpp"

using namespace fitob;
using namespace std;

MeshBase::MeshBase(const Domain* dom , const string& name): name_(name) , meshDomain_(dom) {
}

/*
 * FitobMultigridFGBase.cpp
 *
 *  Created on: Jan 3, 2011
 *      Author: benk
 */

#include "FitobMultigridFGBase.hpp"
#include "src/pdesolver/FitobMultigridFG.hpp"
#include "src/pdesolver/FitobMultigridFG_WB.hpp"

using namespace fitob;
using namespace std;


MultigridFGBase::MultigridFGBase(const FullGridBase* fullgrid) {
	  // todo: source out here code which is common
}


MultigridFGBase::MultigridFGBase(const MultigridFGBase *mgfg , bool coarse) {
	  // todo: source out here code which is common
}

MultigridFGBase* MultigridFGBase::createMultigridFG(FullGridBase* fg){

	switch (fg->getGridType()){
	case (GRID_WITH_BOUNDARY):
		return (new MultigridFG(fg));
	case (GRID_WITHOUT_BOUNDARY):
		return (new MultigridFG_WB(fg));
	}

	return (MultigridFG*)NULL;
}

MultigridFGBase* MultigridFGBase::createMultigridFG(MultigridFGBase* mgfg , bool corse){

	switch (mgfg->getMultigrigFGType()){
	case (MG_FG_WITH_BOUNDARY):
		return (new MultigridFG( dynamic_cast<MultigridFG*>(mgfg) ,corse));
	case (MG_FG_WITHOUT_BOUNDARY):
		return (new MultigridFG_WB(dynamic_cast<MultigridFG_WB*>(mgfg) ,corse));
	}

	return (MultigridFG*)NULL;
}

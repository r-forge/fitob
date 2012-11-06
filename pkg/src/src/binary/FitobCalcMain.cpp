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
 * FitobCalcMain.cpp
 *
 *  Created on: Aug 18, 2010
 *      Author: benk
 */

#include "src/scripteval/FitobCalculator.hpp"
#include "src/mesh/FitobFullGrid.hpp"
#include <iostream>

#include <R.h>
#include <Rdefines.h> 

/** binary to price one script */

extern "C" 
{

SEXP mainScriptPrice(SEXP XMLFile, SEXP ScriptFile){

    // MPI init
    fitob::FITOB_MPI_Init( NULL , NULL );
    std::string XMLFileName(CHAR(STRING_ELT(XMLFile, 0)));
    std::string ScriptFileName(CHAR(STRING_ELT(ScriptFile, 0)));

    fitob::FitobCalculator calc( XMLFileName , ScriptFileName );
    // calculate the script
    double res = calc.evaluateScript();
    SEXP ret;
    PROTECT(ret = NEW_NUMERIC(1));
    NUMERIC_POINTER(ret)[0] = res;
    std::cout << ::std::setprecision( 12 ) << " RESULT = " << res << std::endl;
    // MPI finalize in case
    fitob::FITOB_MPI_Finalize( );
    // return the price
    return(ret);
}

SEXP mainScriptPriceFG(SEXP XMLFile, SEXP ScriptFile, SEXP levelIn) {

    // MPI init in case
    fitob::FITOB_MPI_Init( NULL , NULL );
    std::string XMLFileName(CHAR(STRING_ELT(XMLFile, 0)));
    std::string ScriptFileName(CHAR(STRING_ELT(ScriptFile, 0)));
    fitob::FitobCalculator calc( XMLFileName , ScriptFileName );

    // compute the script
    double res = 0; 
    int plevel_in = INTEGER(levelIn)[0];
    // call the pricing 
    boost::shared_ptr<fitob::FullGrid> fg = calc.evaluateScript_FG(res, plevel_in);
    std::cout << ::std::setprecision( 12 ) << " RESULT = " << res << std::endl;

    const boost::shared_ptr<fitob::ScriptModel> &scriptModel = calc.getScriptModel();

    int dim = fg->dim();
    int meshL = (int) fg->unknVect().size();
    int axisLength = (int) fg->getScalingAxis(0).size();

    SEXP retPrice, retDim, retMeshLength, retMesh, retAxisLength, retMeshAxis, retAxisNames, retList;
    PROTECT(retPrice = NEW_NUMERIC(1));
    NUMERIC_POINTER(retPrice)[0] = res;

    PROTECT(retDim = NEW_NUMERIC(1));
    NUMERIC_POINTER(retDim)[0] = dim;

    PROTECT(retList = allocVector(VECSXP, 7));
    PROTECT(retMeshLength = allocVector(INTSXP, 1));
    PROTECT(retMesh = allocVector(REALSXP, meshL));
    PROTECT(retAxisLength = allocVector(INTSXP, 1));
    PROTECT(retMeshAxis = allocMatrix(REALSXP, dim , axisLength) );
    PROTECT(retAxisNames = allocVector(STRSXP, dim ));

    INTEGER(retMeshLength)[0] = meshL;
    INTEGER(retAxisLength)[0] = axisLength;

    SET_VECTOR_ELT(retList, 0, retPrice);
    SET_VECTOR_ELT(retList, 1, retDim);
    SET_VECTOR_ELT(retList, 2, retMeshLength);
    SET_VECTOR_ELT(retList, 3, retMesh);
    SET_VECTOR_ELT(retList, 4, retAxisLength);
    SET_VECTOR_ELT(retList, 5, retMeshAxis);
    SET_VECTOR_ELT(retList, 6, retAxisNames);

    int i,j;
    // copy the mesh's vector
    for (i = 0; i < meshL; i++){
      REAL(retMesh)[i] = fg->unknVect()[i];
    }

    // copy the axis
    for (i = 0; i < dim; i++){
      for (j = 0; j < axisLength; j++){
          REAL(retMeshAxis)[i + dim*j] = fg->getScalingAxis(i)[j];
      }
    }
    // copy the axis names
    char tmpName[15];
    for (i = 0; i < dim; i++){
    	// todo: the domain is unfortunately null here ... we should change this
    	//fg->domain()->localToGlobalIndex(i);
        const std::string &axisName = scriptModel->getVariable( i + 2 )->getVariableName();
        axisLength = (axisName.size() > 13) ? 13 : axisName.size();
        for (j = 0 ; j < axisLength ; j++){
        	tmpName[j] = axisName[j];
        }
        tmpName[j] = '\0';
        SET_STRING_ELT(retAxisNames, i, mkChar(tmpName));
     }
    // MPI finalize in case
    fitob::FITOB_MPI_Finalize( );

    return(retList);
}

} // extern C

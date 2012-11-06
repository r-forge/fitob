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
 * FitobCalcMonteCarlo.cpp
 *
 *  Created on: Mar 22, 2011
 *      Author: benk
 */
#include "src/scripteval/FitobCalculator.hpp"
#include "src/mesh/FitobFullGrid.hpp"
#include <iostream>

#include <R.h>
#include <Rdefines.h> 

extern "C" 
{

SEXP mainScriptMCPrice(SEXP XMLFile, SEXP ScriptFile)
{
    // MPI init
    fitob::FITOB_MPI_Init( NULL , NULL );
    std::string XMLFileName(CHAR(STRING_ELT(XMLFile, 0)));
    std::string ScriptFileName(CHAR(STRING_ELT(ScriptFile, 0)));

    fitob::FitobCalculator calc( XMLFileName , ScriptFileName );
    // calculate the script
    DVector res = calc.evaluateScript_MonteCarlo();
    std::cout << ::std::setprecision( 12 ) << " RESULT[0] = " << res[0] << std::endl;
    SEXP ret;
    PROTECT(ret = NEW_NUMERIC(1));
    NUMERIC_POINTER(ret)[0] = res[0];
    // MPI finalize in case
    fitob::FITOB_MPI_Finalize( );
    // return the price
    return(ret);
}

}

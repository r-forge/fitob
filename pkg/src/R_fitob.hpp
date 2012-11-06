/*
 * R_fitob.hpp
 *
 *  Created on: Oct 23, 2012
 *      Author: benk
 */

#ifndef R_FITOB_HPP_
#define R_FITOB_HPP_

#include <R.h>
#include <Rdefines.h> 

extern "C" {

SEXP mainScriptPrice(SEXP XMLFile, SEXP ScriptFile);

SEXP mainScriptPriceFG(SEXP XMLFile, SEXP ScriptFile, SEXP levelIn);

SEXP mainScriptMCPrice(SEXP XMLFile, SEXP ScriptFile);

}

#endif /* R_FITOB_HPP_ */

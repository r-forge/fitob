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
 * FitobConvergenceAnaly.hpp
 *
 *  Created on: Aug 7, 2010
 *      Author: benk
 */

#ifndef FITOBCONVERGENCEANALY_HPP_
#define FITOBCONVERGENCEANALY_HPP_

#include "src/scripteval/FitobCalculator.hpp"
#include "src/utils/fitobdefs.hpp"

namespace fitob {

  using namespace std;

/** Class to analyze the converge of one script */

class ConvergenceAnaly : public VerbClass {

  public:

	  /** Ctor of the analyzer tool
	   * @param configurationXMLFile
	   * @param scriptFile
	   * @param lowLevel [in] the lowest level to calculate
	   * @param highLevel [in] the highest level and the reference value */
	ConvergenceAnaly(const string& configurationXMLFile,
	                 const string& scriptFile,
	                 int lowLevel = -1 , int highLevel = -1 );

	virtual ~ConvergenceAnaly() {;}

	/** do the convergence analysis
	 * @param prices [out]
	 * @param conv [out]
	 * @param norm_L2 [out]
	 * @param norm_Inf [out]
	 * @param rateL2 [out]
	 * @param rate_Inf [out]
	 * @param rate_pointw [out] */
	void doAnalysisConvergence( DVector& prices , DVector& conv,
			                    DVector& norm_L2 , DVector& norm_Inf,
			                    DVector& rate_L2 , DVector& rate_Linf, DVector& rate_pointw);


  private:

	/** the lowes level to do the convergence analysis*/
	int lowLevel_;

	/** high level, which will be the reference value */
	int highLevel_;

	/** conf File path + name */
	string confXMLFile_;

	/** script file path + */
	string sciprtFile_;

	/** if we have dimension adaptive run */
	bool dimensionAdaptive_;
};

}

#endif /* FITOBCONVERGENCEANALY_HPP_ */

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
 * FitobGridPlotter.hpp
 *
 *  Created on: Jul 7, 2010
 *      Author: benk
 */

#ifndef FITOBGRIDPLOTTER_HPP_
#define FITOBGRIDPLOTTER_HPP_

#include "src/utils/fitobdefs.hpp"
#include "src/utils/FitobXMLConfiguration.hpp"
#include "src/evalcontext/FitobMeshContext.hpp"

namespace fitob {

using namespace std;

    class ExpressionBasis;
    class FitobCalculator;

/** class to plot a grid <br>
 * only most the first two active dimensions will be taken in consideration , for the
 * rest the average value will be considered */
class GridPlotter : public VerbClass{

public:

	/** Ctor
	 * @param is the configuration object */
	GridPlotter(const boost::shared_ptr<XMLConfiguration>& config);

	/** empty Dtor*/
	virtual ~GridPlotter() {;}

	/** The function which should be called by the operators and based on the the
	 * configurations, grid and domain it plots the current mesh */
	void plot(const MeshContext* meshContext , const string& filename ,
			const ExpressionBasis* expr = 0 , const FitobCalculator* calc = 0);

	/** gets the flag of the plotting*/
	inline bool getActivePlotting() const { return activePlotting_;}

	/** sets the flag of the plotting*/
	inline void setActivePlotting( bool b) const { activePlotting_ = b; }

private :

	/** */
	void plotVTK(const MeshContext* meshContext , const string& filename,
			const ExpressionBasis* expr = 0 , const FitobCalculator* calc = 0);

	/** */
	void plotMAT(const MeshContext* meshContext , const string& filename,
			const ExpressionBasis* expr = 0 , const FitobCalculator* calc = 0);

	/** */
	void plotMAT_grad(const MeshContext* meshContext , const string& filename,
			const ExpressionBasis* expr = 0 , const FitobCalculator* calc = 0);

	/** */
	void plotGNU(const MeshContext* meshContext , const string& filename ,
			const ExpressionBasis* expr = 0 , const FitobCalculator* calc = 0);

	/** */
	void plotGridSpecific(const MeshContext* meshContext , const string& filename);


	int couter_;

	int couter_gsp_;

	int res_x_;

	int res_y_;

	int axisIndex0_;

	int axisIndex1_;

	// if plotting is enabled or not
	mutable bool activePlotting_;

	/** if on the top of the normal plotting we do also grid specific ploting*/
	bool makeGridSpecificPlotting_;

	/** the directory where to dump all the output files */
	static const string outDir_;

};

}

#endif /* FITOBGRIDPLOTTER_HPP_ */

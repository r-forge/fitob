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
 * FitobFactorModelGeneral.hpp
 *
 *  Created on: Aug 7, 2012
 *      Author: benk
 */

#ifndef FITOBFACTORMODELGENERAL_HPP_
#define FITOBFACTORMODELGENERAL_HPP_

#include "src/utils/fitobdefs.hpp"
#include "src/diffusionmodel/FitobFactorModelBase.hpp"
#include "src/utils/FitobXMLConfiguration.hpp"
#include "src/expressions/FitobExpressionBasis.hpp"
#include "FitobModelCollection.hpp"

namespace fitob {

/** General model that the user can define in a general way \br
 *  S1 = MODEL(0.02*S1,0.02*S,0.4*S1);
 *  S2 = MODEL(0.02*S1-S2,0.4*S2);
 *  r = INTEREST(0.03*r,SQRT(0.4),1.0);
 *  */
class FactorModelGeneral : public FactorModelBase {
public:

	FactorModelGeneral(const Variable* variable , int factorIndex, ExpressionBasis* drift,
			 ExpressionBasis* convec, ExpressionBasis* sigma, ExpressionBasis* discount);

	FactorModelGeneral(const FactorModelGeneral* copyFact , int factorIndex);

	virtual ~FactorModelGeneral();

	/** see parent class */
	virtual double convectionCoef(const DVector &vars) const;

	/** see parent class */
	virtual double driftCoef(const DVector &vars) const;

	/** see parent class */
	virtual double diffusionCoef(const DVector &vars) const;

	/** see parent class */
	virtual void forwardEstimation( double initSize,
			                        double timeStep,
			                        double enlargementFactor,
			                        const DVector &averageVarValues,
			                        double &endSize ) const;

	/** see parent class */
	virtual void simulateForward( int nrPaths , double actTime, double timeStep , int nrGlobalVar ,
			DVector& InValues ,
			DVector& OutValues ,
			const DVector& randN ) const;

	/** see parent class */
	virtual double discountFactor( const DVector &vars , double t1, double t2) const;

private:

	/**  the mu fron the SDE */
	const ExpressionBasis* drift_;

	/**  the mu fron the PDE */
	const ExpressionBasis* convec_;

	/**  the sigma fron the PDE and SDE*/
	const ExpressionBasis* sigma_;

	/** the discount factor if this is an interest rate*/
	const ExpressionBasis* discount_;

};

}

#endif /* FITOBFACTORMODELGENERAL_HPP_ */

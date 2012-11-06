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
 * FitobRegressionController.hpp
 *
 *  Created on: Jul 7, 2011
 *      Author: benk
 */

#ifndef FITOBREGRESSIONCONTROLLER_HPP_
#define FITOBREGRESSIONCONTROLLER_HPP_

#include "FitobRegularization.hpp"

namespace fitob {

  using namespace std;

  // forward declaration
  class MCStep;
  class FitobCalculator;

/** class which controlls the regression process. */
class RegressionController {

public:

	/** empty ctor , no private field */
	RegressionController( ) {;}

	/** tells if the regression can be trigerred by an expression (or it has been already)*/
	static bool canRegressionBeTriggerd() { return canRegressionBeTriggerd_;}

	/** this should be called before a MC step is made */
	static void triggerableRegression( const MCStep* mcStep , const Domain* dom , const FitobCalculator* calc) {
		canRegressionBeTriggerd_ = true;  mcStep_= mcStep;  dom_ = dom;  calc_ = calc;
	}

	/** static function which is called before a Monte-Carlo backward or forward step is maid <br>
	 * During this step one expression like EXPECT(X) can be called and then */
	static void startRegression( const ExpectedExpression* expectExpr , const ExpressionBasis* expEval );

	/** it tells the static Regression, that the regularization results can be dumped <br>
	 * in the same way */
	static void doneRegression();

	/** return the mesh which has been generated */
	static inline const MeshBase* getMesh() { return regClass_->getMesh(); }

	static inline Regularization* getRegClass() {return regClass_;}

private:

	/** flag shows if the regression can be trigerred */
	static bool canRegressionBeTriggerd_;

	/** flags shows if the regression can be triggered */
	static bool isRegressionTriggerd_;

	/** static pointer to the class itself, which will do the regression */
	static Regularization* regClass_;

	/** */
	static const MCStep* mcStep_;

	/** */
	static const Domain* dom_;

	/** */
	static const FitobCalculator* calc_;

};

}


#endif /* FITOBREGRESSIONCONTROLLER_HPP_ */

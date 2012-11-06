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
 * FitobThetaOperator.hpp
 *
 *  Created on: Mar 1, 2010
 *      Author: benk
 */

#ifndef FITOBTHETAOPERATOR_HPP_
#define FITOBTHETAOPERATOR_HPP_

#include "FitobOperatorBasis.hpp"
#include "src/expressions/FitobExpressionBasis.hpp"

namespace fitob {

  using namespace std;

  class ThetaOperator: public fitob::OperatorBasis {
  public:

	ThetaOperator(ExpressionBasis* thetaExpression);

	virtual ~ThetaOperator();

	/** function to do forward estimation */
	virtual void forwardEstimation( boost::ptr_vector<MeshContext>& contextStack ,
			                        FitobCalculator* calc ,
			                        int& stackIndex ,
			                        double& timeStamp );

	/** function for diffusion axis size (enlargement) estimation */
	virtual void forwardEstimation_DiffusionEnlarement(
			                          boost::ptr_vector<MeshContext>& contextStack ,
			                          FitobCalculator* calc ,
			                          DVector& gradValues,
			                          int globalVariableIndex,
				                      int& stackIndex ,
			                          double timeStamp ) const;

	/** function to do backward calculation */
	virtual void backwardCalculation( boost::ptr_vector<MeshContext>& contextStack ,
			                          FitobCalculator* calc ,
				                      int& stackIndex ,
			                          double& timeStamp );

	/** function for simulate the scenarios forward */
	virtual void forwardMCSimulation( boost::ptr_vector<MeshContext>& contextStack ,
									  FitobCalculator* calc , MCMachine* MC ,
									  int& stackIndex , double& timeStamp ) ;

	/** function for evaluate the scenarios forward */
	virtual void forwardMCEvaluation( boost::ptr_vector<MeshContext>& contextStack ,
									  FitobCalculator* calc , MCMachine* MC ,
									  int& stackIndex , double& timeStamp ) ;

	/** function for evaluate the scenarios backward */
	virtual void backwardMCEvaluation( boost::ptr_vector<MeshContext>& contextStack ,
									  FitobCalculator* calc , MCMachine* MC ,
									  int& stackIndex , double& timeStamp ) ;

	string toString() { return operationName_ + " ThetaTime: " +
		 thetaExpression_->toString();}

  private:

	/** Function calculates the scaled normal distribution points
	 * @param level [in] level of the resolution (the nr of returned points)
	 * @param stdFaktor [in] standard deviation factor which multiplies the [-1,1] interval
	 * @param scaleIndex [in] internal variable to try different scaling out
	 * @param output [out] vector contains the normal distribution scaling */
	static void returnScaledNormalDistribution(
			  int level ,
			  double stdFactor ,
			  int scaleIndex ,
			  DVector& output);

	/** */
	ExpressionBasis* thetaExpression_;

	/** */
	double thetaTime_;

	/** We might limit the theta time step */
	double maxThetaTime_;

	/** If limit the theta step this will be grater than one */
	int nr_Theta_;
};

}

#endif /* FITOBTHETAOPERATOR_HPP_ */

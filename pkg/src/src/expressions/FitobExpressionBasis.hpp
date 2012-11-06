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
 * FitobExpressionBasis.hpp
 *
 *  Created on: Feb 25, 2010
 *      Author: benk
 */

#ifndef FITOBEXPRESSIONBASIS_HPP_
#define FITOBEXPRESSIONBASIS_HPP_

#include  <boost/shared_ptr.hpp>
#include  <boost/utility.hpp>

#include  "src/utils/fitobdefs.hpp"
#include  "src/evalcontext/FitobMeshContext.hpp"
#include  "src/evalcontext/FitobEvaluable.hpp"

//typedef shared_ptr<ExpressionBasis> ExpSP;

namespace fitob{

   using namespace std;

   /** forward declaration of the Fitob calculator*/
   class FitobCalculator;

  /** This class is the superclass for all expressions <br>
   * Expressions are used in many operations , e.g IF exp ..., A = exp*/
  class ExpressionBasis : public boost::noncopyable , public Evaluable , public VerbClass{
  public:

	ExpressionBasis();

	ExpressionBasis(const string& name);

	virtual ~ExpressionBasis();

	/** if it is constant across all the mesh than this should be true*/
    virtual inline bool isConstantExpression(const MeshContext& context) const {
    	FITOB_ERROR_EXIT("ExpressionBasis::isConstantExpression , must be overwritten! ");
    }

	/** if it is constant across all the mesh than this should be true*/
    virtual inline double eval(const DVector& globalCoordonates) const {
    	FITOB_ERROR_EXIT("ExpressionBasis::eval , must be overwritten! ");
    }

	/** more performance oriented evaluation , to avoid to much function calls
	 * @param globalCoordonates [in] vector of input coordinates
	 * @param resVector [out] result vector */
	virtual void eval(const std::vector<DVector>& globalCoordonates , DVector& resVector) const {
    	FITOB_ERROR_EXIT("ExpressionBasis::eval 2 , must be overwritten! ");
	}

    /** */
    virtual inline string toString() const {
    	FITOB_ERROR_EXIT("ExpressionBasis::toString , must be overwritten! ");
    }

    /** returns the global Fitob Calculator which is needed in some cases for evaluation*/
    static const FitobCalculator* getFitobCalculator() { return fitobCalc_;}

    /** sets the Fitob Calculator for the evaluation */
    static void setFitobCalculator(FitobCalculator* fitCalc) { fitobCalc_ = fitCalc; }

  protected:

     /** The name of the expression */
     string   name_;

  private:

     /** the global Fitob calculator which should be used in case of evaluation*/
     static FitobCalculator* fitobCalc_;

  };
}

#endif /* FITOBEXPRESSIONBASIS_HPP_ */

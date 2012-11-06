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
 * FitobModelCollection.hpp
 *
 *  Created on: Apr 13, 2010
 *      Author: benk
 */

#ifndef FITOBMODELCOLLECTION_HPP_
#define FITOBMODELCOLLECTION_HPP_

#include "src/utils/fitobdefs.hpp"
#include "src/utils/FitobXMLConfiguration.hpp"
#include "src/diffusionmodel/FitobFactorModelBase.hpp"
#include "src/scriptmodel/FitobScriptModel.hpp"
#include "src/variables/FitobVariable.hpp"


namespace fitob {

using namespace std;

   /** Class containing the array of models which form the collection
    * of factors which influence the product price. <br>
    * Each model is responsible for the forward estimation, and to provide
    * the coefficients value of the PDE */
class ModelCollection : public VerbClass {
public:

	/** The Ctor which creates the array of models */
	ModelCollection( const boost::shared_ptr<XMLConfiguration> &ptr , boost::shared_ptr<ScriptModel> &scriptModel);

	/** empty Dtor, the smart pointers will delete itself*/
	virtual ~ModelCollection() {;}

	/** returns the nr of Factors */
	inline int nrFactors() const { return nrFactors_;}

	/** Here we just return one constant reference, the input is the factor index <br>
	 * */
	inline const FactorModelBase& getModel(int factorIndex) const {
		// Todo: check if index is not out of boundary
		return modelContainer_[factorIndex];
	}

	/** return the global constant risk free interest rate
	 * @param vars the global variables */
	inline double r(const DVector &vars) const {
		// return eighter the constant value or the value of the variable
		return (r_index_ < 0 )? r_ : vars[intRateVariable_->getGlobalIndex()];
	}

	/** returns the standard deviation factor for forward estimation <br>
	 * */
	double inline stdFactor(int factIndex) const { return standardDeviationFactors_[factIndex]; }

	/** returns the correlation of two factors */
	inline double corr(int f1, int f2) const {
		return correlations_[f1*nrFactors_ + f2];
	}

	/** returns the correlation matrix */
	inline const DVector& getCorrelationMatrix() const {
		return correlations_;
	}

	/** return the type as a string of one factor */
	inline const string& getFactorType(int factIndex) const {
		return factorTypes_[factIndex];
	}

	/** If the interest rate is a variable this will be a positive number, <br>
	 * it returns -1 otherwise (when the interest rate is constant over the whole period of time)*/
	inline int getInterestRateGlobalIndex() const { return intRateVariable_->getGlobalIndex(); }

	/** returns the discount factor of this model collection. <br> Since the interest rate can be also a factor
	 * the discount might vary locally in the mesh. for this reason we need an expression */
	inline double getDiscountFactor(const DVector &vars, double t1, double t2) const {
		// if we have constant interest rate then just return the exp(-rT) value
		if (r_index_ < 0){
			return ::exp(- r_* (t2-t1));
		} else {
		// in case of interest rate which is a model then just call the method of the model
			return modelContainer_[ r_index_ ].discountFactor(vars, t1, t2);
		}
	}

	/** this functions completes the information of the model collection by using the
	 * information from the scriptModel. In the Ctor the scriptModel is not parsed yet.
	 * e.g. model definitions are now available
	 * @param ptr [IN] the XML config file
	 * @param scriptModel [IN] the parsed script model */
	void completeScriptInformation(const boost::shared_ptr<XMLConfiguration> &ptr, boost::shared_ptr<ScriptModel> &scriptModel);

private:

	/** create one Factor model based on the type
	 * @param typeName string with the type of the factor
	 * @param globalIndex the index of this factor as variable in the global var list
	 * @param factorIndex the index in the factor list
	 * @param ptr, the XMl configuration object */
	FactorModelBase* makeFactorModel(
			const string& typeName ,
			const Variable* var,
			int factorIndex,
			boost::shared_ptr<ScriptModel> &scriptModel,
			const boost::shared_ptr<XMLConfiguration> &ptr );

	/** number of factors */
	int nrFactors_;

	/** Vector holding the pointers to the */
	boost::ptr_vector<FactorModelBase> modelContainer_;

	/** the correlation matrix (the correlation matrix is model independent) <br>
	 * on the diagonal there are the volatility coefficients */
	DVector correlations_;

	/** Here we store the factors needed for forward estimation*/
	DVector standardDeviationFactors_;

	/** the name of the factors (as variables in the script) <br>
	 * Not all variables will appear in the script (e.g. underlying volatility factor)
	 * those have to be added to the global variable list*/
	std::vector<string> factorNames_;

	/** the vector contains the strings which specify the types of each factor*/
	std::vector<string> factorTypes_;

	/** The risk free interest rate, which is global <br>
	 * each model has this value internally, but this value is distributed to other entities */
	double r_;

	/** if the interest rate is coupled then return the variable, this is the factor index of that variable
	 * (but hte global index might change, for global index use intRateVariable_)
	 * else it is -1 -> interest rate is constant*/
	int r_index_;

	/** variable which is the interest rate*/
	const Variable* intRateVariable_;

	/** used for tracing the previous factor in the case of the Heston underlying volatility process */
	static int lastGlobalIndex_;

	/** used for tracing the previous factor in the case of the Heston underlying volatility process */
	static int lastFactorIndex_;

};
}
#endif /* FITOBMODELCOLLECTION_HPP_ */

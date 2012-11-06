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
 * FitobCalculator.hpp
 *
 *  Created on: Apr 15, 2010
 *      Author: benk
 */

#ifndef FITOBCALCULATOR_HPP_
#define FITOBCALCULATOR_HPP_

#include "src/utils/fitobdefs.hpp"

#include "src/evalcontext/FitobDomain.hpp"
#include "src/parser/FitobScriptParser.hpp"
#include "src/scriptmodel/FitobScriptModel.hpp"
#include "src/utils/FitobXMLConfiguration.hpp"
#include "src/diffusionmodel/FitobModelCollection.hpp"
#include "src/mesh/FitobGridFactory.hpp"
#include "src/mesh/FitobGridPlotter.hpp"
#include "src/mesh/FitobFullGrid.hpp"
#include "src/pdesolver/FitobSolverBase.hpp"
#include "src/montecarlo/FitobMCMachine.hpp"

namespace fitob{

using namespace std;

  /** Main class to do the complete calculation, price the product based on the
   * information from the script and from the configuration file. <br>
   * The calculation is mainly composed out of two steps: <br>
   * - forward calculation: in forward manner (as it is in the script specified)
   * we evaluate (some) operators. Here per each operator evaluation we store one
   * object so called Domain. The Domain not just represents the are which needs to be covered
   * but also should store the grading of each axis.<br>
   *  At the beginning of the calculations we must know the required maximal level
   *  so that each axis will have at least this resolution <br>
   *  <br>
   *  - backward calculation: based on the information stored in the forward calculation
   *  we create mesh (grid) based on the Domains stored, and then we only have to evaluate
   *  two operators, the Da operator (shift the grid) and the Theta operator (combi+solve PDE). */
  class FitobCalculator : public VerbClass{
  public:

	/** Ctor
	 * @param configurationXMLFile the XMl configuration file
	 * @param scriptFile the file which contains the script */
  FitobCalculator(const string& configurationXMLFile,
		          const string& scriptFile ,
		          const int inputLevel = -1,
		          const int convergenceRunNr = -1);

  virtual ~FitobCalculator() { ; }

  /** returns the price of this financial product */
  double evaluateScript();

  /** returns the price as a grid of this financial product */
  boost::shared_ptr<FullGrid> evaluateScript_FG(double& retVal, int level_in = 4);

  /** returns the price of this financial product calculated with Monte-Carlo*/
  DVector evaluateScript_MonteCarlo();

  /** does the forward estimation
   * todo: should be private later*/
  void forwardEstimation();

  /** does the backward calculation
   * todo: should be private later*/
  void backwardCalculation();

  /** returns reference to the XML configuration */
  inline const boost::shared_ptr<XMLConfiguration>& getXMLConfiguration() const { return xmlconfiguration_; }

  /** returns reference to the ScriptModel */
  inline const boost::shared_ptr<ScriptModel>& getScriptModel() const { return scriptmodel_; }

  /** returns reference to the ModelCollection */
  inline const boost::shared_ptr<ModelCollection>& getModelCollection() const { return modelcollection_; }

  /** returns the grid factory object*/
  inline const boost::shared_ptr<GridFactory>& getGridFactory() const {return gridFactory_;}

  /** returns the PDE solver */
  inline const boost::shared_ptr<SolverBase>& getSolver() const {return pdesolver_;}

  /** returns the grid plotter */
  inline const boost::shared_ptr<GridPlotter> getPlotter() const { return gridPlotter_;}

  /** returns the point at which the price can be evaluated*/
  inline const Domain* getEvalDomain() const {return evalDomain_.get();}

  /** returns the starting domain */
  inline const Domain* getStartDomain() const { return startDomain_.get(); }

  /** returns the domain which is needed for the norm evaluation */
  inline const Domain* getNormEvalDomain() const { return normMeasureDomain_.get(); }

  /** function below are for testing purposes (testing) result of forward calculation--- */
  inline const MeshContext& getLastContextFromStack() const { return contextStack_[contextStack_.size()-1];}

  /** returns the flag weather this */
  inline const bool dimensionAdaptiveComp() const { return dimensionAdaptive_; }

  /** returns the level vector in case of dimension adaptive case */
  inline const IVector& dimAdaptiveLevels() const { return dimensionAdaptiveLevels_; }

  private:

    /** script parser */
    boost::shared_ptr<ScriptParser> scriptparser_;

    /** XMl configuration which should be global */
    boost::shared_ptr<XMLConfiguration> xmlconfiguration_;

    /** model of the script */
    boost::shared_ptr<ScriptModel> scriptmodel_;

    /** model collection of the factors which influence the product */
    boost::shared_ptr<ModelCollection> modelcollection_;

    /** evaluation domain, which is a point in the N dimensional space*/
    boost::shared_ptr<Domain> evalDomain_;

    /** the starting domain */
    boost::shared_ptr<Domain> startDomain_;

    /** the norm measuring domain */
    boost::shared_ptr<Domain> normMeasureDomain_;

    /** context (Domain) stack */
    boost::ptr_vector<MeshContext> contextStack_;

    /** object to create mesh for the next domain (from the contextStack_)*/
    boost::shared_ptr<GridFactory> gridFactory_;

    /** the PDE solver for the theta operator*/
    boost::shared_ptr<SolverBase> pdesolver_;

    /** grid plotter for debugging purposes*/
    boost::shared_ptr<GridPlotter> gridPlotter_;

    /** Monte-Carlo simulation Engine */
    boost::shared_ptr<MCMachine> mCMachine_;

    /** time horizon for the whole script */
    double endTime_;

    /** how many domains we have on the stack */
	int stackIndex_;

	/** weather we have different levels in different dimensions */
	bool dimensionAdaptive_;

	/** in case of dimension adaptive computation we store the levels in each dimension*/
	IVector dimensionAdaptiveLevels_;

	/** the level for norm measuring */
	static int norm_level_measure;
};
}
#endif /* FITOBCALCULATOR_HPP_ */

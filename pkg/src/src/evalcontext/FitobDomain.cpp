/*
 * FitobDomain.cpp
 *
 *  Created on: Feb 26, 2010
 *      Author: benk
 */

#include <boost/lexical_cast.hpp>

#include "FitobDomain.hpp"

using namespace fitob;
using namespace std;


Domain::Domain(
		const IVector &axisLevel ,
		const std::vector<DVector> &axisGrading ,
		int maxRequiredLevel,
		int nrExportVariables
		) {
	nrExportVariable_ = nrExportVariables;
    setVerb(0);
	unsigned int nrAxis = axisLevel.size();
	int nrRealAxis = 0;

	maxRequiredLevel_ = 0;
	axisGrading_.resize(nrAxis);
	axisLevel_.resize(nrAxis);
	global_to_local_IndexMapping_.resize(nrAxis);
	local_to_global_IndexMapping_.resize(nrAxis);
	axisAverageValue_.resize(nrAxis+nrExportVariable_+1);
	maxGradAxisSize_ = 0;
	for (unsigned int i = 0 ; i < nrAxis ; i++)
	{
		axisLevel_[i] = axisLevel[i];
		FITOB_OUT_LEVEL3(verb(),"Domain Ctor ,i:" << i << " ,axisLevel_[i]:" << axisLevel_[i] );
		unsigned int gradSize = axisGrading[i].size()-1;
		// we test if the difference of the end of the axis is different from zero
		if (fabs(axisGrading[i][gradSize] - axisGrading[i][0]) > FITOB_NUMERICALZERO )
		{ // axis is not constant
			global_to_local_IndexMapping_[i] = nrRealAxis;
			local_to_global_IndexMapping_[nrRealAxis] = i;
			axisGrading_[i].resize(axisGrading[i].size());
			double avrg = 0.0;
			for (unsigned int jj = 0 ; jj < axisGrading[i].size() ; jj++ ){
				axisGrading_[i][jj] = axisGrading[i][jj];
				avrg += axisGrading[i][jj];
			}
			axisAverageValue_[i+(nrExportVariable_+1)] = avrg / (double)axisGrading[i].size();
			maxRequiredLevel_ = (axisLevel[i] > maxRequiredLevel_ ) ? axisLevel[i] : maxRequiredLevel_;
			nrRealAxis++;
		}
		else
		{ // axis is constant
			axisLevel_[i] = 0;
			global_to_local_IndexMapping_[i] = nrRealAxis;
			axisGrading_[i].resize(axisGrading[i].size());
			double avrg = 0.0;
			for (unsigned int jj = 0 ; jj < axisGrading[i].size() ; jj++ ){
				axisGrading_[i][jj] = axisGrading[i][jj];
				avrg += axisGrading[i][jj];
			}
			axisAverageValue_[i+(nrExportVariable_+1)] = avrg / (double)axisGrading[i].size();
		}
		maxGradAxisSize_ = (maxGradAxisSize_ < axisGrading_[i].size()) ? axisGrading_[i].size() : maxGradAxisSize_;
	}
	local_to_global_IndexMapping_.resize(nrRealAxis);
	nrImportVariable_ = nrAxis;
	nrAxis_ = nrRealAxis;

	// calculate the evaluation point
	evalPoint_.resize(this->nrRealAxis());
	for (int ind = 0 ; ind < this->nrRealAxis() ; ind++ ){
		// get the global coordinate for one local coordinate
		int globCoord = this->localToGlobalIndex(ind);
		const DVector& scaling = getGradedAxis(globCoord);
		// the middle point
		evalPoint_[ind] = scaling[ (scaling.size()) / 2 ];
	}
}


Domain::Domain( const IVector &axisLevel ,
			    const DVector &minAxisValues ,
			    const DVector &maxAxisValues ,
			    int nrExportVariables ,
			    const DVector* evalPoint ) {
	nrExportVariable_ = nrExportVariables;
	setVerb(0);
	unsigned int nrAxis = axisLevel.size();
	int nrRealAxis = 0;

	axisGrading_.resize(nrAxis);
	axisLevel_.resize(nrAxis);
	maxRequiredLevel_ = 0;
	maxGradAxisSize_ = 0;
	axisAverageValue_.resize(nrAxis+(nrExportVariable_+1));
	global_to_local_IndexMapping_.resize(nrAxis);
	local_to_global_IndexMapping_.resize(nrAxis);

	for (unsigned int i = 0 ; i < nrAxis ; i++)
	{
		axisLevel_[i] = axisLevel[i];
		FITOB_OUT_LEVEL3(verb(),"Domain Ctor ,i:" << i << " ,axisLevel_[i]:" << axisLevel_[i] );
		if (fabs(maxAxisValues[i] - minAxisValues[i]) > FITOB_NUMERICALZERO )
		{ // axis is not constant
			global_to_local_IndexMapping_[i] = nrRealAxis;
			local_to_global_IndexMapping_[nrRealAxis] = i;
			axisGrading_[i].resize(2);
			// todo: now we do here linear division we could do something more centered
			//       this remains as a future improvement
			axisGrading_[i].resize(powerTwo[axisLevel[i]]+1 , 0.0);
			for (int jj = 0; jj <= powerTwo[axisLevel[i]] ; jj++){
			     axisGrading_[i][jj] = minAxisValues[i] + ((double)jj/(double)(powerTwo[axisLevel[i]]))*(maxAxisValues[i]-minAxisValues[i]);
			}
			axisAverageValue_[i+(nrExportVariable_+1)] = (maxAxisValues[i] + minAxisValues[i]) / 2.0;
			maxRequiredLevel_ = (axisLevel[i] > maxRequiredLevel_ ) ? axisLevel[i] : maxRequiredLevel_;
			nrRealAxis++;
		}
		else
		{ // axis is constant
			axisLevel_[i] = 0;
			global_to_local_IndexMapping_[i] = nrRealAxis;
			axisGrading_[i].resize(1);
			axisGrading_[i][0] = minAxisValues[i];
			axisAverageValue_[i+(nrExportVariable_+1)] = minAxisValues[i];
		}
		maxGradAxisSize_ = (maxGradAxisSize_ < axisGrading_[i].size()) ? axisGrading_[i].size() : maxGradAxisSize_;
	}
	local_to_global_IndexMapping_.resize(nrRealAxis);
	nrImportVariable_ = nrAxis;
	nrAxis_ = nrRealAxis;

	// calculate the evaluation point
	if (evalPoint == NULL){
		evalPoint_.resize(this->nrRealAxis());
		for (int ind = 0 ; ind < this->nrRealAxis() ; ind++ ){
			// get the global coordinate for one local coordinate
			int globCoord = this->localToGlobalIndex(ind);
			// the middle point
			evalPoint_[ind] = (minAxisValues[globCoord-(nrExportVariable_+1)]+maxAxisValues[globCoord-(nrExportVariable_+1)])/2;
		}
	}
	else
	{
		evalPoint_ = (*evalPoint);
	}
}


Domain::Domain( const Domain& dom ){
	// here we just copy one object
	setVerb(0);
	FITOB_OUT_LEVEL3(verb(),"Domain copy Ctor 1 ");
	maxRequiredLevel_ = dom.maxRequiredLevel_;
	nrImportVariable_ = dom.nrImportVariable_;
	maxGradAxisSize_ = dom.maxGradAxisSize_;
	nrExportVariable_ = dom.nrExportVariable_;
	nrAxis_ = dom.nrAxis_;
	evalPoint_ = dom.evalPoint_;

	axisGrading_.resize(dom.axisGrading_.size());
	for (unsigned int jj = 0 ; jj < axisGrading_.size() ; jj++)
	{
		axisGrading_[jj].resize(dom.axisGrading_[jj].size());
	    for (unsigned int pp = 0 ; pp < dom.axisGrading_[jj].size() ; pp++)
		     axisGrading_[jj][pp] = dom.axisGrading_[jj][pp];
	}

	axisLevel_.resize( dom.axisLevel_.size());
	for (unsigned int jj = 0 ; jj < dom.axisLevel_.size() ; jj++)
		axisLevel_[jj] = dom.axisLevel_[jj];
	axisAverageValue_.resize( dom.axisAverageValue_.size() );
	for (unsigned int jj = 0 ; jj < dom.axisAverageValue_.size() ; jj++)
		axisAverageValue_[jj] = dom.axisAverageValue_[jj];
	local_to_global_IndexMapping_.resize( dom.local_to_global_IndexMapping_.size() );
	for (unsigned int jj = 0 ; jj < dom.local_to_global_IndexMapping_.size() ; jj++)
		local_to_global_IndexMapping_[jj] = dom.local_to_global_IndexMapping_[jj];
	global_to_local_IndexMapping_.resize( dom.global_to_local_IndexMapping_.size() );
	for (unsigned int jj = 0 ; jj < dom.global_to_local_IndexMapping_.size() ; jj++)
		global_to_local_IndexMapping_[jj] = dom.global_to_local_IndexMapping_[jj];
}


Domain::Domain( const Domain* dom ){
	setVerb(0);
	// here we just copy one object
	FITOB_OUT_LEVEL3(verb(),"Domain copy Ctor 2 ");
	maxRequiredLevel_ = dom->maxRequiredLevel_;
	nrImportVariable_ = dom->nrImportVariable_;
	maxGradAxisSize_ = dom->maxGradAxisSize_;
	nrExportVariable_ = dom->nrExportVariable_;
	nrAxis_ = dom->nrAxis_;
	evalPoint_ = dom->evalPoint_;

	axisGrading_.resize(dom->axisGrading_.size());
	for (unsigned int jj = 0 ; jj < axisGrading_.size() ; jj++)
	{
		axisGrading_[jj].resize(dom->axisGrading_[jj].size());
	    for (unsigned int pp = 0 ; pp < dom->axisGrading_[jj].size() ; pp++)
		     axisGrading_[jj][pp] = dom->axisGrading_[jj][pp];
	}

	axisLevel_.resize( dom->axisLevel_.size());
	for (unsigned int jj = 0 ; jj < dom->axisLevel_.size() ; jj++)
		axisLevel_[jj] = dom->axisLevel_[jj];
	axisAverageValue_.resize( dom->axisAverageValue_.size() );
	for (unsigned int jj = 0 ; jj < dom->axisAverageValue_.size() ; jj++)
		axisAverageValue_[jj] = dom->axisAverageValue_[jj];
	local_to_global_IndexMapping_.resize( dom->local_to_global_IndexMapping_.size() );
	for (unsigned int jj = 0 ; jj < dom->local_to_global_IndexMapping_.size() ; jj++)
		local_to_global_IndexMapping_[jj] = dom->local_to_global_IndexMapping_[jj];
	global_to_local_IndexMapping_.resize( dom->global_to_local_IndexMapping_.size() );
	for (unsigned int jj = 0 ; jj < dom->global_to_local_IndexMapping_.size() ; jj++)
		global_to_local_IndexMapping_[jj] = dom->global_to_local_IndexMapping_[jj];
}


void Domain::init()
{
	// we reinitialize the internal data structures
	// after a potential shift operation
	unsigned int nrAxis = axisLevel_.size();
	int nrRealAxis = 0;

	global_to_local_IndexMapping_.resize(nrAxis);
	local_to_global_IndexMapping_.resize(nrAxis);
	maxGradAxisSize_ = 0;
	for (unsigned int i = 0 ; i < nrAxis ; i++)
	{
		unsigned int gradSize = axisGrading_[i].size()-1;
		// we test if the difference of the end of the axis is different from zero
		if (fabs(axisGrading_[i][gradSize] - axisGrading_[i][0]) > FITOB_NUMERICALZERO )
		{ // axis is not constant
			global_to_local_IndexMapping_[i] = nrRealAxis;
			local_to_global_IndexMapping_[nrRealAxis] = i;
			nrRealAxis++;
		}
		else
		{ // axis is constant
			axisLevel_[i] = 0;
			global_to_local_IndexMapping_[i] = nrRealAxis;
		}
		maxGradAxisSize_ = (maxGradAxisSize_ < axisGrading_[i].size()) ? axisGrading_[i].size() : maxGradAxisSize_;
	}
	local_to_global_IndexMapping_.resize(nrRealAxis);
	nrImportVariable_ = nrAxis;
	nrAxis_ = nrRealAxis;

	// calculate the evaluation point
	evalPoint_.resize(this->nrRealAxis());
	for (int ind = 0 ; ind < this->nrRealAxis() ; ind++ ){
		// get the global coordinate for one local coordinate
		int globCoord = this->localToGlobalIndex(ind);
		const DVector& scaling = getGradedAxis(globCoord);
		// the middle point
		evalPoint_[ind] = scaling[ (scaling.size()) / 2 ];
	}
}


Domain* Domain::getNonScaledDomain() const {
	DVector axisMin(axisGrading_.size());
	DVector axisMax(axisGrading_.size());

    // loop over global indexes
	for (unsigned int i = 0 ; i < axisGrading_.size() ; i++){
		axisMin[i] = getGradedAxisMin(i+(nrExportVariable_+1));
		axisMax[i] = getGradedAxisMax(i+(nrExportVariable_+1));
	}

	// create the non scaled domain
	return ( new Domain( axisLevel_ , axisMin , axisMax , nrExportVariable_ , &(evalPoint_) ) );
}


void Domain::setAxis(DVector& newAxis , int axisGlobalIndex ) {
	// set the new grading
	axisGrading_[axisGlobalIndex - (nrExportVariable_+1) ] = newAxis;
	// calculate average
	double avrg = 0.0;
	for (unsigned int jj = 0 ; jj < newAxis.size() ; jj++ ){
		avrg += newAxis[jj];
	}
	axisAverageValue_[axisGlobalIndex] =  avrg / (double) newAxis.size();
	axisLevel_[axisGlobalIndex-(nrExportVariable_+1)] = (newAxis.size() > 2) ? this->getMaximalLevel() : 0 ;
	// initialize internal data structures
	init();
}


void Domain::applyDimAdaptLevels(const IVector& adaptLevels){

	// this method applies dimension adaptivity to the domain

   int adaptvectSize = adaptLevels.size();

   for (int jj = 0 ; jj < nrAxis_ ; jj++){
	   if (jj >= adaptvectSize ) break;
       if ( axisLevel_[jj] != adaptLevels[jj]){

    	   int globalIndex = local_to_global_IndexMapping_[jj];
    	   // the index in the vector can be only positive
    	   int offs = fitob::powerTwo[axisLevel_[globalIndex] - adaptLevels[jj]];

           DVector newGrading(fitob::powerTwo[adaptLevels[jj]]+1 , 0.0);
           // copy the grading vector (we size down the grading vector)
           double avrg = 0.0;
           for (unsigned int p = 0 ; p < newGrading.size() ; p++)
           {
        	   newGrading[p] = axisGrading_[globalIndex][p*offs];
        	   avrg += newGrading[p];
        	   //FITOB_OUT_LEVEL3( 5 , "[" << p << "]" << newGrading[p] << ",");
           }
           //FITOB_OUT_LEVEL3( 5 , "Domain::applyDimAdaptLevels gi:" << globalIndex << " , "<< axisLevel_[globalIndex] << " to " << adaptLevels[jj]);
           axisGrading_[globalIndex] = newGrading;
           axisAverageValue_[globalIndex+(nrExportVariable_+1)] = avrg/(double)newGrading.size();
    	   axisLevel_[globalIndex] = adaptLevels[jj];
       }
   }
   // initialize internal data structures
   init();
}


const string Domain::toString() const{

	string retString;

	retString = retString + " Domain print , maxRequiredLevel:" + boost::lexical_cast<std::string>(maxRequiredLevel_) + "\n";
	for (unsigned int i = 0 ; i < axisLevel_.size() ; i++)
	{
		retString = retString + " level:" + boost::lexical_cast<std::string>(axisLevel_[i]) + " Grading:" + "\n";
		for (unsigned int ii = 0 ; ii < axisGrading_[i].size() ; ii++)
			retString = retString + boost::lexical_cast<std::string>(axisGrading_[i][ii]) + ",";
		retString = retString + " \n";
	}
	retString = retString + " Local to global:";
	for (unsigned int ii = 0 ; ii < local_to_global_IndexMapping_.size() ; ii++)
		retString = retString + boost::lexical_cast<std::string>(local_to_global_IndexMapping_[ii]) + ",";
	retString = retString + " \n";
	retString = retString + " Global to local:";
	for (unsigned int ii = 0 ; ii < global_to_local_IndexMapping_.size() ; ii++)
		retString = retString + boost::lexical_cast<std::string>(global_to_local_IndexMapping_[ii]) + ",";
	retString = retString + " \n";
	return retString;
}

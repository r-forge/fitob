/*
 * FitobFullGrid_WB.cpp
 *
 *  Created on: Jan 1, 2011
 *      Author: benk
 */

#include "FitobFullGrid_WB.hpp"
#include "src/operators/FitobOperatorSequence.hpp"
#include <boost/lexical_cast.hpp>

using namespace fitob;
using namespace std;

FullGrid_WB::FullGrid_WB(const Domain* dom) :
   FullGridBase( dom , "FullGrid_WB" ) {

	setVerb(0);

	actualLevels_.resize(domain()->nrRealAxis());
	// here we set the actual levels of the full grid from the domain
	FITOB_OUT_LEVEL1(verb(),"FullGrid_WB query levels , real axis:" << domain()->nrRealAxis() );
	FITOB_OUT_LEVEL3(verb()," FullGrid_WB::init , domain:" << domain()->toString() );
	for (int ii = 0 ; ii < domain()->nrRealAxis() ; ii++){
		FITOB_OUT_LEVEL3(verb(),"FullGrid_WB get levels , ii: " << ii << " , globIndex:" << domain()->localToGlobalIndex(ii));
		FITOB_OUT_LEVEL3(verb(),"FullGrid_WB" << ii << " , level:" << domain()->getAxisLevel( domain()->localToGlobalIndex(ii)) );
		actualLevels_[ii] = domain()->getAxisLevel( domain()->localToGlobalIndex(ii));
	}
	// call init
	init();

	FITOB_OUT_LEVEL3(verb(),"FullGrid_WB created");
}

FullGrid_WB::FullGrid_WB(const Domain* dom, IVector& specialLevels) :
  FullGridBase( dom , "FullGrid_WB" ) {
	actualLevels_ = specialLevels;
	// call init
	init();
}

FullGrid_WB::FullGrid_WB(const Domain* dom , IVector& specialLevels , DVector* values) :
  FullGridBase( dom , "FullGrid_WB" ) {
	actualLevels_ = specialLevels;
	// take the input vector as our vector
	values_ = values; values_created_ = false;
	// call init
	init();
}

void FullGrid_WB::init(){

	setVerb(0);

	nrAxis_ = domain()->nrRealAxis();
	nrPointsPerAxis_.resize(domain()->nrRealAxis());
	axisOffset_.resize(domain()->nrRealAxis());
	axisGrading_.resize(domain()->nrRealAxis());

	FITOB_OUT_LEVEL1(verb()," FullGrid_WB::init START");
	FITOB_OUT_LEVEL2(verb()," FullGrid_WB::init Domain:" << domain()->toString());
	// calculate offsets
	int offset = 1;
	for (int ii = 0 ; ii < domain()->nrRealAxis() ; ii++){
		axisOffset_[ii] = offset;
		nrPointsPerAxis_[ii] = fitob::powerTwo[actualLevels_[ii]] - 1;
		FITOB_OUT_LEVEL2(verb()," FullGrid_WB::init , ii: " << ii << " , axisOffset_[ii]" << axisOffset_[ii] );
		FITOB_OUT_LEVEL2(verb()," FullGrid_WB::init , ii: " << ii << " , actualLevels_[ii]:" << actualLevels_[ii] );
		FITOB_OUT_LEVEL2(verb()," FullGrid_WB::init , ii: " << ii << " , offset:" << offset << " , nrPointsPerAxis_[ii]:" << nrPointsPerAxis_[ii]);
		offset = offset * nrPointsPerAxis_[ii];
	}
	// this must be the total length of the vector
	nrTotalPoints_ = offset;
	FITOB_OUT_LEVEL2(verb()," FullGrid_WB::init , nrTotalPoints_:" << nrTotalPoints_);

	// allocate the vector for all the points, and set to zero
	//values_ = new double[nrTotalPoints_];
	if (values_ == 0){
		values_ = new DVector();
		(*values_).resize(nrTotalPoints_, 0.0);
	} else {
	    FITOB_ERROR_TEST( ( (nrTotalPoints_ == (int)(*values_).size()) || ((*values_).size() == 0) ) ,
	    		"If Input vector was given size must match , values.size(): " << (*values_).size() << " nrTotalPoints_" << nrTotalPoints_ );
	}

	FITOB_OUT_LEVEL3(verb()," FullGrid_WB::init , copy axis scaling domain()->nrRealAxis(): " << domain()->nrRealAxis() );
	// now create the grading of the
	IVector tmp_mult( domain()->nrRealAxis() );
	// calculate the multiplication factor
	for (int ii = 0 ; ii < domain()->nrRealAxis() ; ii++){
		int index = (domain()->getAxisLevel((domain()->localToGlobalIndex(ii))) - actualLevels_[ii]);
		FITOB_ERROR_TEST( index >= 0 , " FullGrid_WB::init index >= 0 , the proposed level is probably higher then in the Domain");
		tmp_mult[ii] = fitob::powerTwo[ index ] ;
		axisGrading_[ii].resize(nrPointsPerAxis_[ii]+2); // we add two since here we do not have boundary points
		FITOB_OUT_LEVEL3(verb()," FullGrid_WB::init , ii: " << ii << " , domain()->getAxisLevel( (domain()->localToGlobalIndex(ii)):"
				<< domain()->getAxisLevel( (domain()->localToGlobalIndex(ii)) )  << " , actualLevels_[ii]:" << actualLevels_[ii]);
		FITOB_OUT_LEVEL3(verb()," FullGrid_WB::init , ii: " << ii << " , tmp_mult[ii]:"
				<< tmp_mult[ii] << ", actualLevels_[ii]:" << actualLevels_[ii] << ", index:" << index);
	}

	// copy the axis grading
	for (int ii = 0 ; ii < domain()->nrRealAxis() ; ii++){
		 //if (verb() > 3) std::cout << "Axis scaling ii:" << ii << " | ";
         for (int jj = 0 ; jj < nrPointsPerAxis_[ii]+2 ; jj++){
        	 axisGrading_[ii][jj] =
        		domain()->getGradedAxis(domain()->localToGlobalIndex(ii))[ jj * tmp_mult[ii]];
        	 //if (verb() > 3) std::cout << axisGrading_[ii][jj] << " , ";
         }
         //if (verb() > 3) std::cout << std::endl;
	}
	FITOB_OUT_LEVEL1(verb()," FullGrid_WB::init END");
}

FullGrid_WB::~FullGrid_WB() {
   // the only task is to delete the dynamic allocated array
   //if (values_ != NULL){
	    //delete[] values_;
   //}
}

double FullGrid_WB::eval(const DVector& globalCoords) const
{
	DVector localCoords(nrAxis_);
	DVector intersect(nrAxis_);
	IVector minIndex(nrAxis_);
	IVector maxIndex(nrAxis_);
	IVector localCoord(nrAxis_);
	std::vector<IVector*> aIndex(2);
	std::vector<DVector> intersectD2(2);

	intersectD2[0].resize(nrAxis_); intersectD2[1].resize(nrAxis_);

	int middle = 0;
	domain()->globalToLocal( globalCoords , localCoords);

	//FITOB_OUT_LEVEL3(verb()," FullGrid_WB::eval START ");
	// find the position on each axis
	for (int ii = 0 ; ii < nrAxis_ ; ii++){
		minIndex[ii] = 0;
		maxIndex[ii] = nrPointsPerAxis_[ii]+1;

		// this should have a complexity of O(logN) or O(L), L-> level
		// bisection search, since we assume that the array is ordered
		//FITOB_OUT_LEVEL3(verb()," FullGrid_WB::eval , globalCoords[globalI]:" << globalCoords[domain()->localToGlobalIndex(ii)]);
		//FITOB_OUT_LEVEL3(verb()," FullGrid_WB::eval , localCoords[ii]:" << localCoords[ii]);
		for (; ;){
		  middle = (minIndex[ii] + maxIndex[ii])/2;
		  //FITOB_OUT_LEVEL3(verb()," FullGrid_WB::eval , ii: " << ii << ", middle:" << middle <<
				  //", axisGrading_[ii][middle]:" << axisGrading_[ii][middle]);
		  if (localCoords[ii] >= axisGrading_[ii][middle] ){
			  minIndex[ii] = middle;
			  //FITOB_OUT_LEVEL3(verb()," minIndex[ii] = " << minIndex[ii]);
		  }
		  else{
			  maxIndex[ii] = middle;
			  //FITOB_OUT_LEVEL3(verb()," maxIndex[ii] = " << maxIndex[ii]);
		  }
          // break the for cycle when we are close enough
		  if ((maxIndex[ii] - minIndex[ii]) < 2) break;
		}

		// calculate the linear
		intersect[ii] = (axisGrading_[ii][maxIndex[ii]] - localCoords[ii]) /
				(axisGrading_[ii][maxIndex[ii]] - axisGrading_[ii][minIndex[ii]]);
		//FITOB_OUT_LEVEL3(verb()," FullGrid_WB::eval , ii: " << ii << ", localCoords[ii]:" << localCoords[ii]);
		//FITOB_OUT_LEVEL3(verb()," FullGrid_WB::eval , ii: " << ii << ", intersect[ii]:" << intersect[ii]);
		//FITOB_OUT_LEVEL3(verb()," FullGrid_WB::eval , ii: " << ii << ", axisGrading_[ii][minIndex[ii]]:" << axisGrading_[ii][minIndex[ii]]);
		//FITOB_OUT_LEVEL3(verb()," FullGrid_WB::eval , ii: " << ii << ", axisGrading_[ii][maxIndex[ii]]:" << axisGrading_[ii][maxIndex[ii]]);
		//FITOB_OUT_LEVEL3(verb()," FullGrid_WB::eval , ii: " << ii << ", minIndex[ii]:" << minIndex[ii]);
		//FITOB_OUT_LEVEL3(verb()," FullGrid_WB::eval , ii: " << ii << ", maxIndex[ii]:" << maxIndex[ii]);
	}

	// storing some intermediate results
	aIndex[0] = &(minIndex);
	aIndex[1] = &(maxIndex);
	for (int ii = 0 ; ii < nrAxis_ ; ii++){
		// test if we have to extrapolate left
		if (minIndex[ii] <= 0){
			(*aIndex[0])[ii] = 0;
			(*aIndex[1])[ii] = 1;
			double frac = (axisGrading_[ii][1] - axisGrading_[ii][0]) / (axisGrading_[ii][2] - axisGrading_[ii][1]);
			intersect[ii] = frac*((1-intersect[ii]) - 1);
			intersectD2[0][ii] = 1 - intersect[ii] ;
			intersectD2[1][ii] = intersect[ii];
			//FITOB_OUT_LEVEL3(verb()," FullGrid_WB::eval LEFT, ii: " << ii << ", intersectD2[0][ii]:" << intersectD2[0][ii]);
			//FITOB_OUT_LEVEL3(verb()," FullGrid_WB::eval LEFT, ii: " << ii << ", intersectD2[1][ii]:" << intersectD2[1][ii]);
			continue;
		}
		// test if we have to extrapolate right
		if (maxIndex[ii] >= nrPointsPerAxis_[ii]+1 ){
			(*aIndex[0])[ii] = nrPointsPerAxis_[ii] - 2;
			(*aIndex[1])[ii] = nrPointsPerAxis_[ii] - 1;
			int nrPoints = nrPointsPerAxis_[ii]+1;
			double frac = (axisGrading_[ii][nrPoints] - axisGrading_[ii][nrPoints-1]) / (axisGrading_[ii][nrPoints-1] - axisGrading_[ii][nrPoints-2]);
			intersect[ii] = 1 + (1-intersect[ii])*frac;
			intersectD2[0][ii] = 1 - intersect[ii];
			intersectD2[1][ii] = intersect[ii];
			//FITOB_OUT_LEVEL3(verb()," FullGrid_WB::eval RIGHT, ii: " << ii << ", intersectD2[0][ii]:" << intersectD2[0][ii]);
			//FITOB_OUT_LEVEL3(verb()," FullGrid_WB::eval RIGHT, ii: " << ii << ", intersectD2[1][ii]:" << intersectD2[1][ii]);
			continue;
		}
		minIndex[ii] = minIndex[ii]-1;
		maxIndex[ii] = maxIndex[ii]-1;
		intersectD2[0][ii] = intersect[ii];
		intersectD2[1][ii] = 1-intersect[ii];
	}

	// once we have the indexes now do the interpolation with tensor product
    // N-linear interpolation or extrapolation when it is outside the domain
	double val = 0.0;
	for (int ii=0; ii < fitob::powerTwo[nrAxis_] ; ii++){
		 double baseVal = 1;
		 int i = 0;
		 int tmp_val = ii;
		 int vv = 0;
		 // here we use a more efficient version , without branches
		 // make N-liniear interpolation , or extrapolation
         for (int jj = 0 ; jj < nrAxis_ ; jj++){
        	 vv = (tmp_val & 1); // tmp % 2 ;
          	 baseVal = baseVal * intersectD2[vv][jj];
           	 localCoord[jj] = (*aIndex[vv])[jj];
           	 i += (*aIndex[vv])[jj]* axisOffset_[jj];
        	 tmp_val = tmp_val >> 1; //tmp_val / 2; // tmp_val >> 1;
           	 //FITOB_OUT_LEVEL3(verb()," FullGrid_WB::eval , ii:" << ii << " , jj:" << jj);
           	 //FITOB_OUT_LEVEL3(verb()," FullGrid_WB::eval , (*aIndex[vv])[ii]:" << (*aIndex[vv])[jj]);
           	 //FITOB_OUT_LEVEL3(verb()," FullGrid_WB::eval , vv:" << vv << " , tmp_val:" << tmp_val << ", i:" << i );
         }
         // multiply the basis function value with the coefficient
         //FITOB_ERROR_TEST( i < nrTotalPoints_ , "FullGrid_WB::eval , out of bound");
         //FITOB_OUT_LEVEL3(verb(),"Point i:" << i << " values_[i]:" << values_[i] << ", baseVal:" << baseVal);
         val += baseVal* (*values_)[i];
	}
    // just return the value
	//FITOB_OUT_LEVEL3(verb()," FullGrid_WB::eval END: " << val);
	return val;
}

void FullGrid_WB::eval(const std::vector<DVector>& globalCoordonates , DVector& resVector) const
{
   for (unsigned int ii = 0 ; ii < globalCoordonates.size() ; ii++){
	   // we just call the simple version of this function
	   // todo: this could be optimized by copy paste the body of the function above
	   resVector[ii] = this->eval(globalCoordonates[ii]);
   }
}

void FullGrid_WB::setValues(const Evaluable* func , const DVector& globalCoords)
{
#if defined(FITOB_OPENMP)
#pragma omp parallel
{
#endif
	IVector axisInd(nrAxis_);
	DVector localCoor(nrAxis_);
	DVector tmpGlCoo  = globalCoords;
	int tmp_vect;
#if defined(FITOB_OPENMP)
#pragma omp for schedule(static)
#endif
	for (int pind = 0 ; pind < nrTotalPoints_ ; pind++){
		// get the index
		tmp_vect = pind;
		for (int ii = nrAxis_-1 ; ii >= 0 ; ii--){
			axisInd[ii] = (tmp_vect / axisOffset_[ii]);
			tmp_vect = (tmp_vect % axisOffset_[ii]);
			localCoor[ii] = axisGrading_[ii][ 1 + axisInd[ii] ];
			FITOB_OUT_LEVEL3(verb(),"FullGrid_WB::setValues ii:" << ii <<
					", localCoor[ii]:" << localCoor[ii] << ", "<<axisInd[ii] << ", axisOffset_[ii]:" << axisOffset_[ii]);
		}
        // transform to global coordinates
        domain()->localToGlobal(localCoor,tmpGlCoo);
        // eval the function
        // TODO: we could do this for work sets, since this initiates the evaluation process
        tmpGlCoo[1] = (*values_)[pind];
        (*values_)[pind] = func->eval(tmpGlCoo);
        FITOB_OUT_LEVEL3(verb(),"FullGrid_WB::setValues pind:" << pind << " , values_[pind]:" << (*values_)[pind]);
	}
#if defined(FITOB_OPENMP)
}
#endif
}

void FullGrid_WB::applyConstraints(const OperatorSequence* constraintOpSeq_ ,
		                        const DVector& globalCoords )
{
#if defined(FITOB_OPENMP)
#pragma omp parallel
{
#endif
	IVector axisInd(nrAxis_);
	DVector localCoor(nrAxis_);
	DVector tmpGlCoo  = globalCoords;
	int tmp_vect;
#if defined(FITOB_OPENMP)
#pragma omp for schedule(static)
#endif
	for (int pind = 0 ; pind < nrTotalPoints_ ; pind++){
		// get the index
		tmp_vect = pind;
		for (int ii = nrAxis_-1 ; ii >= 0 ; ii--){
			axisInd[ii] = (tmp_vect / axisOffset_[ii]);
			tmp_vect = (tmp_vect % axisOffset_[ii]);
			localCoor[ii] = axisGrading_[ii][ 1 + axisInd[ii] ];
		}
        // transform to global coordinates
        domain()->localToGlobal(localCoor,tmpGlCoo);
        tmpGlCoo[1] = (*values_)[pind];
        // eval the function
        // here call the restriction for one point
        // todo: it would be more efficient if we would call the method for one couple of global coordinates
        constraintOpSeq_->applyConstraints(tmpGlCoo);
        // the result should be in the first vector element
        (*values_)[pind] = tmpGlCoo[1];
	}
#if defined(FITOB_OPENMP)
}
#endif
}


string FullGrid_WB::toStringFG() const {
	string ret = " FullGrid_WB , dim: " + boost::lexical_cast<std::string>(actualLevels_.size());
	for ( unsigned int ii = 0 ; ii < actualLevels_.size() ; ii++){
		ret = ret + " | L:" +  boost::lexical_cast<std::string>(actualLevels_[ii]) +
				" , NrP:" +  boost::lexical_cast<std::string>(nrPointsPerAxis_[ii]) +
				" , OFs:" +  boost::lexical_cast<std::string>(axisOffset_[ii]);
	}
    return ret;
}

// --------- plotting for debugging purposes --------
void FullGrid_WB::plotGrid(const string& filename , int count)
{
	   const string completeName = filename + boost::lexical_cast<std::string>(count)+".m";
	   DVector globalCoord = domain()->getAverage();

	   if (axisGrading_.size() >= 2)
	   {
		  DVector result( axisGrading_[0].size() * axisGrading_[1].size() );
		  // loop and evaluate points
	      for ( unsigned int ii = 0 ; ii < axisGrading_[0].size() ; ii++){
	    	  for ( unsigned int jj = 0 ; jj < axisGrading_[1].size() ; jj++){
	    		  // todo: this will work only in 2D
	    		 globalCoord[2] = axisGrading_[0][ii];
	    		 globalCoord[3] = axisGrading_[1][jj];
	    		 // eval the grid at those positions
	    	     result[ii*axisGrading_[1].size() + jj] = this->eval(globalCoord);
	    	  }
	      }

	      ofstream myfile;
	      myfile.open(completeName.c_str());
	      myfile << "X = [ " << axisGrading_[0][0];
	      for ( unsigned int ii = 1 ; ii < axisGrading_[0].size() ; ii++){
	   	     myfile << " , " << axisGrading_[0][ii];
	      }
	      myfile << "]; \n ";
	      myfile << "Y = [ " << axisGrading_[1][0];
	      for ( unsigned int ii = 1 ; ii < axisGrading_[1].size() ; ii++){
	   	     myfile << " , " << axisGrading_[1][ii];
	      }
	      myfile << "]; \n ";
	      myfile << "res = [ " << result[0];
	      for ( unsigned int ii = 1 ; ii < axisGrading_[0].size()*axisGrading_[1].size() ; ii++) {
	    	  if ( (ii % axisGrading_[1].size()) == 0){
	    		  myfile << " ; " << result[ii];
	    	  }
	    	  else{
	   	         myfile << " , " << result[ii];
	    	  }
	      }
	      myfile << "]; \n ";
	      myfile << "[x,y]=meshgrid(Y,X);\n";
	      myfile << " surf(x,y,res); \n ";
	      myfile.close();
	   }
}

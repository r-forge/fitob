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
 * fitobdefs.hpp
 *
 *  Created on: Feb 26, 2010
 *      Author: benk
 */

#ifndef FITOBDEFS_HPP_
#define FITOBDEFS_HPP_

#include <boost/assert.hpp>
#include <boost/array.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <boost/utility.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include <boost/math/complex/atan.hpp>

#include <math.h>
#include <iostream>
#include <stack>
#include <functional>
#include <string>
#include <cassert>
#include <fstream>
#include <vector>

#include "src/utils/fitobOpenMP.hpp"
#include "src/utils/fitobMPI.hpp"

/** -------- Vector definitions ------------ */

typedef std::vector<double> DVector;

typedef std::vector<int> IVector;

/** ------- Assert definition ------ */
#define FITOB_ASSERT(val) BOOST_ASSERT(val);

#define NO_STD_ERROR

/** -------- */

namespace fitob{
   /** vector storing the powers of two */
   const int powerTwo[23] = { 1 , 2 , 4 , 8 , 16 , 32 , 64 , 128 , 256 , 512 , 1024 , 2048 , 4096 , 8192 , 16384 , 65536 , 131072 , 262144 , 524288, 1048576 ,
		                      2097152 , 4194304 , 8388608};

   /** class which has a verbosity attribute*/
   class VerbClass{
   public:
	   void setVerb(int v) { verb_ = v; }
	   inline short int verb() const { return verb_;}
   private:
	   short int verb_;
   };
}

//NumericalZero
#define    FITOB_NUMERICALZERO      1e-14

#define    FITOB_HALF_PI            1.57079632679490

#define    FITOB_PI                  3.14159265358979


#define isZero(value) bool(fabs(value)<1e-12);


#ifdef NO_STD_ERROR
/** Unconditionally write out */
#if defined(FITOB_MPI)
#define FITOB_OUT(str) {};
#else
#define FITOB_OUT(str) {};
#endif

#define FITOB_ERROR_TEST(test,str) {};

#define FITOB_ERROR_EXIT(str) {};

#define FITOB_WARNING_TEST(test,str) {};

#define FITOB_WARNING(str) {};

#define FITOB_ERROR_MSG(str) {};

#else
	/** Unconditionally write out */
	#if defined(FITOB_MPI)
	#define FITOB_OUT(str) { \
			std::cout << FITOB_MPI_Comm_rank() << ":" << str << std::endl; };
	#else
	#define FITOB_OUT(str) { \
			std::cout << str << std::endl; };
	#endif


	#define FITOB_ERROR_TEST(test,str) { \
	    if (!(test)) std::cout << std::endl << "ERROR: " << str << std::endl; \
	    BOOST_ASSERT(test); } \


	#define FITOB_ERROR_EXIT(str) { \
	    std::cout << std::endl << "ERROR: " << str << std::endl; \
	    BOOST_ASSERT(false); } \

	#define FITOB_WARNING_TEST(test,str) { \
	    if (!(test)) std::cout << std::endl << "WARNING: " << str << std::endl; \
	    } \


	#define FITOB_WARNING(str) { \
	    std::cout << std::endl << "WARNING: " << str << std::endl; \
	     } \

	#define FITOB_ERROR_MSG(str) { \
	    std::cout << std::endl << "ERROR: " << str << std::endl; } \

#endif

/** Conditionally write out */
#define FITOB_OUT_LEVEL(limit,level,str) \
	   if (limit < level ) FITOB_OUT(str);


/** Conditionally write out */
#define FITOB_OUT_LEVEL1(level,str) \
	   if (level > 1) FITOB_OUT(str);


#define FITOB_OUT_LEVEL2(level,str) \
	   if (level > 2) FITOB_OUT(str);


#define FITOB_OUT_LEVEL3(level,str) \
	   if (level > 3) FITOB_OUT(str);

#define FITOB_OUT_LEVEL4(level,str) \
	   if (level > 4) FITOB_OUT(str);

#define FITOB_OUT_LEVEL5(level,str) \
	   if (level > 5) FITOB_OUT(str);

#define FITOB_OUT_LEVEL6(level,str) \
	   if (level > 6) FITOB_OUT(str);

#define FITOB_OUT_LEVEL7(level,str) \
	   if (level > 7) FITOB_OUT(str);

// ---------- all the utility functions are in the fitob namespace -------------

namespace fitob{

  /** the maximum double value */
  inline double FITOB_DMAX(double v1 , double v2)   { return (v1<v2)?v2:v1; }

  /** the maximum int value */
  inline int FITOB_IMAX(int v1 , int v2)   { return (v1<v2)?v2:v1; }

  /** the minimum double value */
  inline double FITOB_DMIN(double v1 , double v2)   { return (v1>v2)?v2:v1; }

  /** the minimum int value */
  inline int FITOB_IMIN(int v1 , int v2)   { return (v1>v2)?v2:v1; }

  /** calculates the difference between two vectors and stores the result in the first one*/
  inline void vect_diff(DVector* v1 , const DVector* v2)
  {
	FITOB_ERROR_TEST( v1->size() == v2->size() , " L2Norm , size do not match v1->size():"
			<< v1->size() << " , v2->size():" << v2->size());
	for (unsigned int i = 0; i < v1->size() ; i++){
		v1->at(i) = v1->at(i) - v2->at(i);
	}
  }

  /** calculates the L2 norm compared to an other vector*/
  inline double l2_norm(const DVector* v1 , const DVector* v2)
  {
	FITOB_ERROR_TEST( v1->size() == v2->size() , " L2Norm , size do not match , v1->size():"
			<< v1->size() << " , v2->size():" << v2->size());
	double diff = 0.0 , tmp;
	//std::cout << " L2 norm , diff vect: ";
	for (unsigned int i = 0; i < v1->size() ; i++){
		tmp = v1->at(i) - v2->at(i);
		diff = diff + tmp*tmp;
		//std::cout << tmp << ",";
	}
	//std::cout << std::endl;
	// return norm
	return sqrt(diff/double(v1->size()));
  }

  /** calculates the L2 norm of one vector*/
  inline double l2_norm(const DVector* v1)
  {
	double diff = 0.0;
	for (unsigned int i = 0; i < v1->size() ; i++){
		diff = diff + v1->at(i)*v1->at(i);
	}
	return sqrt(diff/double(v1->size()));
  }

  /** calculates the Inf norm compared to an other vector*/
  inline double inf_norm(const DVector* v1 , const DVector* v2)
  {
	FITOB_ERROR_TEST( v1->size() == v2->size() , " InfNorm , size do not match v1->size():"
			<< v1->size() << " , v2->size():" << v2->size());
	double diff = 0.0 , tmp;
	//std::cout << " Inf norm , diff vect: ";
	for (unsigned int i = 0; i < v1->size() ; i++){
		tmp = fabs(v1->at(i) - v2->at(i) );
		diff = (tmp > diff) ? tmp : diff;
		//std::cout << tmp << ",";
	}
	//std::cout << std::endl;
	return diff;
  }

  /** calculates the Inf norm of one vector*/
  inline double inf_norm(const DVector* v1 )
  {
	double diff = 0.0 , tmp;
	for (unsigned int i = 0; i < v1->size() ; i++){
		tmp = fabs(v1->at(i));
		diff = (tmp > diff) ? tmp : diff;
	}
	return diff;
  }
}

#endif /* FITOBDEFS_HPP_ */

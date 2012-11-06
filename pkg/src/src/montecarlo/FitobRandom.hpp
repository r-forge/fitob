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
 * FitobRandom.hpp
 *
 *  Created on: Mar 15, 2011
 *      Author: benk
 */

#ifndef FITOBRANDOM_HPP_
#define FITOBRANDOM_HPP_

#include "src/utils/fitobdefs.hpp"
#include "src/utils/FitobXMLConfiguration.hpp"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>


namespace fitob {

/** Class to generate normaly distributed psehudo numbers */
class Random : public VerbClass{
public:

	/** Ctor for the random number generator */
	Random( const XMLConfiguration* xmlConfig );

	virtual ~Random() {;}

	/** sets the seed for the random generator */
	void setSeed(int seed);

	/** returns "nrRandVar" nr normally distributed random variables */
	void getRandomN( int nrRandVar , DVector& randOut ) ;

private:

	/** random engine */
	boost::mt19937 engine_;

	/** normal number distribution*/
	boost::normal_distribution<> normalDistribution_;

	boost::variate_generator< boost::mt19937& , boost::normal_distribution<> > variateGenrator_;

};

}

#endif /* FITOBRANDOM_HPP_ */

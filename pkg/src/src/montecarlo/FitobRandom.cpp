/*
 * FitobRandom.cpp
 *
 *  Created on: Mar 15, 2011
 *      Author: benk
 */

#include "FitobRandom.hpp"

using namespace fitob;
using namespace std;

Random::Random( const XMLConfiguration* xmlConfig ) : engine_() , normalDistribution_(0.0,1.0) ,
		variateGenrator_( engine_ , normalDistribution_ ){

	setVerb(6);

	// set the seed from the XML file
	if ( xmlConfig!= 0 ){
		int seed = xmlConfig->getIntConfiguration("thetaconfigurations.montecarlo.RANDOM_GENERATOR.<xmlattr>.seed");
		if (seed >= 0) engine_.seed(seed);
	}
}


void Random::setSeed(int seed) {
	// set the seed for the engine
	engine_.seed(seed);
}


void Random::getRandomN( int nrRandVar , DVector& randOut ) {
	// generate nrRandVar piece of random numbers
	for (int n = 0 ; n < nrRandVar ; n++ ){
		randOut[n] = variateGenrator_();
		//FITOB_OUT_LEVEL2(verb()," RandN : " << randOut[n] );
	}
}

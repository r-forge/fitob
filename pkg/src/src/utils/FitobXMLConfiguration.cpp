/*
 * FitobXMLConfiguration.cpp
 *
 *  Created on: Feb 26, 2010
 *      Author: benk
 */

#include "FitobXMLConfiguration.hpp"

#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <boost/lexical_cast.hpp>

using namespace fitob;
using namespace std;

XMLConfiguration::XMLConfiguration(const string& configFileName) {

  setVerb(0);
  // Load the XML file into the property tree. If reading fails
  // (cannot open file, parse error), an exception is thrown.
  read_xml(configFileName, pt_ , 0x2 ); // means no comment

}

// return the number of children
int XMLConfiguration::nrXMLNodes(const string& pathInTree) const{
	int counter = 0;
	// todo: improve this later
    BOOST_FOREACH(const ptree::value_type &v , pt_.get_child(pathInTree))
    {   FITOB_OUT_LEVEL3( verb(), v.first );
    	counter++;
    }
    return counter;
}

int XMLConfiguration::getIntConfiguration(
		                const string& pathInTree_first,
		                const int indexChild ,
		                const string& pathInTree_second) const{
	int counter = 0;
	BOOST_FOREACH(const ptree::value_type &v , pt_.get_child(pathInTree_first))
    {
    	if (indexChild <= counter)
    	  {
    		// take the tree , and ask for the second path inside the tree
    		const ptree p = v.second;
    		FITOB_OUT_LEVEL3(verb(), " FOUND SUBTREE" << v.first);
    		return p.get(pathInTree_second, -1);
    		break;
    	  }
    	else
    	{counter++;}
    }
    return -1; //Wall
}

double XMLConfiguration::getDoubleConfiguration(
		                       const string& pathInTree_first,
		                       const int indexChild ,
		                       const string& pathInTree_second) const{
	int counter = 0;
	BOOST_FOREACH(const ptree::value_type &v , pt_.get_child(pathInTree_first))
    {
    	if (indexChild <= counter)
    	  {
    		// take the tree , and ask for the second path inside the tree
    		const ptree p = v.second;
    		FITOB_OUT_LEVEL3( verb() , " FOUND SUBTREE" << v.first);
    		return p.get(pathInTree_second, -1.0);
    		break;
    	  }
    	else
    	{counter++;}
    }
    return -1.0; //Wall
}

const string XMLConfiguration::getStringConfiguration(
		                            const string& pathInTree_first,
		                            const int indexChild ,
		                            const string& pathInTree_second) const{
	int counter = 0;
	BOOST_FOREACH(const ptree::value_type &v , pt_.get_child(pathInTree_first))
    {
    	if (indexChild <= counter)
    	  {
    		// take the tree , and ask for the second path inside the tree
    		const ptree p = v.second;
    		FITOB_OUT_LEVEL3( verb() , " FOUND SUBTREE" << v.first );
    		return p.get(pathInTree_second, "");
    		break;
    	  }
    	else
    	{counter++;}
    }
    return ""; //Wall
}


/** splits one string in different substrings using the delimiter char "delim" */
void XMLConfiguration::splitString(
		const std::string &s, char delim, std::vector<std::string> &elems) const {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
}

/** converts one array of strings into numbers (Integers)*/
void XMLConfiguration::vectStringToIArray(
	   const std::vector<std::string> &elems , IVector& intArr) const {
	   using boost::lexical_cast;
	   using boost::bad_lexical_cast;
       for (unsigned int i = 0 ; i < elems.size() ; i++ ){
    	   //FITOB_OUT_LEVEL3( 4 , "XMLConfiguration::vectStringToIArray:" << elems[i] );
    	   intArr.push_back(lexical_cast<int>(elems[i]));
       }
}

/** converts one array of strings into numbers (Integers)*/
void XMLConfiguration::vectStringToDArray(
		const std::vector<std::string> &elems , DVector& doubleArr) const {
	 using boost::lexical_cast;
	 using boost::bad_lexical_cast;
	 for (unsigned int i = 0 ; i < elems.size() ; i++ ){
  	     //FITOB_OUT_LEVEL3( 4 , "XMLConfiguration::vectStringToDArray:" << elems[i] );
		 doubleArr.push_back(lexical_cast<double>(elems[i]));
	 }
}

void XMLConfiguration::getIntVectorConfiguration(
		                       const string& pathInTree,
		                       char delim,
		                       IVector& outVect) const {

	string tmp = getStringConfiguration(pathInTree);
	std::vector<std::string> tmp_strs;
	splitString(tmp, delim, tmp_strs);
	vectStringToIArray( tmp_strs , outVect );
}

void XMLConfiguration::getDoubleVectorConfiguration(
		                          const string& pathInTree,
		                          char delim,
		                          DVector& outVect) const {

	string tmp = getStringConfiguration(pathInTree);
	std::vector<std::string> tmp_strs;
	splitString(tmp, delim, tmp_strs);
	vectStringToDArray( tmp_strs , outVect );

}

void  XMLConfiguration::getIntVectorConfiguration(
		                       const string& pathInTree_first,
                               const int indexChild ,
                               const string& pathInTree_second,
		                       char delim,
		                       IVector& outVect) const {

	string tmp = getStringConfiguration(pathInTree_first,indexChild,pathInTree_second);
	std::vector<std::string> tmp_strs;
	splitString(tmp, delim, tmp_strs);
	vectStringToIArray( tmp_strs , outVect );

}

void  XMLConfiguration::getDoubleVectorConfiguration(
		                          const string& pathInTree_first,
                                  const int indexChild ,
                                  const string& pathInTree_second,
		                          char delim,
		                          DVector& outVect) const {

	string tmp = getStringConfiguration(pathInTree_first,indexChild,pathInTree_second);
	std::vector<std::string> tmp_strs;
	splitString(tmp, delim, tmp_strs);
	vectStringToDArray( tmp_strs , outVect );

}

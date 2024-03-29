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
 * AbstractCombiGrid.hpp
 *
 *  Created on: Feb 23, 2011
 *      Author: benk
 */

#ifndef ABSTRACTCOMBIGRID_HPP_
#define ABSTRACTCOMBIGRID_HPP_

#include "combigrid/utils/combigrid_ultils.hpp"
#include "combigrid/fullgrid/CombiFullGrid.hpp"
#include "combigrid/combigridkernel/CombiGridKernel.hpp"
#include "combigrid/combischeme/CombiSchemeBasis.hpp"
#include "combigrid/domain/CombiGridDomain.hpp"

namespace combigrid{

	/** full grid type definition */
    typedef FullGrid< double > FullGridD;

    /** combi kernel type definition */
    typedef CombiGridKernel< double > CombiGridKernelD;

	/** The virtual class which defines the interface of the combi grid. <br>
	 * This interface can be implemented by several implementations, serial, OpenMP or even MPI */
	class AbstractCombiGrid {

	public:

    	/** Ctor
    	 * @param combischeme combi schme of the combi grid
    	 * @param hasBoundaryPts array of flag to indicate weather we have boundary points in the dimension*/
		AbstractCombiGrid(const CombiSchemeBasis* combischeme ,
				const std::vector<bool>& hasBoundaryPts ) :gridDomain_(0){
			// create the combi kernel which has the non initialized
			combischeme_ = new CombiSchemeBasis(*combischeme);
			combikernel_ = new CombiGridKernelD( combischeme_ , hasBoundaryPts );
		}

    	/** Ctor
    	 * @param combischeme combi schme of the combi grid
    	 * @param hasBoundaryPts array of flag to indicate weather we have boundary points in the dimension*/
		AbstractCombiGrid(const  CombiSchemeBasis* combischeme ,
				bool hasBoundaryPts = true ) : gridDomain_(0){
			// create the combi kernel which has the non initialized
			std::vector<bool> boundaryTmp( combischeme->getDim() , hasBoundaryPts);
			combischeme_ = new CombiSchemeBasis(*combischeme);
			combikernel_ = new CombiGridKernelD( combischeme_ , boundaryTmp );
		}

		/** Dtor which deletes the kernel */
		~AbstractCombiGrid(){
			// delete the kernel, and with it all the full grids
			delete combikernel_;
		}

		/** sets the domain of the combi grid, this sets globaly for the combi grid,
		 *  and not for each full grid */
		void setDomain( GridDomain* gridDomain ) const { gridDomain_ = gridDomain; }

		/** sets the domain of all the full grids, this is the correct way for extrapolation */
		virtual void setDomainAllFG( GridDomain* gridDomain ) const = 0;

		/** returns the domain of the combi grid */
		const GridDomain* getDomain() const { return gridDomain_; }

		/** returns the array, size of the dimension, for each dimensions indicates weather there are
		 * boundary points */
		const std::vector<bool>& getBoundaryFlags() const { return combikernel_->getBoundaryFlags(); }

		/** create the actual vector for the full grids. <br>
		 * This is be different for serial and parallel implementations */
		virtual void createFullGrids() = 0;

		/** returns one full grid . <br>
		 * In case of the serial implementation this will be the whole
		 * combination scheme, but in the
		 * @param i is the local index of the full grid */
		virtual FullGridD* getFullGrid( int i ) = 0;

		/** same as the previous method but with const environment */
		virtual const FullGridD* getFullGrid( int i ) const = 0;

		/** get the number of full grids. <br>
		 * In serial case this is simple but in parallel (MPI) case this might be only a part
		 * of the combi scheme */
		virtual int getNrFullGrid() const = 0;

		/** evaluate the combi grid at one specified point
		 * @param coords , the coordinates on the unit form [0,1]^D */
		virtual double eval( std::vector<double>& coords ) const = 0;

		/** evaluate the combi grid at one specified point. Buffered evaluation
		 * which might be faster than one evaluation point.
		 * @param coords , the coordinates on the unit form [0,1]^D
		 * @param results , the result vector */
		virtual void eval( std::vector< std::vector<double> >& coords , std::vector<double>& results ) const = 0;

		/** return the combi scheme */
		inline const CombiSchemeBasis* getCombiScheme() const { return combischeme_; }

		/** return the combi kernel, needed for the combination of the full grids */
		inline const CombiGridKernelD* getCombiKernel() const { return combikernel_; }
	protected:



		/** pointer to the combi scheme which is set from the constructor argument */
		CombiSchemeBasis* combischeme_;

		/** the combi kernel which stores the full grids */
		CombiGridKernelD* combikernel_;

		/** grid domain for the combi grid */
		mutable GridDomain* gridDomain_;

	};
}


#endif /* ABSTRACTCOMBIGRID_HPP_ */

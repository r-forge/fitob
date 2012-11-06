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
 * fitobOpenMP.hpp
 *
 *  Created on: Dec 29, 2010
 *      Author: benk
 */

#ifndef FITOBOPENMP_HPP_
#define FITOBOPENMP_HPP_

// define this but only for coding purposes
//#ifndef FITOB_OPENMP
//#define FITOB_OPENMP
//#endif
//#ifndef FITOB_OPENMP_COMBI
//#define FITOB_OPENMP_COMBI
//#endif

namespace fitob{
#ifdef FITOB_OPENMP
// ---------------------- CODE WITH OPENMP -------------------------

#include <omp.h>

// if the combi is not parallel then we set the solver(MG) for parallel
// at the same time (now) both can not work simultaneously
#ifndef FITOB_OPENMP_COMBI
#define FITOB_OPENMP_SOLVER
#endif

      /** get the number of the thread */
      inline int FITOB_OPENMP_GET_THREAD_NUM() { return omp_get_thread_num(); }

      /** set the number of threads
       * This method should be only called at the very beginning
       * (since all the object create themselves based on the the number of threads , create temporary buffer)*/
      inline void FITOB_OPENMP_SET_THREAD_NUM(int nrThread) { omp_set_num_threads(nrThread); }

      /** get the maximal number of threads*/
      inline int FITOB_OPENMP_GET_MAX_THREADS() { return omp_get_max_threads(); }

#else

// no parallel Gauss-Seidl or Combi solving

// ---------------------- CODE WITHOUT OPENMP -------------------------

      /** get the number of the thread */
      inline int FITOB_OPENMP_GET_THREAD_NUM() { return 0; }

      /** set the number of threads */
      inline void FITOB_OPENMP_SET_THREAD_NUM(int nrThread) { /* NOP */ }

      /** get the maximal number of threads*/
      inline int FITOB_OPENMP_GET_MAX_THREADS() { return 1; }

#endif
}

#endif /* FITOBOPENMP_HPP_ */

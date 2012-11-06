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
 * fitobMPI.hpp
 *
 *  Created on: Jan 15, 2012
 *      Author: benk
 */

#ifndef FITOBMPI_HPP_
#define FITOBMPI_HPP_

#include <boost/assert.hpp>


namespace fitob{
#ifdef FITOB_MPI

// ---------------------- CODE WITH MPI -------------------------

#include <mpi.h>

   inline  int  FITOB_MPI_Init( int *argc, char ***argv ){
	   int provided = 0;
	   int ret = MPI_Init_thread( argc, argv , MPI_THREAD_MULTIPLE, &provided);
	   if( provided != MPI_THREAD_MULTIPLE) {
		   std::cout << std::endl << "ERROR: MPI_Init_thread , provided should be " << MPI_THREAD_MULTIPLE <<" , but it is " << provided << std::endl;
		   BOOST_ASSERT(false);
	   }
	   return ret;
   }

   inline  int  FITOB_MPI_Finalize( ) {
	   return MPI_Finalize( );
   }


   inline  int  FITOB_MPI_Comm_size() {
	   int ssi;
	   MPI_Comm_size(MPI_COMM_WORLD, &ssi);
	   return ssi;
   }

   inline  int  FITOB_MPI_Comm_rank() {
	   int rankMy;
	   MPI_Comm_rank(MPI_COMM_WORLD, &rankMy);
	   return rankMy;
   }

   inline  int  FITOB_MPI_Send_Int(void *buf, int count, int dest )
   {
	   return MPI_Send( buf , count, MPI_INT, dest, 0 , MPI_COMM_WORLD );
   }

   inline  int  FITOB_MPI_Send_Long(void *buf, int count, int dest )
   {
	   return MPI_Send( buf, count, MPI_LONG, dest, 0 , MPI_COMM_WORLD );
   }

   inline  int  FITOB_MPI_Send_Double(void *buf, int count, int dest )
   {
	   return MPI_Send( buf, count, MPI_DOUBLE, dest, 0 , MPI_COMM_WORLD );
   }

   inline  int  FITOB_MPI_Recv_Int(void *buf, int count, int source ){
	   MPI_Status status;
       return MPI_Recv( buf, count, MPI_INT , source, 0 , MPI_COMM_WORLD , &status);
   }

   inline  int  FITOB_MPI_Recv_Long(void *buf, int count, int source ){
	   MPI_Status status;
       return MPI_Recv( buf, count, MPI_LONG , source, 0 , MPI_COMM_WORLD , &status);
   }

   inline  int  FITOB_MPI_Recv_Double(void *buf, int count, int source ){
	   MPI_Status status;
       return MPI_Recv( buf, count, MPI_DOUBLE , source, 0 , MPI_COMM_WORLD , &status);
   }

   inline  int  FITOB_MPI_Recv_Int_AnySource(void *buf, int count, int& source ){
	   MPI_Status status;
	   int ret = MPI_Recv( buf, count, MPI_INT , MPI_ANY_SOURCE , 0 , MPI_COMM_WORLD , &status);
	   source = status.MPI_SOURCE;
	   return ret;
   }

   inline  int  FITOB_MPI_Recv_Long_AnySource(void *buf, int count, int& source ){
	   MPI_Status status;
       int ret = MPI_Recv( buf, count, MPI_LONG , MPI_ANY_SOURCE , 0 , MPI_COMM_WORLD , &status);
       source = status.MPI_SOURCE;
	   return ret;
   }

   inline  int  FITOB_MPI_Recv_Double_AnySource(void *buf, int count, int& source ){
	   MPI_Status status;
	   int ret = MPI_Recv( buf, count, MPI_DOUBLE , MPI_ANY_SOURCE , 0 , MPI_COMM_WORLD , &status);
	   source = status.MPI_SOURCE;
	   return ret;
   }

   inline  int  FITOB_MPI_Bcast_Double( void *buffer , int count , int root ){
	   return MPI_Bcast( buffer, count, MPI_DOUBLE , root, MPI_COMM_WORLD );
   }


   inline  int  FITOB_MPI_Barrier(){
	   return MPI_Barrier( MPI_COMM_WORLD );
   }

#else


// ---------------------- CODE WITHOUT MPI -------------------------

   inline  int  FITOB_MPI_Init( int *argc, char ***argv ) { return 0;}

   inline  int  FITOB_MPI_Finalize( ) { return 0;}

   inline  int  FITOB_MPI_Comm_size( ) { return 1;}

   inline  int  FITOB_MPI_Comm_rank( ) { return 0;}

   inline  int  FITOB_MPI_Send_Int(void *buf, int count, int dest ) { return 0;}

   inline  int  FITOB_MPI_Send_Long(void *buf, int count, int dest ){ return 0;}

   inline  int  FITOB_MPI_Send_Double(void *buf, int count, int dest ){ return 0;}

   inline  int  FITOB_MPI_Recv_Int(void *buf, int count, int source ){ return 0;}

   inline  int  FITOB_MPI_Recv_Long(void *buf, int count, int source ){ return 0;}

   inline  int  FITOB_MPI_Recv_Double(void *buf, int count, int source ){ return 0;}

   inline  int  FITOB_MPI_Recv_Int_AnySource(void *buf, int count, int& source ){ return 0;}

   inline  int  FITOB_MPI_Recv_Long_AnySource(void *buf, int count, int& source ){ return 0;}

   inline  int  FITOB_MPI_Recv_Double_AnySource(void *buf, int count, int& source ){ return 0;}

   inline  int  FITOB_MPI_Bcast_Double( void *buffer , int count , int root ){ return 0; }

   inline  int  FITOB_MPI_Barrier(){ return 0;}
#endif
}

#endif /* FITOBMPI_HPP_ */

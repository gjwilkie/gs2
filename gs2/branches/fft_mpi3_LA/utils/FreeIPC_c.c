/*

FreeIPC - A library to allow Fortran programs that use MPI
          to easily access shared memory within multicore nodes
Date - 23rd Oct Version 0.0
Copyright(C) 2009 Ian Bush

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 3.0 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


This library is released under LGPL. The full text of LGPL may be found
at http://www.gnu.org/licenses/lgpl-3.0.txt

In particular, the terms of the GNU LGPL imply that you CAN use this library with
proprietary code with some restrictions (you can link, but not use the code directly). 
Please consult the above URL for exact details.

*/

/*
  The shared memory facilities are provide by use of the standard System V
  IPC facilities. See 
  
  http://www.opengroup.org/onlinepubs/009695399/functions/xsh_chap02_07.html#tag_02_07
  
  for details
*/


/* These routines are mostly very thin wrappers for System V IPC routines.
The main work occurs in the Fortran. The reason these wrappers are being used
( instead we might call the routines direct from Fortran ) is that Fortran 
interop with C does define what you need to use to be interoperable with such
things like size_t, key_t, pid_t etc. Therefore call the wrappers using
basic int and long and then let the prototypes in scope weave their magic
to make sure we get the correct types */

/* Apparently this is required to stop some warning from the gnu compiler */
#if !defined __USE_SVID && !defined __USE_XOPEN && __GNUC__ >= 2
#define _SVID_SOURCE
#define _XOPEN_SOURCE 600
#endif

#define USE_POSIX_SHM
#define DBG 1

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include <sys/types.h>

#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <semaphore.h>
#include <string.h>
#include <errno.h>





void  fipc_get_errval      ( int *eacces, int *eexist, int *einval, int *enfile, 
                             int *enoent, int *enomem, int *enospc );
int fipc_get_seg           ( int fcomm, long size, long disp, void **shmaddr, int *shmid );
int   fipc_crit_start      ( int semid, int rank );
int   fipc_crit_end        ( int semid, int rank );
int   fipc_sizeof_c_int    ( void );
int   fipc_sizeof_c_long   ( void );
int   fipc_sizeof_c_double ( void );
int   fipc_sizeof_c_complex( void );

#include "mpi.h"

#define SEGMENT   1
#define SEMAPHORE 2

#define INVALID_ACTION  2
#define SANITY_FAILED   3
#define ERROR_ELSEWHERE 4


void  fipc_get_errval ( int *eacces, int *eexist, int *einval, int *enfile, 
                        int *enoent, int *enomem, int *enospc ) {

  /* Get the values of the macros used for marking Sys V errors */

  *eacces = - EACCES;
  *eexist = - EEXIST;
  *einval = - EINVAL;
  *enfile = - ENFILE;
  *enoent = - ENOENT;
  //  *enomem = - ENOMEM;
  *enospc = - ENOSPC;

}




int fipc_get_seg( int fcomm, long size, long disp, void **shmaddr, int *shmid ) {

  /* 
     Attempt to create a shared memory segment of size size bytes. If create
     is set create a new one. If exclusive is set create it in exclusive mode.
     The segment will have permissions given by perms

     On a positive value means success, a negative one an error. On error the value
     corresponds to the negative of the errno set by shmget
  */

  MPI_Aint seg;
  MPI_Win win;
  void *mem;
  int rank, d;

  MPI_Comm comm = MPI_Comm_f2c(fcomm);


  MPI_Comm_rank(comm, &rank);

  if (rank == 0)
     seg = size;
   else
     seg = 0;

  int ierr = MPI_Win_allocate_shared(seg, (int) disp, MPI_INFO_NULL, comm, &mem, &win);
 
  if( ierr != MPI_SUCCESS ) {
    int len;
    char *msg;
    fprintf(stderr,"allocate shared failed \n");
    MPI_Error_string(ierr, msg, &len);
    fprintf(stderr,"%s\n",msg);
    return ierr;
  }



  ierr = MPI_Win_shared_query(win, 0, &seg , &d, shmaddr);
  if( ierr != MPI_SUCCESS ) {
    int len;
    char *msg;
    fprintf(stderr,"shared query failed \n");
    MPI_Error_string(ierr, msg, &len);
    fprintf(stderr,"%s\n",msg);
    return ierr;
  }

 // printf("get seg %d %ld %d %ld %p %p\n", rank, seg, d, disp, shmaddr, mem);

  *shmid = win;
  return 0;

}


int fipc_crit_start( int key, int rank ) {

  //sem_init_keys();
  return MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, key);  
}

int fipc_crit_end( int key, int rank ) {
  return MPI_Win_unlock( 0, key);
}

/* Stupid functions to work around problems on stupid 
   broken Cray XT4 series */

int fipc_sizeof_c_int( void ) {

  return sizeof( int );

}

int fipc_sizeof_c_long( void ) {

  return sizeof( long );

}

int fipc_sizeof_c_double( void ) {

  return sizeof( double );

}


int fipc_sizeof_c_complex( void ) {

  return 2 * sizeof( double );

}


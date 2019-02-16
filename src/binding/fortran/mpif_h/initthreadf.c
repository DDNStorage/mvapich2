/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*  
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 *
 * This file is automatically generated by buildiface 
 * DO NOT EDIT
 */
#include "mpi_fortimpl.h"


/* Begin MPI profiling block */
#if defined(USE_WEAK_SYMBOLS) && !defined(USE_ONLY_MPI_NAMES) 
#if defined(HAVE_MULTIPLE_PRAGMA_WEAK)
extern FORT_DLL_SPEC void FORT_CALL MPI_INIT_THREAD( MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread__( MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread( MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread_( MPI_Fint *, MPI_Fint *, MPI_Fint * );

#if defined(F77_NAME_UPPER)
#pragma weak MPI_INIT_THREAD = PMPI_INIT_THREAD
#pragma weak mpi_init_thread__ = PMPI_INIT_THREAD
#pragma weak mpi_init_thread_ = PMPI_INIT_THREAD
#pragma weak mpi_init_thread = PMPI_INIT_THREAD
#elif defined(F77_NAME_LOWER_2USCORE)
#pragma weak MPI_INIT_THREAD = pmpi_init_thread__
#pragma weak mpi_init_thread__ = pmpi_init_thread__
#pragma weak mpi_init_thread_ = pmpi_init_thread__
#pragma weak mpi_init_thread = pmpi_init_thread__
#elif defined(F77_NAME_LOWER_USCORE)
#pragma weak MPI_INIT_THREAD = pmpi_init_thread_
#pragma weak mpi_init_thread__ = pmpi_init_thread_
#pragma weak mpi_init_thread_ = pmpi_init_thread_
#pragma weak mpi_init_thread = pmpi_init_thread_
#else
#pragma weak MPI_INIT_THREAD = pmpi_init_thread
#pragma weak mpi_init_thread__ = pmpi_init_thread
#pragma weak mpi_init_thread_ = pmpi_init_thread
#pragma weak mpi_init_thread = pmpi_init_thread
#endif



#elif defined(HAVE_PRAGMA_WEAK)

#if defined(F77_NAME_UPPER)
extern FORT_DLL_SPEC void FORT_CALL MPI_INIT_THREAD( MPI_Fint *, MPI_Fint *, MPI_Fint * );

#pragma weak MPI_INIT_THREAD = PMPI_INIT_THREAD
#elif defined(F77_NAME_LOWER_2USCORE)
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread__( MPI_Fint *, MPI_Fint *, MPI_Fint * );

#pragma weak mpi_init_thread__ = pmpi_init_thread__
#elif !defined(F77_NAME_LOWER_USCORE)
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread( MPI_Fint *, MPI_Fint *, MPI_Fint * );

#pragma weak mpi_init_thread = pmpi_init_thread
#else
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread_( MPI_Fint *, MPI_Fint *, MPI_Fint * );

#pragma weak mpi_init_thread_ = pmpi_init_thread_
#endif

#elif defined(HAVE_PRAGMA_HP_SEC_DEF)
#if defined(F77_NAME_UPPER)
#pragma _HP_SECONDARY_DEF PMPI_INIT_THREAD  MPI_INIT_THREAD
#elif defined(F77_NAME_LOWER_2USCORE)
#pragma _HP_SECONDARY_DEF pmpi_init_thread__  mpi_init_thread__
#elif !defined(F77_NAME_LOWER_USCORE)
#pragma _HP_SECONDARY_DEF pmpi_init_thread  mpi_init_thread
#else
#pragma _HP_SECONDARY_DEF pmpi_init_thread_  mpi_init_thread_
#endif

#elif defined(HAVE_PRAGMA_CRI_DUP)
#if defined(F77_NAME_UPPER)
#pragma _CRI duplicate MPI_INIT_THREAD as PMPI_INIT_THREAD
#elif defined(F77_NAME_LOWER_2USCORE)
#pragma _CRI duplicate mpi_init_thread__ as pmpi_init_thread__
#elif !defined(F77_NAME_LOWER_USCORE)
#pragma _CRI duplicate mpi_init_thread as pmpi_init_thread
#else
#pragma _CRI duplicate mpi_init_thread_ as pmpi_init_thread_
#endif

#elif defined(HAVE_WEAK_ATTRIBUTE)
#if defined(F77_NAME_UPPER)
extern FORT_DLL_SPEC void FORT_CALL MPI_INIT_THREAD( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("PMPI_INIT_THREAD")));
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread__( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("PMPI_INIT_THREAD")));
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread_( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("PMPI_INIT_THREAD")));
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("PMPI_INIT_THREAD")));

#elif defined(F77_NAME_LOWER_2USCORE)
extern FORT_DLL_SPEC void FORT_CALL MPI_INIT_THREAD( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_init_thread__")));
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread__( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_init_thread__")));
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread_( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_init_thread__")));
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_init_thread__")));

#elif defined(F77_NAME_LOWER_USCORE)
extern FORT_DLL_SPEC void FORT_CALL MPI_INIT_THREAD( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_init_thread_")));
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread__( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_init_thread_")));
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread_( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_init_thread_")));
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_init_thread_")));

#else
extern FORT_DLL_SPEC void FORT_CALL MPI_INIT_THREAD( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_init_thread")));
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread__( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_init_thread")));
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread_( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_init_thread")));
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_init_thread")));

#endif
#endif /* HAVE_PRAGMA_WEAK */
#endif /* USE_WEAK_SYMBOLS */
/* End MPI profiling block */


/* These definitions are used only for generating the Fortran wrappers */
#if defined(USE_WEAK_SYMBOLS) && defined(USE_ONLY_MPI_NAMES)
#if defined(HAVE_MULTIPLE_PRAGMA_WEAK)
extern FORT_DLL_SPEC void FORT_CALL MPI_INIT_THREAD( MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread__( MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread( MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread_( MPI_Fint *, MPI_Fint *, MPI_Fint * );

#if defined(F77_NAME_UPPER)
#pragma weak mpi_init_thread__ = MPI_INIT_THREAD
#pragma weak mpi_init_thread_ = MPI_INIT_THREAD
#pragma weak mpi_init_thread = MPI_INIT_THREAD
#elif defined(F77_NAME_LOWER_2USCORE)
#pragma weak MPI_INIT_THREAD = mpi_init_thread__
#pragma weak mpi_init_thread_ = mpi_init_thread__
#pragma weak mpi_init_thread = mpi_init_thread__
#elif defined(F77_NAME_LOWER_USCORE)
#pragma weak MPI_INIT_THREAD = mpi_init_thread_
#pragma weak mpi_init_thread__ = mpi_init_thread_
#pragma weak mpi_init_thread = mpi_init_thread_
#else
#pragma weak MPI_INIT_THREAD = mpi_init_thread
#pragma weak mpi_init_thread__ = mpi_init_thread
#pragma weak mpi_init_thread_ = mpi_init_thread
#endif
#elif defined(HAVE_WEAK_ATTRIBUTE)
#if defined(F77_NAME_UPPER)
extern FORT_DLL_SPEC void FORT_CALL MPI_INIT_THREAD( MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread__( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("MPI_INIT_THREAD")));
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread_( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("MPI_INIT_THREAD")));
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("MPI_INIT_THREAD")));

#elif defined(F77_NAME_LOWER_2USCORE)
extern FORT_DLL_SPEC void FORT_CALL MPI_INIT_THREAD( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_init_thread__")));
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread__( MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread_( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_init_thread__")));
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_init_thread__")));

#elif defined(F77_NAME_LOWER_USCORE)
extern FORT_DLL_SPEC void FORT_CALL MPI_INIT_THREAD( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_init_thread_")));
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread__( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_init_thread_")));
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread_( MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_init_thread_")));

#else
extern FORT_DLL_SPEC void FORT_CALL MPI_INIT_THREAD( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_init_thread")));
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread__( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_init_thread")));
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread_( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_init_thread")));
extern FORT_DLL_SPEC void FORT_CALL mpi_init_thread( MPI_Fint *, MPI_Fint *, MPI_Fint * );

#endif
#endif

#endif

/* Map the name to the correct form */
#ifndef MPICH_MPI_FROM_PMPI
#if defined(USE_WEAK_SYMBOLS)
#if defined(HAVE_MULTIPLE_PRAGMA_WEAK)
/* Define the weak versions of the PMPI routine*/
#ifndef F77_NAME_UPPER
extern FORT_DLL_SPEC void FORT_CALL PMPI_INIT_THREAD( MPI_Fint *, MPI_Fint *, MPI_Fint * );
#endif
#ifndef F77_NAME_LOWER_2USCORE
extern FORT_DLL_SPEC void FORT_CALL pmpi_init_thread__( MPI_Fint *, MPI_Fint *, MPI_Fint * );
#endif
#ifndef F77_NAME_LOWER_USCORE
extern FORT_DLL_SPEC void FORT_CALL pmpi_init_thread_( MPI_Fint *, MPI_Fint *, MPI_Fint * );
#endif
#ifndef F77_NAME_LOWER
extern FORT_DLL_SPEC void FORT_CALL pmpi_init_thread( MPI_Fint *, MPI_Fint *, MPI_Fint * );

#endif

#if defined(F77_NAME_UPPER)
#pragma weak pmpi_init_thread__ = PMPI_INIT_THREAD
#pragma weak pmpi_init_thread_ = PMPI_INIT_THREAD
#pragma weak pmpi_init_thread = PMPI_INIT_THREAD
#elif defined(F77_NAME_LOWER_2USCORE)
#pragma weak PMPI_INIT_THREAD = pmpi_init_thread__
#pragma weak pmpi_init_thread_ = pmpi_init_thread__
#pragma weak pmpi_init_thread = pmpi_init_thread__
#elif defined(F77_NAME_LOWER_USCORE)
#pragma weak PMPI_INIT_THREAD = pmpi_init_thread_
#pragma weak pmpi_init_thread__ = pmpi_init_thread_
#pragma weak pmpi_init_thread = pmpi_init_thread_
#else
#pragma weak PMPI_INIT_THREAD = pmpi_init_thread
#pragma weak pmpi_init_thread__ = pmpi_init_thread
#pragma weak pmpi_init_thread_ = pmpi_init_thread
#endif /* Test on name mapping */

#elif defined(HAVE_WEAK_ATTRIBUTE)
#if defined(F77_NAME_UPPER)
extern FORT_DLL_SPEC void FORT_CALL pmpi_init_thread__( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("PMPI_INIT_THREAD")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_init_thread_( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("PMPI_INIT_THREAD")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_init_thread( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("PMPI_INIT_THREAD")));

#elif defined(F77_NAME_LOWER_2USCORE)
extern FORT_DLL_SPEC void FORT_CALL PMPI_INIT_THREAD( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_init_thread__")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_init_thread_( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_init_thread__")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_init_thread( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_init_thread__")));

#elif defined(F77_NAME_LOWER_USCORE)
extern FORT_DLL_SPEC void FORT_CALL PMPI_INIT_THREAD( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_init_thread_")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_init_thread__( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_init_thread_")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_init_thread( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_init_thread_")));

#else
extern FORT_DLL_SPEC void FORT_CALL PMPI_INIT_THREAD( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_init_thread")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_init_thread__( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_init_thread")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_init_thread_( MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_init_thread")));

#endif /* Test on name mapping */
#endif /* HAVE_MULTIPLE_PRAGMA_WEAK */
#endif /* USE_WEAK_SYMBOLS */

#ifdef F77_NAME_UPPER
#define mpi_init_thread_ PMPI_INIT_THREAD
#elif defined(F77_NAME_LOWER_2USCORE)
#define mpi_init_thread_ pmpi_init_thread__
#elif !defined(F77_NAME_LOWER_USCORE)
#define mpi_init_thread_ pmpi_init_thread
#else
#define mpi_init_thread_ pmpi_init_thread_
#endif /* Test on name mapping */

#ifdef F77_USE_PMPI
/* This defines the routine that we call, which must be the PMPI version
   since we're renaming the Fortran entry as the pmpi version.  The MPI name
   must be undefined first to prevent any conflicts with previous renamings. */
#undef MPI_Init_thread
#define MPI_Init_thread PMPI_Init_thread 
#endif

#else

#ifdef F77_NAME_UPPER
#define mpi_init_thread_ MPI_INIT_THREAD
#elif defined(F77_NAME_LOWER_2USCORE)
#define mpi_init_thread_ mpi_init_thread__
#elif !defined(F77_NAME_LOWER_USCORE)
#define mpi_init_thread_ mpi_init_thread
/* Else leave name alone */
#endif


#endif /* MPICH_MPI_FROM_PMPI */

/* Prototypes for the Fortran interfaces */
#include "fproto.h"
FORT_DLL_SPEC void FORT_CALL mpi_init_thread_ ( MPI_Fint *v1, MPI_Fint *v2, MPI_Fint *ierr ){
    mpirinitf_(); MPIR_F_NeedInit = 0;
    *ierr = MPI_Init_thread( 0, 0, (int)*v1, v2 );
}

/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

/* Copyright (c) 2001-2018, The Ohio State University. All rights
 * reserved.
 *
 * This file is part of the MVAPICH2 software package developed by the
 * team members of The Ohio State University's Network-Based Computing
 * Laboratory (NBCL), headed by Professor Dhabaleswar K. (DK) Panda.
 *
 * For detailed copyright and licensing information, please refer to the
 * copyright file COPYRIGHT in the top level MVAPICH2 directory.
 *
 */

#include "mpiimpl.h"
#include "coll_shmem.h"
#ifdef MRAIL_GEN2_INTERFACE
#include <cr.h>
#endif

MPIR_T_PVAR_ULONG2_COUNTER_DECL_EXTERN(MV2, mv2_coll_barrier_pairwise);
MPIR_T_PVAR_ULONG2_COUNTER_DECL_EXTERN(MV2, mv2_coll_barrier_shmem);

MPIR_T_PVAR_ULONG2_COUNTER_DECL_EXTERN(MV2, mv2_coll_barrier_pairwise_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL_EXTERN(MV2, mv2_coll_barrier_pairwise_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL_EXTERN(MV2, mv2_coll_barrier_pairwise_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL_EXTERN(MV2, mv2_coll_barrier_pairwise_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL_EXTERN(MV2, mv2_coll_barrier_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL_EXTERN(MV2, mv2_coll_barrier_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL_EXTERN(MV2, mv2_coll_barrier_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL_EXTERN(MV2, mv2_coll_barrier_count_recv);

static int MPIR_Pairwise_Barrier_MV2(MPID_Comm * comm_ptr, MPIR_Errflag_t *errflag)
{

    int size, rank;
    int d, dst, src;
    int mpi_errno = MPI_SUCCESS;

    MPIR_T_PVAR_COUNTER_INC(MV2, mv2_coll_barrier_pairwise, 1);

    size = comm_ptr->local_size;
    /* Trivial barriers return immediately */
    if (size == 1)
        return MPI_SUCCESS;

    rank = comm_ptr->rank;

    /*  N2_prev = greatest power of two < size of Comm  */
    int N2_prev = comm_ptr->dev.ch.gpof2;
    int surfeit = size - N2_prev;

    /* Perform a combine-like operation */
    if (rank < N2_prev) {
        if (rank < surfeit) {
            /* get the fanin letter from the upper "half" process: */
            dst = N2_prev + rank;
            MPIR_PVAR_INC(barrier, pairwise, recv, 0, MPI_BYTE);
            mpi_errno = MPIC_Recv(NULL, 0, MPI_BYTE, dst, MPIR_BARRIER_TAG,
                                     comm_ptr, MPI_STATUS_IGNORE, errflag);
        }

        /* combine on embedded N2_prev power-of-two processes */
        for (d = 1; d < N2_prev; d <<= 1) {
            dst = (rank ^ d);
            MPIR_PVAR_INC(barrier, pairwise, send, 0, MPI_BYTE);
            MPIR_PVAR_INC(barrier, pairwise, recv, 0, MPI_BYTE);
            mpi_errno =
                MPIC_Sendrecv(NULL, 0, MPI_BYTE, dst, MPIR_BARRIER_TAG, NULL,
                                 0, MPI_BYTE, dst, MPIR_BARRIER_TAG, comm_ptr,
                                 MPI_STATUS_IGNORE, errflag);
        }

        /* fanout data to nodes above N2_prev... */
        if (rank < surfeit) {
            dst = N2_prev + rank;
            MPIR_PVAR_INC(barrier, pairwise, send, 0, MPI_BYTE);
            mpi_errno = MPIC_Send(NULL, 0, MPI_BYTE, dst, MPIR_BARRIER_TAG,
                                     comm_ptr, errflag);
        }
    } else {
        /* fanin data to power of 2 subset */
        src = rank - N2_prev;
        MPIR_PVAR_INC(barrier, pairwise, send, 0, MPI_BYTE);
        MPIR_PVAR_INC(barrier, pairwise, recv, 0, MPI_BYTE);
        mpi_errno = MPIC_Sendrecv(NULL, 0, MPI_BYTE, src, MPIR_BARRIER_TAG,
                                     NULL, 0, MPI_BYTE, src, MPIR_BARRIER_TAG,
                                     comm_ptr, MPI_STATUS_IGNORE, errflag);
    }

    return mpi_errno;

}

static int MPIR_shmem_barrier_MV2(MPID_Comm * comm_ptr, MPIR_Errflag_t *errflag)
{

    int mpi_errno = MPI_SUCCESS;

    MPI_Comm shmem_comm = MPI_COMM_NULL, leader_comm = MPI_COMM_NULL;
    MPID_Comm *shmem_commptr = NULL, *leader_commptr = NULL;
    int local_rank = -1, local_size = 0;
    int total_size, shmem_comm_rank;

    MPIR_T_PVAR_COUNTER_INC(MV2, mv2_coll_barrier_shmem, 1);
    MPIR_T_PVAR_COUNTER_INC(MV2, mv2_num_shmem_coll_calls, 1);
    shmem_comm = comm_ptr->dev.ch.shmem_comm;
    leader_comm = comm_ptr->dev.ch.leader_comm;

    total_size = comm_ptr->local_size;
    shmem_comm = comm_ptr->dev.ch.shmem_comm;

    MPID_Comm_get_ptr(shmem_comm, shmem_commptr);
    local_rank = shmem_commptr->rank;
    local_size = shmem_commptr->local_size;
    shmem_comm_rank = shmem_commptr->dev.ch.shmem_comm_rank;
    leader_comm = comm_ptr->dev.ch.leader_comm;
    MPID_Comm_get_ptr(leader_comm, leader_commptr);

    if (local_size > 1) {
        MPIDI_CH3I_SHMEM_COLL_Barrier_gather(local_size, local_rank,
                                             shmem_comm_rank);
    }

    if ((local_rank == 0) && (local_size != total_size)) {
        mpi_errno = MPIR_Pairwise_Barrier_MV2(leader_commptr, errflag);
    }

    if (local_size > 1) {
        MPIDI_CH3I_SHMEM_COLL_Barrier_bcast(local_size, local_rank,
                                            shmem_comm_rank);
    }

    return mpi_errno;

}

/* This is the default implementation of the barrier operation.  The
   algorithm is:
   
   Algorithm: MPI_Barrier

   We use pairwise exchange with recursive doubling algorithm 
   described in:
   R. Gupta, V. Tipparaju, J. Nieplocha and D.K. Panda,
   "Efficient Barrier using Remote Memory Operations on VIA-Based Clusters",
   IEEE Cluster Computing, 2002

   Possible improvements: 

   End Algorithm: MPI_Barrier

   This is an intracommunicator barrier only!
*/

/* not declared static because it is called in ch3_comm_connect/accept */
#undef FUNCNAME
#define FUNCNAME MPIR_Barrier_intra_MV2
#undef FCNAME
#define FCNAME MPL_QUOTE(FUNCNAME)
int MPIR_Barrier_intra_MV2(MPID_Comm * comm_ptr, MPIR_Errflag_t *errflag)
{
    int size;
    int mpi_errno = MPI_SUCCESS;
    int mpi_errno_ret = MPI_SUCCESS;

    size = comm_ptr->local_size;
    /* Trivial barriers return immediately */
    if (size == 1)
        return MPI_SUCCESS;

    /* Only one collective operation per communicator can be active at any
       time */
    MPIDU_ERR_CHECK_MULTIPLE_THREADS_ENTER(comm_ptr);

    if (mv2_enable_shmem_collectives && mv2_enable_shmem_barrier
        && comm_ptr->dev.ch.shmem_coll_ok == 1) {

        mpi_errno = MPIR_shmem_barrier_MV2(comm_ptr, errflag);

    } else {

        mpi_errno = MPIR_Pairwise_Barrier_MV2(comm_ptr, errflag);
    }

    if (mpi_errno) {
        /* for communication errors, just record the error but continue */
        *errflag = MPIR_ERR_GET_CLASS(mpi_errno);
        MPIR_ERR_SET(mpi_errno, MPI_ERR_OTHER, "**fail");
        MPIR_ERR_ADD(mpi_errno_ret, mpi_errno);
    }

    MPIDU_ERR_CHECK_MULTIPLE_THREADS_EXIT(comm_ptr);

    return mpi_errno;
}

#undef FUNCNAME
#define FUNCNAME MPIR_Barrier_MV2
#undef FCNAME
#define FCNAME MPL_QUOTE(FUNCNAME)
int MPIR_Barrier_MV2(MPID_Comm * comm_ptr, MPIR_Errflag_t *errflag)
{
    int mpi_errno = MPI_SUCCESS;
    mpi_errno = MPIR_Barrier_intra_MV2(comm_ptr, errflag);
    if (mpi_errno)
        MPIR_ERR_POP(mpi_errno);

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

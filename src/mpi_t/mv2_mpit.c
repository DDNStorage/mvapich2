/* Copyright (c) 2001-2018, The Ohio State University. All rights
 * reserved.
 *
 * This file is part of the MVAPICH2 software package developed by the
 * team members of The Ohio State University's Network-Based Computing
 * Laboratory (NBCL), headed by Professor Dhabaleswar K. (DK) Panda.
 *
 * For detailed copyright and licensing information, please refer to the
 * copyright file COPYRIGHT in the top level MVAPICH2 directory.
 */
#include <mpichconf.h>
#include <mpiimpl.h>

/*
 * Profile IB channel-manager
 */
MPIR_T_PVAR_ULONG_COUNTER_DECL(MV2, mv2_ibv_channel_ctrl_packet_count);
MPIR_T_PVAR_ULONG_COUNTER_DECL(MV2, mv2_ibv_channel_out_of_order_packet_count);
MPIR_T_PVAR_ULONG_COUNTER_DECL(MV2, mv2_ibv_channel_exact_recv_count);

/*
 * Profile RDMA_FP connections and packets
 */
MPIR_T_PVAR_ULONG_COUNTER_DECL(MV2, mv2_rdmafp_ctrl_packet_count);
MPIR_T_PVAR_ULONG_COUNTER_DECL(MV2, mv2_rdmafp_out_of_order_packet_count);
MPIR_T_PVAR_ULONG_COUNTER_DECL(MV2, mv2_rdmafp_exact_recv_count);

/*
 * Track vbuf usage
 */
MPIR_T_PVAR_ULONG_COUNTER_DECL(MV2, mv2_vbuf_allocated);
MPIR_T_PVAR_ULONG_COUNTER_DECL(MV2, mv2_vbuf_freed);
MPIR_T_PVAR_ULONG_LEVEL_DECL(MV2, mv2_vbuf_available);
MPIR_T_PVAR_ULONG_COUNTER_DECL(MV2, mv2_ud_vbuf_allocated);
MPIR_T_PVAR_ULONG_COUNTER_DECL(MV2, mv2_ud_vbuf_freed);
MPIR_T_PVAR_ULONG_LEVEL_DECL(MV2, mv2_ud_vbuf_available);

/*
 * Track smp memory usage
 */
MPIR_T_PVAR_ULONG_COUNTER_DECL(MV2, mv2_smp_eager_sent);
MPIR_T_PVAR_ULONG_COUNTER_DECL(MV2, mv2_smp_rndv_sent);
MPIR_T_PVAR_ULONG_COUNTER_DECL(MV2, mv2_smp_eager_received);
MPIR_T_PVAR_ULONG_COUNTER_DECL(MV2, mv2_smp_rndv_received);
MPIR_T_PVAR_ULONG_LEVEL_DECL(MV2, mv2_smp_eager_total_buffer);
MPIR_T_PVAR_ULONG_LEVEL_DECL(MV2, mv2_smp_rndv_total_buffer);
MPIR_T_PVAR_ULONG_LEVEL_DECL(MV2, mv2_smp_eager_avail_buffer);
MPIR_T_PVAR_ULONG_LEVEL_DECL(MV2, mv2_smp_rndv_avail_buffer);
MPIR_T_PVAR_ULONG_LEVEL_DECL(MV2, mv2_smp_eager_buffer_max_use);
MPIR_T_PVAR_ULONG_LEVEL_DECL(MV2, mv2_smp_rndv_buffer_max_use);

/*
 * Track the registration cache hits and misses
 */
MPIR_T_PVAR_ULONG_COUNTER_DECL(MV2, mv2_reg_cache_hits);
MPIR_T_PVAR_ULONG_COUNTER_DECL(MV2, mv2_reg_cache_misses);

/*
 * Count progress engine polling
 */
MPIR_T_PVAR_ULONG_COUNTER_DECL(MV2, mpit_progress_poll);
MPIR_T_PVAR_ULONG_COUNTER_DECL(MV2, mv2_smp_read_progress_poll);
MPIR_T_PVAR_ULONG_COUNTER_DECL(MV2, mv2_smp_write_progress_poll);
MPIR_T_PVAR_ULONG_COUNTER_DECL(MV2, mv2_smp_read_progress_poll_success);
MPIR_T_PVAR_ULONG_COUNTER_DECL(MV2, mv2_smp_write_progress_poll_success);

/*
 * Count number of shared-memory collective calls
 */
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_num_shmem_coll_calls);

/*
 * Count 2-level communicator creation requests
 */
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_num_2level_comm_requests);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_num_2level_comm_success);

/*
 * Count MVAPICH Broadcast algorithms used
 */
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_binomial);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_scatter_doubling_allgather);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_scatter_ring_allgather);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_scatter_ring_allgather_shm);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_shmem);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_knomial_internode);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_knomial_intranode);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_mcast_internode);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_pipelined);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_binomial_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_scatter_for_bcast_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_scatter_doubling_allgather_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_scatter_ring_allgather_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_scatter_ring_allgather_shm_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_knomial_internode_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_knomial_intranode_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_mcast_internode_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_pipelined_zcpy_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_tune_inter_node_helper_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_inter_node_helper_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_binomial_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_scatter_for_bcast_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_scatter_doubling_allgather_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_scatter_ring_allgather_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_scatter_ring_allgather_shm_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_knomial_internode_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_knomial_intranode_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_mcast_internode_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_pipelined_zcpy_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_tune_inter_node_helper_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_inter_node_helper_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_binomial_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_scatter_for_bcast_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_scatter_doubling_allgather_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_scatter_ring_allgather_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_scatter_ring_allgather_shm_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_knomial_internode_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_knomial_intranode_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_mcast_internode_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_pipelined_zcpy_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_tune_inter_node_helper_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_inter_node_helper_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_binomial_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_scatter_for_bcast_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_scatter_doubling_allgather_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_scatter_ring_allgather_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_scatter_ring_allgather_shm_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_knomial_internode_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_knomial_intranode_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_mcast_internode_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_pipelined_zcpy_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_tune_inter_node_helper_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_inter_node_helper_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_bcast_count_recv);


/*
 * Count MVAPICH Alltoall algorithms used
 */
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_inplace);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_bruck);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_rd);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_sd);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_pw);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_inplace_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_bruck_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_sd_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_pw_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_intra_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_inplace_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_bruck_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_sd_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_pw_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_intra_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_inplace_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_bruck_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_sd_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_pw_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_intra_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_inplace_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_bruck_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_sd_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_pw_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_intra_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_count_recv);

/*
 * Count MVAPICH Alltoallv algorithms used
 */ 
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoallv_pw); 
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoallv_intra_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoallv_intra_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoallv_intra_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoallv_intra_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoallv_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoallv_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoallv_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoallv_count_recv);

/*
 * Count MVAPICH Alltoall Cuda algorithms used
 */
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_cuda_intra_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_cuda_intra_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_cuda_intra_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_cuda_intra_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_cuda_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_cuda_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_cuda_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_alltoall_cuda_count_recv);

/*
 * Count MVAPICH Allreduce algorithms used
 */
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allreduce_sharp);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allreduce_shm_rd);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allreduce_shm_rs);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allreduce_shm_intra);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allreduce_intra_p2p);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allreduce_2lvl);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allreduce_shmem);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allreduce_mcast);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allreduce_reduce_scatter_allgather_colls);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allreduce_pt2pt_rd_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allreduce_pt2pt_rs_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allreduce_pt2pt_rd_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allreduce_pt2pt_rs_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allreduce_pt2pt_rd_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allreduce_pt2pt_rs_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allreduce_pt2pt_rd_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allreduce_pt2pt_rs_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allreduce_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allreduce_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allreduce_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allreduce_count_recv);

/*
 * Count MVAPICH Allgather algorithms used
 */
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_rd_allgather_comm);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_rd);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_bruck);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_ring);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_direct);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_directspread);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_gather_bcast);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_2lvl_nonblocked);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_2lvl_ring_nonblocked);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_2lvl_direct);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_2lvl_ring);

MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_rd_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_bruck_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_ring_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_direct_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_directspread_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_2lvl_ring_nonblocked_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_2lvl_direct_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_2lvl_ring_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_rd_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_bruck_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_ring_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_direct_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_directspread_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_2lvl_ring_nonblocked_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_2lvl_direct_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_2lvl_ring_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_rd_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_bruck_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_ring_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_direct_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_directspread_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_2lvl_ring_nonblocked_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_2lvl_direct_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_2lvl_ring_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_rd_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_bruck_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_ring_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_direct_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_directspread_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_2lvl_ring_nonblocked_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_2lvl_direct_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_2lvl_ring_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_count_recv);

/*
 * Count MVAPICH AllGather Cuda algorithms used
 */
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_cuda_intra_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_cuda_intra_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_cuda_intra_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_cuda_intra_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_cuda_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_cuda_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_cuda_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgather_cuda_count_recv);


/*
 * Count MVAPICH Gather algorithms used
 */
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gather_pt2pt);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gather_direct);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gather_direct_blk);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gather_two_level_direct);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gather_limic_scheme_pt_pt);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gather_limic_scheme_pt_linear);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gather_limic_scheme_linear_pt);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gather_limic_scheme_linear_linear);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gather_limic_scheme_single_leader);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gather_intra_node_limic);

MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gather_direct_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gather_direct_blk_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gather_two_level_direct_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gather_direct_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gather_direct_blk_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gather_two_level_direct_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gather_direct_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gather_direct_blk_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gather_two_level_direct_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gather_direct_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gather_direct_blk_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gather_two_level_direct_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gather_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gather_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gather_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gather_count_recv);

/*
 * Count MVAPICH Reduce Scatter algorithms used
 */
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_scatter_noncomm);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_scatter_basic);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_scatter_rec_halving);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_scatter_pairwise);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_scatter_non_comm);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_scatter_noncomm_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_scatter_basic_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_scatter_rec_halving_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_scatter_pairwise_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_scatter_non_comm_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_scatter_noncomm_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_scatter_basic_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_scatter_rec_halving_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_scatter_pairwise_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_scatter_non_comm_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_scatter_noncomm_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_scatter_basic_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_scatter_rec_halving_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_scatter_pairwise_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_scatter_non_comm_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_scatter_noncomm_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_scatter_basic_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_scatter_rec_halving_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_scatter_pairwise_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_scatter_non_comm_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_scatter_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_scatter_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_scatter_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_scatter_count_recv);

/*
 * Count MVAPICH Scatter algorithms used
 */
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_mcast);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_binomial);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_direct);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_direct_blk);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_two_level_binomial);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_two_level_direct);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_inter);

/*
 * Count of messages and total bytes communicated by Send and Recv messages in MVAPICH Scatter algorithms
 */
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_mcast_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_mcast_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_binomial_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_binomial_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_direct_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_direct_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_direct_blk_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_direct_blk_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_two_level_binomial_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_two_level_binomial_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_two_level_direct_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_two_level_direct_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_inter_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_inter_bytes_recv);

MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_mcast_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_mcast_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_binomial_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_binomial_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_direct_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_direct_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_direct_blk_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_direct_blk_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_two_level_binomial_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_two_level_binomial_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_two_level_direct_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_two_level_direct_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_inter_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_inter_count_recv);

MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scatter_count_recv);

/*
 * Count MVAPICH Reduce algorithms used
 */
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_binomial);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_redscat_gather);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_shmem);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_knomial);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_zcpy);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_binomial_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_redscat_gather_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_two_level_helper_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_knomial_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_zcpy_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_binomial_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_redscat_gather_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_two_level_helper_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_knomial_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_zcpy_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_binomial_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_redscat_gather_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_two_level_helper_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_knomial_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_zcpy_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_binomial_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_redscat_gather_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_two_level_helper_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_knomial_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_zcpy_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_reduce_count_recv);

/*
 * Count MVAPICH Gatherv algorithms used
 */
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gatherv_algo);

MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gatherv_default_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gatherv_default_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gatherv_default_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gatherv_default_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gatherv_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gatherv_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gatherv_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_gatherv_count_recv);

/*
 * Count MVAPICH Allgatherv algorithms used
 */
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgatherv_rec_doubling);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgatherv_bruck);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgatherv_ring);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgatherv_ring_cyclic);

MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgatherv_rec_doubling_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgatherv_bruck_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgatherv_ring_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgatherv_rec_doubling_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgatherv_bruck_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgatherv_ring_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgatherv_rec_doubling_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgatherv_bruck_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgatherv_ring_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgatherv_rec_doubling_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgatherv_bruck_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgatherv_ring_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgatherv_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgatherv_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgatherv_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_allgatherv_count_recv);

/*
 * Count MVAPICH Barrier algorithms used
 */
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_barrier_pairwise);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_barrier_shmem);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_barrier_pairwise_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_barrier_pairwise_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_barrier_pairwise_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_barrier_pairwise_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_barrier_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_barrier_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_barrier_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_barrier_count_recv);

/*
 * Count MVAPICH Exscan algorithms used
 */
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_exscan_algo);

/*
 * Count MVAPICH Scan algorithms used
 */
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_scan_algo);

/*
 * Count MVAPICH Iallreduce algorithms used
 */
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_iallreduce_sharp);

/*
 * Count number of UD retransmissions
 */
MPIR_T_PVAR_ULONG_COUNTER_DECL(MV2, rdma_ud_retransmissions);

/*
 * Count MVAPICH iScatter algorithms used
 */
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_iscatter_binomial_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_iscatter_binomial_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_iscatter_binomial_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_iscatter_binomial_count_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_iscatter_bytes_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_iscatter_bytes_recv);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_iscatter_count_send);
MPIR_T_PVAR_ULONG2_COUNTER_DECL(MV2, mv2_coll_iscatter_count_recv);


void
MPIT_REGISTER_MV2_VARIABLES (void)
{
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_num_2level_comm_requests,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "Number of 2-level comm creation requests");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_num_2level_comm_success,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "Number of successful 2-level comm creations");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_num_shmem_coll_calls,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "Number of times MV2 shared-memory collective calls were invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG,
            mpit_progress_poll,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "CH3 RDMA progress engine polling count");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG,
            mv2_smp_read_progress_poll,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "CH3 SMP read progress engine polling count");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG,
            mv2_smp_write_progress_poll,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "CH3 SMP write progress engine polling count");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG,
            mv2_smp_read_progress_poll_success,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "Unsucessful CH3 SMP read progress engine polling count");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG,
            mv2_smp_write_progress_poll_success,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "Unsucessful CH3 SMP write progress engine polling count");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG,
            rdma_ud_retransmissions,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "CH3 RDMA UD retransmission count");
			
	/* BEGIN: Register PVARs for Bcast algorithms */
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_binomial,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times MV2 binomial bcast algorithm  was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_scatter_doubling_allgather,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times MV2 scatter+double allgather bcast algorithm "
            "was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_scatter_ring_allgather,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "Number of times MV2 scatter+ring allgather bcast algorithm "
            "was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_scatter_ring_allgather_shm,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times MV2 scatter+ring allgather shm bcast "
            "algorithm was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_shmem,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times MV2 shmem bcast algorithm was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_knomial_internode,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times MV2 knomial internode bcast algorithm "
            "was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_knomial_intranode,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times MV2 knomial intranode bcast algorithm "
            "was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_mcast_internode,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times MV2 mcast internode bcast algorithm "
            "was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_pipelined,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times MV2 pipelined bcast algorithm was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_binomial_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by binomial algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_scatter_for_bcast_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by scatter for bcast algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_scatter_doubling_allgather_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by scatter doubling allgather algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_scatter_ring_allgather_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by scatter ring allgather algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_scatter_ring_allgather_shm_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by scatter ring allgather shm algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_knomial_internode_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by knomial internode algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_knomial_intranode_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by knomial intranode algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_mcast_internode_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by mcast internode algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_pipelined_zcpy_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by pipelined zcpy algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_tune_inter_node_helper_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by tune inter node helper algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_inter_node_helper_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by inter node helper algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_binomial_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by binomial algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_scatter_for_bcast_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by scatter for bcast algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_scatter_doubling_allgather_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by scatter doubling allgather algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_scatter_ring_allgather_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by scatter ring allgather algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_scatter_ring_allgather_shm_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by scatter ring allgather shm algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_knomial_internode_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by knomial internode algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_knomial_intranode_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by knomial intranode algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_mcast_internode_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by mcast internode algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_pipelined_zcpy_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by pipelined zcpy algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_tune_inter_node_helper_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by tune inter node helper algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_inter_node_helper_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by inter node helper algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_binomial_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by binomial algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_scatter_for_bcast_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by scatter for bcast algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_scatter_doubling_allgather_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by scatter doubling allgather algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_scatter_ring_allgather_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by scatter ring allgather algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_scatter_ring_allgather_shm_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by scatter ring allgather shm algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_knomial_internode_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by knomial internode algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_knomial_intranode_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by knomial intranode algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_mcast_internode_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by mcast internode algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_pipelined_zcpy_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by pipelined zcpy algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_tune_inter_node_helper_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by tune inter node helper algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_inter_node_helper_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by inter node helper algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_binomial_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by binomial algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_scatter_for_bcast_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by scatter for bcast algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_scatter_doubling_allgather_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by scatter doubling allgather algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_scatter_ring_allgather_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by scatter ring allgather algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_scatter_ring_allgather_shm_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by scatter ring allgather shm algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_knomial_internode_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by knomial internode algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_knomial_intranode_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by knomial intranode algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_mcast_internode_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by mcast internode algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_pipelined_zcpy_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by pipelined zcpy algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_tune_inter_node_helper_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by tune inter node helper algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_inter_node_helper_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by inter node helper algorithm of bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by bcast collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_bcast_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by bcast collective");

    /* End: Register PVARs for Bcast algorithms */
 
    /* BEGIN: Register PVARs for alltoall algorithms */
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_inplace,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times MV2 in-place alltoall algorithm was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_bruck,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times MV2 brucks alltoall algorithm was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_rd,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times MV2 recursive-doubling alltoall algorithm was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_sd,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times MV2 scatter-destination alltoall algorithm was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_pw,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times MV2 pairwise alltoall algorithm was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_inplace_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by inplace algorithm of alltoall collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_bruck_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by bruck algorithm of alltoall collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_sd_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by sd algorithm of alltoall collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_pw_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by pw algorithm of alltoall collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_inplace_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by inplace algorithm of alltoall collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_bruck_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by bruck algorithm of alltoall collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_sd_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by sd algorithm of alltoall collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_pw_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by pw algorithm of alltoall collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_inplace_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by inplace algorithm of alltoall collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_bruck_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by bruck algorithm of alltoall collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_sd_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by sd algorithm of alltoall collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_pw_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by pw algorithm of alltoall collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_inplace_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by inplace algorithm of alltoall collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_bruck_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by bruck algorithm of alltoall collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_sd_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by sd algorithm of alltoall collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_pw_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by pw algorithm of alltoall collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by alltoall collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by alltoall collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by alltoall collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by alltoall collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_intra_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by intra algorithm of alltoall collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_intra_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by intra algorithm of alltoall collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_intra_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by intra algorithm of alltoall collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_intra_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by intra algorithm of alltoall collective");
    /* End: Register PVARs for alltoall algorithms */
 
    /* BEGIN: Register PVARs for alltoallv algorithms */
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoallv_pw,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times MV2 pairwise alltoallv algorithm was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoallv_intra_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by intra algorithm of alltoallv collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoallv_intra_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by intra algorithm of alltoallv collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoallv_intra_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by intra algorithm of alltoallv collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoallv_intra_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by intra algorithm of alltoallv collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoallv_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by alltoallv collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoallv_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by alltoallv collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoallv_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by alltoallv collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoallv_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by alltoallv collective");
    /* End: Register PVARs for alltoallv algorithms */
	
	/* Begin: Register PVARs for alltoall cuda algorithms */
	MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_cuda_intra_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by intra algorithm of alltoall_cuda collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_cuda_intra_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by intra algorithm of alltoall_cuda collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_cuda_intra_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by intra algorithm of alltoall_cuda collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_cuda_intra_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by intra algorithm of alltoall_cuda collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_cuda_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by alltoall_cuda collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_cuda_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by alltoall_cuda collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_cuda_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by alltoall_cuda collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_alltoall_cuda_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by alltoall_cuda collective");
	/* End: Register PVARs for alltoall cuda algorithms */
 
    /* BEGIN: Register PVARs for allreduce algorithms */
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allreduce_sharp,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times MV2 sharp allreduce algorithm was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allreduce_shm_rd,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times MV2 shm rd allreduce algorithm was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allreduce_shm_rs,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times MV2 shm rs allreduce algorithm was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allreduce_shm_intra,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times MV2 shm intra allreduce algorithm was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allreduce_intra_p2p,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times MV2 intra p2p allreduce algorithm was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allreduce_2lvl,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times MV2 two-level allreduce algorithm was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allreduce_shmem,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times MV2 shmem allreduce algorithm was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allreduce_mcast,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times MV2 multicast-based allreduce algorithm was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allreduce_reduce_scatter_allgather_colls,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times MV2 red-scat allga colls based allreduce algorithm was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allreduce_pt2pt_rd_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by pt2pt rd algorithm of allreduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allreduce_pt2pt_rs_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by pt2pt rs algorithm of allreduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allreduce_pt2pt_rd_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by pt2pt rd algorithm of allreduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allreduce_pt2pt_rs_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by pt2pt rs algorithm of allreduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allreduce_pt2pt_rd_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by pt2pt rd algorithm of allreduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allreduce_pt2pt_rs_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by pt2pt rs algorithm of allreduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allreduce_pt2pt_rd_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by pt2pt rd algorithm of allreduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allreduce_pt2pt_rs_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by pt2pt rs algorithm of allreduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allreduce_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by allreduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allreduce_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by allreduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allreduce_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by allreduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allreduce_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by allreduce collective");
	/* End: Register PVARs for AllReduce algorithms */
	
    /* Begin: Register PVARs for Allgather algorithms */
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_rd_allgather_comm,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times optimzied recursive doubling Allgather was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_rd,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times recursive doubling Allgather was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_bruck,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times bruck Allgather was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_ring,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times ring Allgather was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_direct,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times direct Allgather was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_directspread,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times Allgather direct spread was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_gather_bcast,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times gather bcast Allgather was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_2lvl_nonblocked,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times Allgather for non-block process mapping was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_2lvl_ring_nonblocked,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times ring Allgather for non-block process mapping was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_2lvl_direct,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times direct Allgather for non-block process mapping was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_2lvl_ring,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times two level ring Allgather was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_rd_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by rd algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_bruck_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by bruck algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_ring_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by ring algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_direct_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by direct algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_directspread_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by directspread algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_2lvl_ring_nonblocked_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by 2lvl ring nonblocked algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_2lvl_direct_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by 2lvl direct algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_2lvl_ring_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by 2lvl ring algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_rd_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by rd algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_bruck_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by bruck algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_ring_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by ring algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_direct_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by direct algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_directspread_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by directspread algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_2lvl_ring_nonblocked_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by 2lvl ring nonblocked algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_2lvl_direct_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by 2lvl direct algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_2lvl_ring_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by 2lvl ring algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_rd_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by rd algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_bruck_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by bruck algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_ring_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by ring algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_direct_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by direct algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_directspread_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by directspread algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_2lvl_ring_nonblocked_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by 2lvl ring nonblocked algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_2lvl_direct_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by 2lvl direct algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_2lvl_ring_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by 2lvl ring algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_rd_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by rd algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_bruck_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by bruck algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_ring_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by ring algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_direct_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by direct algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_directspread_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by directspread algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_2lvl_ring_nonblocked_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by 2lvl ring nonblocked algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_2lvl_direct_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by 2lvl direct algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_2lvl_ring_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by 2lvl ring algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by default algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by default algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by default algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by default algorithm of allgather collective");			
    /* End: Register PVARs for Allgather algorithms */
	
	/* Begin: Register PVARs for AllGather CUDA algorithms */
	 MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_cuda_intra_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by cuda intra algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_cuda_intra_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by cuda intra algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_cuda_intra_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by cuda intra algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_cuda_intra_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by cuda intra algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_cuda_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by cuda algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_cuda_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by cuda algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_cuda_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by cuda algorithm of allgather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgather_cuda_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by cuda algorithm of allgather collective");
	/* End: Register PVARs for AllGather CUDA algorithms */	
	
    /* Begin: Register PVARs for Gather algorithms */
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gather_pt2pt,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times point to point Gather was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gather_direct,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times direct gather was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gather_direct_blk,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times direct blk gather was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gather_two_level_direct,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times two level direct gather was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gather_limic_scheme_pt_pt,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times point to point limic gather was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gather_limic_scheme_pt_linear,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times pt_linear gather was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gather_limic_scheme_linear_pt,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times linear_pt gather was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gather_limic_scheme_linear_linear,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times linear limic gather was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gather_limic_scheme_single_leader,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times single leader limic gather was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gather_intra_node_limic,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times intra-node limic gather was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gather_direct_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by direct algorithm of gather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gather_direct_blk_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by direct blk algorithm of gather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gather_two_level_direct_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by two level direct algorithm of gather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gather_direct_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by direct algorithm of gather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gather_direct_blk_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by direct blk algorithm of gather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gather_two_level_direct_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by two level direct algorithm of gather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gather_direct_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by direct algorithm of gather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gather_direct_blk_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by direct blk algorithm of gather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gather_two_level_direct_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by two level direct algorithm of gather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gather_direct_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by direct algorithm of gather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gather_direct_blk_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by direct blk algorithm of gather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gather_two_level_direct_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by two level direct algorithm of gather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gather_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by gather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gather_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by gather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gather_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by gather collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gather_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by gather collective");			
    /* End: Register PVARs for Gather algorithms */
	
    /* Begin: Register PVARs for Reduce Scatter algorithms */
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_scatter_basic,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times reduce scatter basic was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_scatter_rec_halving,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times recursive halving reduce scatter was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_scatter_pairwise,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times pairwise reduce scatter was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_scatter_non_comm,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times non_comm reduce scatter was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_scatter_noncomm,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages noncomm by reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_scatter_noncomm_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by scatter noncomm algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_scatter_basic_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by scatter basic algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_scatter_rec_halving_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by scatter rec halving algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_scatter_pairwise_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by scatter pairwise algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_scatter_non_comm_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by scatter non comm algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_scatter_noncomm_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by scatter noncomm algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_scatter_basic_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by scatter basic algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_scatter_rec_halving_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by scatter rec halving algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_scatter_pairwise_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by scatter pairwise algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_scatter_non_comm_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by scatter non comm algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_scatter_noncomm_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by scatter noncomm algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_scatter_basic_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by scatter basic algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_scatter_rec_halving_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by scatter rec halving algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_scatter_pairwise_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by scatter pairwise algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_scatter_non_comm_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by scatter non comm algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_scatter_noncomm_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by scatter noncomm algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_scatter_basic_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by scatter basic algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_scatter_rec_halving_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by scatter rec halving algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_scatter_pairwise_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by scatter pairwise algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_scatter_non_comm_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by scatter non comm algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_scatter_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by scatter algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_scatter_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by scatter algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_scatter_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by scatter algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_scatter_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by scatter algorithm of reduce collective");
    /* End: Register PVARs for Reduce Scatter algorithms */
	
    /* Begin: Register PVARs for Scatter algorithms */
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_mcast,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times mcast scatter was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_binomial,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times binomial scatter was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_direct,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times direct scatter was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_direct_blk,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times direct blk scatter was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_two_level_binomial,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times two level binomial scatter was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_two_level_direct,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times two level direct scatter was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_inter,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times Intercommunicator scatter was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_mcast_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of bytes received by mcast algorithm of scatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_mcast_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of bytes sent by mcast algorithm of scatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_binomial_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of bytes received by binomial algorithm of scatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_binomial_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of bytes send by binomial algorithm of scatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_direct_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of bytes received by direct algorithm of scatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_direct_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of bytes send by direct algorithm of scatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_direct_blk_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of bytes received by direct block algorithm of scatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_direct_blk_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of bytes send by direct block algorithm of scatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_two_level_direct_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of bytes received by two level direct algorithm of scatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_two_level_direct_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of bytes send by two level direct algorithm of scatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_two_level_binomial_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of bytes received by two level binomial algorithm of scatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_two_level_binomial_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of bytes send by two level binomial algorithm of scatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_inter_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of bytes received by inter communicator algorithm of scatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_inter_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of bytes send by inter communicator algorithm of scatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_mcast_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of messages received by mcast algorithm of scatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_mcast_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of messages sent by mcast algorithm of scatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_binomial_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of messages received by binomial algorithm of scatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_binomial_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of messages send by binomial algorithm of scatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_direct_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of messages received by direct algorithm of scatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_direct_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of messages send by direct algorithm of scatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_direct_blk_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of messages received by direct block algorithm of scatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_direct_blk_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of messages send by direct block algorithm of scatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_two_level_direct_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of messages received by two level direct algorithm of scatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_two_level_direct_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of messages send by two level direct algorithm of scatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_two_level_binomial_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of messages received by two level binomial algorithm of scatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_two_level_binomial_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of messages send by two level binomial algorithm of scatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_inter_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of messages received by inter communicator algorithm of scatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scatter_inter_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of messages send by inter communicator algorithm of scatter collective");
	MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
			MV2,
			MPI_UNSIGNED_LONG_LONG,
			mv2_coll_scatter_bytes_send,
			MPI_T_VERBOSITY_USER_BASIC,
			MPI_T_BIND_NO_OBJECT,
			(MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
			"COLLECTIVE", /* category name */
			"Number of bytes send by Scatter");
	MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
			MV2,
			MPI_UNSIGNED_LONG_LONG,
			mv2_coll_scatter_bytes_recv,
			MPI_T_VERBOSITY_USER_BASIC,
			MPI_T_BIND_NO_OBJECT,
			(MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
			"COLLECTIVE", /* category name */
			"Number of bytes receive by Scatter");
	MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
			MV2,
			MPI_UNSIGNED_LONG_LONG,
			mv2_coll_scatter_count_send,
			MPI_T_VERBOSITY_USER_BASIC,
			MPI_T_BIND_NO_OBJECT,
			(MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
			"COLLECTIVE", /* category name */
			"Number of times Send was invoked in Scatter");
	MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
			MV2,
			MPI_UNSIGNED_LONG_LONG,
			mv2_coll_scatter_count_recv,
			MPI_T_VERBOSITY_USER_BASIC,
			MPI_T_BIND_NO_OBJECT,
			(MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
			"COLLECTIVE", /* category name */
			"Number of times Recv was invoked in Scatter");

        /* End: Register PVARs for Scatter algorithms */

    /* Begin: Register PVARs for Reduce algorithms */
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_binomial,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times binomial reduce was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_redscat_gather,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times reduce scatter based reduce was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_shmem,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times shmem reduce was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_knomial,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times knomial reduce was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_zcpy,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times zcpy reduce was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_binomial_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by binomial algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_redscat_gather_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by redscat gather algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_two_level_helper_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by two level helper algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_knomial_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by knomial algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_zcpy_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by zcpy algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_binomial_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by binomial algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_redscat_gather_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by redscat gather algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_two_level_helper_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by two level helper algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_knomial_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by knomial algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_zcpy_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by zcpy algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_binomial_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by binomial algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_redscat_gather_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by redscat gather algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_two_level_helper_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by two level helper algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_knomial_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by knomial algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_zcpy_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by zcpy algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_binomial_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by binomial algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_redscat_gather_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by redscat gather algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_two_level_helper_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by two level helper algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_knomial_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by knomial algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_zcpy_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by zcpy algorithm of reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by reduce collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_reduce_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by reduce collective");
    /* End: Register PVARs for Reduce algorithms */
	
    /* Begin: Register PVARs for Gatherv algorithms */
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gatherv_algo,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times gatherv was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gatherv_default_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by default algorithm of gatherv collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gatherv_default_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by default algorithm of gatherv collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gatherv_default_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by default algorithm of gatherv collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gatherv_default_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by default algorithm of gatherv collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gatherv_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by gatherv collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gatherv_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by gatherv collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gatherv_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by gatherv collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_gatherv_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by gatherv collective");			
    /* End: Register PVARs for Gatherv algorithms */
    /* Begin: Register PVARs for Allgatherv algorithms */
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgatherv_rec_doubling,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times recursive doubling Allgatherv was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgatherv_bruck,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times bruck Allgatherv was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgatherv_ring,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times ring Allgatherv was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgatherv_ring_cyclic,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times cyclic ring Allgatherv was invoked");
    /* End: Register PVARs for Allgatherv algorithms */
    /* Beign: Register PVARs for Barrier algorithms */
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_barrier_pairwise,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times pairwise barrier was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_barrier_shmem,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times shmem barrier was invoked");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_barrier_pairwise_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by pairwise algorithm of barrier collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_barrier_pairwise_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by pairwise algorithm of barrier collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_barrier_pairwise_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by pairwise algorithm of barrier collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_barrier_pairwise_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by pairwise algorithm of barrier collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_barrier_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by barrier collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_barrier_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by barrier collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_barrier_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by barrier collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_barrier_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by barrier collective");
	/* End: Register PVARs for Barrier algorithms */
	
    /* Beign: Register PVARs for allgatherv algorithms */
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgatherv_rec_doubling_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by rec doubling algorithm of allgatherv collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgatherv_bruck_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by bruck algorithm of allgatherv collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgatherv_ring_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by ring algorithm of allgatherv collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgatherv_rec_doubling_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by rec doubling algorithm of allgatherv collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgatherv_bruck_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by bruck algorithm of allgatherv collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgatherv_ring_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by ring algorithm of allgatherv collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgatherv_rec_doubling_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by rec doubling algorithm of allgatherv collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgatherv_bruck_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by bruck algorithm of allgatherv collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgatherv_ring_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by ring algorithm of allgatherv collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgatherv_rec_doubling_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by rec doubling algorithm of allgatherv collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgatherv_bruck_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by bruck algorithm of allgatherv collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgatherv_ring_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by ring algorithm of allgatherv collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgatherv_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by allgatherv collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgatherv_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by allgatherv collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgatherv_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by allgatherv collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_allgatherv_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by allgatherv collective");
    /* End: Register PVARs for AllGatherv algorithms */
	
    /* Begin: Register PVARs for Exscan algorithms */
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_exscan_algo,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times exscan was invoked");
    /* End: Register PVARs for Exscan algorithms */
    /* Begin: Register PVARs for Scan algorithms */
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_scan_algo,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE", /* category name */
            "Number of times scan was invoked");
    /* End: Register PVARs for Scan algorithms */
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG,
            mv2_reg_cache_hits,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "Number of registration cache hits");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG,
            mv2_reg_cache_misses,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "Number of registration cache misses");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG,
            mv2_vbuf_allocated,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "Number of VBUFs allocated");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG,
            mv2_vbuf_freed,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "Number of VBUFs freed");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG,
            mv2_ud_vbuf_allocated,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "Number of UD VBUFs allocated");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG,
            mv2_ud_vbuf_freed,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "Number of UD VBUFs freed");
    MPIR_T_PVAR_LEVEL_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG,
            mv2_vbuf_available,
            0, /* initial value */
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "Number of VBUFs available");
    MPIR_T_PVAR_LEVEL_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG,
            mv2_smp_eager_avail_buffer,
            0, /* initial value */
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "Available number of SMP bytes for eager");
    MPIR_T_PVAR_LEVEL_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG,
            mv2_smp_rndv_avail_buffer,
            0, /* initial value */
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "Available number of SMP bytes for rndv");
    MPIR_T_PVAR_LEVEL_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG,
            mv2_smp_eager_total_buffer,
            0, /* initial value */
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "Total number of SMP bytes for eager");
    MPIR_T_PVAR_LEVEL_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG,
            mv2_smp_rndv_total_buffer,
            0, /* initial value */
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "Total number of SMP bytes for rndv");
    MPIR_T_PVAR_LEVEL_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG,
            mv2_smp_eager_buffer_max_use,
            0, /* initial value */
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "Max number of SMP bytes used for eager");
    MPIR_T_PVAR_LEVEL_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG,
            mv2_smp_rndv_buffer_max_use,
            0, /* initial value */
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "Max number of SMP bytes used for rndv");
    MPIR_T_PVAR_LEVEL_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG,
            mv2_ud_vbuf_available,
            0, /* initial value */
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "Number of UD VBUFs available");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG,
            mv2_smp_eager_sent,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "Number of SMP bytes sent through eager");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG,
            mv2_smp_rndv_sent,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "Number of SMP bytes sent through rendezvous");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG,
            mv2_smp_eager_received,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "Number of SMP bytes received through eager");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG,
            mv2_smp_rndv_received,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "Number of SMP bytes received through rendezvous");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG,
            mv2_rdmafp_ctrl_packet_count,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "Number of RDMA FP CTRL Packet count");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG,
            mv2_rdmafp_out_of_order_packet_count,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "Number of RDMA FP Out of Order Packet count");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG,
            mv2_rdmafp_exact_recv_count,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "Number of RDMA FP Exact Recv count");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG,
            mv2_ibv_channel_ctrl_packet_count,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "Number of IB control packets");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG,
            mv2_ibv_channel_out_of_order_packet_count,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "Number of IB out-of-order packets");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG,
            mv2_ibv_channel_exact_recv_count,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "CH3", /* category name */
            "Number of IB exact receives");
    /* End: Register PVARs for IB algorithms */
	
    /* Begin: Register PVARs for iscatter algorithms */
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_iscatter_binomial_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by binomial algorithm of iscatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_iscatter_binomial_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by binomial algorithm of iscatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_iscatter_binomial_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by binomial algorithm of iscatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_iscatter_binomial_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by binomial algorithm of iscatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_iscatter_bytes_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes send by iscatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_iscatter_bytes_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Number of bytes recv by iscatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_iscatter_count_send,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages send by iscatter collective");
    MPIR_T_PVAR_COUNTER_REGISTER_STATIC(
            MV2,
            MPI_UNSIGNED_LONG_LONG,
            mv2_coll_iscatter_count_recv,
            MPI_T_VERBOSITY_USER_BASIC,
            MPI_T_BIND_NO_OBJECT,
            (MPIR_T_PVAR_FLAG_READONLY | MPIR_T_PVAR_FLAG_CONTINUOUS),
            "COLLECTIVE",
            "Count of messages recv by iscatter collective");
	/* End: Register PVARs for iscatter algorithms */
}

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

#ifndef _IBV_MCAST_H_
#define _IBV_MCAST_H_

#include "mv2_ud.h"

#define MCAST_UMAD_SEND_RETRIES (5)
#define MCAST_UMAD_SEND_TIMEOUT (100)   // milli seconds
#define MCAST_UMAD_RECV_TIMEOUT (-1)    // milli seconds
#define MCAST_DEF_TRANS_ID 0x12345678
#define MCAST_DEF_QKEY (0)
#define MCAST_QP (0xFFFFFF)
#define MCAST_MAX_UMAD_RETRIES (8)
#define SA_CLASS_VERSION (2)


#define GID_FMT "%02x:%02x:%02x:%02x:%02x:%02x:%02x:%02x:%02x:%02x:%02x:%02x:%02x:%02x:%02x:%02x\n"

#define PRINT_DEBUG_GID(gid)    \
do {                            \
    PRINT_DEBUG(DEBUG_MCST_verbose > 2, GID_FMT,                                 \
            (gid)->raw[0], (gid)->raw[1], (gid)->raw[2], (gid)->raw[3],         \
            (gid)->raw[4], (gid)->raw[5], (gid)->raw[6], (gid)->raw[7],         \
            (gid)->raw[8], (gid)->raw[9], (gid)->raw[10], (gid)->raw[11],       \
            (gid)->raw[12], (gid)->raw[13], (gid)->raw[14], (gid)->raw[15]);    \
} while (0);

#define IBV_POST_MCAST_SEND(_v, _mcast_ctx)                     \
do {                                                            \
    int ret;                                                    \
    ret = ibv_post_send(_mcast_ctx->ud_ctx->qp,                 \
                &(_v->desc.u.sr), &(_v->desc.y.bad_sr));        \
    if (ret) {                                                  \
        PRINT_ERROR("Failed to post mcast send errno:%d\n",     \
                                errno);                         \
    }                                                           \
    _mcast_ctx->ud_ctx->send_wqes_avail--;                      \
    while (_mcast_ctx->ud_ctx->send_wqes_avail <= 0) {          \
        MPIDI_CH3I_Progress(FALSE, NULL);                       \
    }                                                           \
}while(0)


typedef struct mcast_init_info_t {
    uint8_t gid[16];
} mcast_init_info_t;

typedef enum {
    MCAST_COMM_INACTIVE,
    MCAST_COMM_INITIALIZING,
    MCAST_COMM_ACTIVE,
    MCAST_COMM_FAILED,
} mcast_comm_status_t;

typedef struct mcast_grp_info_t {
    uint16_t status;
    uint16_t comm_id;
    uint16_t mlid;
    union ibv_gid mgid;
} mcast_grp_info_t;

typedef struct mcast_info {
    struct ibv_ah_attr ah_attr;
    struct ibv_ah *ah;
    // TODO: can we handle differently?
    mcast_init_info_t *init_info;
    mcast_grp_info_t grp_info;
} mcast_info_t;


typedef struct bcast_info_t {
    uint8_t in_recv_progress;
    uint32_t win_head;
    uint32_t win_tail;
    long nack_time;
    long resend_time; 
    message_queue_t send_window;
    message_queue_t recv_window;
    mcast_info_t minfo;         /* inter node mcast info */
} bcast_info_t;

typedef struct mcast_init_elem_t {
    mcast_comm_status_t status;
    int comm_id;
    int init_retries;
    int acks_pending;
    long init_timer;
    unsigned char *init_ack_mask;
    struct mcast_init_elem_t *next;
} mcast_init_elem_t;

typedef struct mcast_context_t {
    mcast_init_elem_t *init_list;
    mv2_ud_ctx_t *ud_ctx;
} mcast_context_t;

extern mcast_context_t *mcast_ctx;

/* function return values */
enum {
    MCAST_SUCCESS,
    MCAST_FAILURE
};

#define MAX_MCAST_FRAGMENT_SIZE (MRAIL_MAX_UD_SIZE - sizeof(MPIDI_CH3_Pkt_mcast_t))

#define MV2_DIFF(start, end) ((int)(start) - (end))
#define IS_MCAST_WINDOW_FULL(a, b) \
    ((MV2_DIFF(a, b) >= mcast_window_size -1) ? 1 : 0)
#define MV2_MCAST_MAX_COMMS (4096)

extern MPID_Comm *comm_table[];

int mv2_mcast_init_bcast_info(bcast_info_t ** bcast_info);
int mv2_setup_multicast(mcast_info_t * minfo, MPID_Comm *);
int mv2_cleanup_multicast(mcast_info_t * minfo, MPID_Comm *);
void mv2_process_mcast_msg(vbuf * v);
void mv2_mcast_flush_sendwin(message_queue_t * q);
void mv2_mcast_handle_nack(MPIDI_CH3_Pkt_mcast_nack_t * p);
void mv2_mcast_send(bcast_info_t * bcast_info, char *buf, int len);
void mv2_mcast_recv(bcast_info_t * bcast_info, char *buf, int len, int root);

void mv2_mcast_process_comm_init_req(mcast_init_elem_t * list);
void mv2_mcast_handle_init_ack(MPIDI_CH3_Pkt_mcast_init_ack_t * p);
int mv2_mcast_progress_comm_ready(MPID_Comm * comm_ptr);
MPID_Comm *mv2_mcast_find_comm(int comm_id);
mv2_ud_ctx_t * mv2_mcast_prepare_ud_ctx();

/* bitmask operation on char bit array */
#define MV2_CHAR_ISBITSET(x,i) ((x[(i)>>3] & (1<<((i)&7))) != 0)
#define MV2_CHAR_SETBIT(x,i)   (x[(i)>>3] |= (1<<((i)&7)))
#define MV2_CHAR_CLRBIT(x,i)   (x[(i)>>3] &= (1<<((i)&7))^0xFF);

/* flags for UMAD usage */

typedef enum {
    SUBN_ADM_COMPMASK_MGID      = (1ULL << 0),
    SUBN_ADM_COMPMASK_PORT_GID  = (1ULL << 1),
    SUBN_ADM_COMPMASK_QKEY      = (1ULL << 2),
    SUBN_ADM_COMPMASK_PKEY      = (1ULL << 7),
    SUBN_ADM_COMPMASK_TCLASS    = (1ULL << 6),
    SUBN_ADM_COMPMASK_SL        = (1ULL << 12),
    SUBN_ADM_COMPMASK_FLOW      = (1ULL << 13),
    SUBN_ADM_COMPMASK_JOIN_STATE = (1ULL << 16),
} subn_adm_component_mask;


#endif /* _IBV_MCAST_H_ */

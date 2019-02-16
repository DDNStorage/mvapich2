/*
 * Copyright (c) 2001-2018, The Ohio State University. All rights
 * reserved.
 *
 * This file is part of the MVAPICH2 software package developed by the
 * team members of The Ohio State University's Network-Based Computing
 * Laboratory (NBCL), headed by Professor Dhabaleswar K. (DK) Panda.
 *
 * For detailed copyright and licensing information, please refer to the
 * copyright file COPYRIGHT in the top level MVAPICH2 directory.
 */
#include "upmi.h"
#include <stdlib.h>
#include <pthread.h>
#include <mpimem.h>
#include <unistd.h>
#include <sys/types.h>

struct PMI_keyval_t;
int _size, _rank, _appnum;
int _singleton_mode = 0;
static int _in_ibarrier = 0;
static int _in_iallgather = 0;
#if defined(HAVE_PMI2_IALLGATHER) && defined(HAVE_PMI2_IALLGATHER_WAIT)
static void * _iallgather_data = NULL;
static size_t _iallgather_data_size = 0;
#endif
pthread_mutex_t upmi_lock;

void UPMI_lock_init(void) {
    pthread_mutex_init(&upmi_lock, NULL);
}

void UPMI_lock_destroy(void) {
    pthread_mutex_destroy(&upmi_lock);
}

void UPMI_lock(void) {
    pthread_mutex_lock(&upmi_lock);
}

void UPMI_unlock(void) {
    pthread_mutex_unlock(&upmi_lock);
}

int UPMI_INIT( int *spawned ) {
    int pmi_ret_val;
    #ifdef USE_PMI2_API
    pmi_ret_val = PMI2_Init( spawned, &_size, &_rank, &_appnum );
    #else
    UPMI_lock_init();
    pmi_ret_val = PMI_Init( spawned );
    #endif
    return pmi_ret_val;
}

int UPMI_INITIALIZED( int *initialized ) { 
    #ifdef USE_PMI2_API
    *initialized = PMI2_Initialized();
    return UPMI_SUCCESS;
    #else
    return PMI_Initialized( initialized );
    #endif
}

int UPMI_FINALIZE( void ) { 
    int pmi_ret_val;
    #ifdef USE_PMI2_API
    pmi_ret_val = PMI2_Finalize();
    #else
    UPMI_lock();
    pmi_ret_val = PMI_Finalize();
    UPMI_unlock();
    UPMI_lock_destroy();
    #endif
    return pmi_ret_val;
}

int UPMI_GET_SIZE( int *size ) { 
    #ifdef USE_PMI2_API
    *size = _size;
    return UPMI_SUCCESS;
    #else
    return PMI_Get_size( size );
    #endif
}

int UPMI_GET_RANK( int *rank ) { 
    #ifdef USE_PMI2_API
    *rank = _rank;
    return UPMI_SUCCESS;
    #else
    return PMI_Get_rank( rank );
    #endif
}

int UPMI_GET_APPNUM( int *appnum ) { 
    int pmi_ret_val;
    #ifdef USE_PMI2_API
    *appnum = _appnum;
    pmi_ret_val = UPMI_SUCCESS;
    #else
    UPMI_lock();
    pmi_ret_val = PMI_Get_appnum( appnum );
    UPMI_unlock();
    #endif
    return pmi_ret_val;
}

int UPMI_GET_UNIVERSE_SIZE( int *size ) { 
    int pmi_ret_val;
    #ifdef USE_PMI2_API
    char name[] = "universeSize";
    int outlen, found;
    PMI2_Info_GetJobAttrIntArray( name, size, sizeof (int), &outlen, &found );
    if( found && outlen==1 ) {
        pmi_ret_val = UPMI_SUCCESS;
    } else {
        pmi_ret_val = UPMI_FAIL;
    }
    #else
    UPMI_lock();
    pmi_ret_val = PMI_Get_universe_size( size );
    UPMI_unlock();
    #endif
    return pmi_ret_val;
}

int UPMI_BARRIER( void ) { 
    int pmi_ret_val;
    #ifdef USE_PMI2_API
    pmi_ret_val = PMI2_KVS_Fence();
    #else
    UPMI_lock();
    pmi_ret_val = PMI_Barrier();
    UPMI_unlock();
    #endif
    return pmi_ret_val;
}

int UPMI_IBARRIER( void ) { 
    int rc;

    UPMI_lock();
    if (!_in_ibarrier) {
        _in_ibarrier = 1;
#ifdef USE_PMI2_API
#   if defined(HAVE_PMI2_KVS_IFENCE) && defined(HAVE_PMI2_KVS_WAIT)
        rc = PMI2_KVS_Ifence();
#   else
        rc = PMI2_KVS_Fence();
#   endif
#else
#   if defined(HAVE_PMI_IBARRIER) && defined(HAVE_PMI_WAIT)
        rc = PMI_Ibarrier();
#   else
        rc = PMI_Barrier();
#   endif
#endif
    } else {
        rc = UPMI_FAIL;
    }
    UPMI_unlock();

    return rc;
}

int UPMI_WAIT( void ) { 
    int rc = UPMI_SUCCESS;

    UPMI_lock();
    if (_in_ibarrier) {
        _in_ibarrier = 0;
#ifdef USE_PMI2_API
#   if defined(HAVE_PMI2_KVS_IFENCE) && defined(HAVE_PMI2_KVS_WAIT)
        rc = PMI2_KVS_Wait();
#   endif
#else
#   if defined(HAVE_PMI_IBARRIER) && defined(HAVE_PMI_WAIT)
        rc = PMI_Wait();
#   endif
#endif
    } else {
        rc = UPMI_SUCCESS;
    }
    UPMI_unlock();

    return rc;
}

int UPMI_IALLGATHER( const char value[] ) {
    int rc;

    UPMI_lock();
    if (!_in_iallgather) {
        _in_iallgather = 1;
#ifdef USE_PMI2_API
#   if defined(HAVE_PMI2_SHMEM_IALLGATHER) \
        && defined(HAVE_PMI2_SHMEM_IALLGATHER_WAIT)
        rc = PMI2_SHMEM_Iallgather(value);
#   elif defined(HAVE_PMI2_IALLGATHER) && defined(HAVE_PMI2_IALLGATHER_WAIT)
        rc = PMI2_Iallgather(value);
        if (UPMI_SUCCESS == rc) {
            if (NULL == _iallgather_data) {
                _iallgather_data_size = _size * PMI2_MAX_VALLEN * sizeof(char);
                _iallgather_data = MPIU_Malloc(_iallgather_data_size);
            }

            if (NULL == _iallgather_data) {
                rc = UPMI_FAIL;
            } else {
                memset(_iallgather_data, 0, _iallgather_data_size);
            }
        }
#   else
        rc = UPMI_FAIL;
#   endif
#endif
    } else {
        rc = UPMI_FAIL;
    }
    UPMI_unlock();

    return rc;
}

int UPMI_IALLGATHER_WAIT( void **buf ) {
    int rc = UPMI_SUCCESS;

    UPMI_lock();
    if (_in_iallgather) {
        _in_iallgather = 0;
#ifdef USE_PMI2_API
#   if defined(HAVE_PMI2_SHMEM_IALLGATHER) \
        && defined(HAVE_PMI2_SHMEM_IALLGATHER_WAIT)
        rc = PMI2_SHMEM_Iallgather_wait(buf);
#   elif defined(HAVE_PMI2_IALLGATHER) && defined(HAVE_PMI2_IALLGATHER_WAIT)
        *buf = _iallgather_data;
        if (NULL != *buf) {
            rc = PMI2_Iallgather_wait(*buf);
        } else {
            rc = UPMI_FAIL;
        }
#   endif
#else
        rc = UPMI_FAIL;
#endif
    } else {
        rc = UPMI_SUCCESS;
    }
    UPMI_unlock();

    return rc;
}

int UPMI_IALLGATHER_FREE( void ) {
    int rc = UPMI_SUCCESS;

    UPMI_lock();
    if (!_in_iallgather) {
#ifdef USE_PMI2_API
#   if defined(HAVE_PMI2_SHMEM_IALLGATHER) \
        && defined(HAVE_PMI2_SHMEM_IALLGATHER_WAIT)
        /* nothing to do */
#   elif defined(HAVE_PMI2_IALLGATHER) && defined(HAVE_PMI2_IALLGATHER_WAIT)
        if (NULL != _iallgather_data) {
            MPIU_Free(_iallgather_data);
        }
#   endif
#endif
    } else {
        rc = UPMI_FAIL;
    }
    UPMI_unlock();

    return rc;
}

int UPMI_ABORT( int exit_code, const char error_msg[] ) { 
    #ifdef USE_PMI2_API
    return PMI2_Abort( 1, error_msg );    //flag = 1, abort all processes
    #else
    return PMI_Abort( exit_code, error_msg );
    #endif
}

int UPMI_KVS_GET_KEY_LENGTH_MAX( int *length ) { 
    #ifdef USE_PMI2_API
    *length = PMI2_MAX_KEYLEN;
    return UPMI_SUCCESS;
    #else
    return PMI_KVS_Get_key_length_max( length );
    #endif
}

int UPMI_KVS_GET_NAME_LENGTH_MAX( int *length ) { 
    #ifdef USE_PMI2_API
    *length = PMI2_MAX_KEYLEN; //TODO is this correct?
    return UPMI_SUCCESS;
    #else
    return PMI_KVS_Get_name_length_max( length );
    #endif
}

int UPMI_KVS_GET_VALUE_LENGTH_MAX( int *length ) { 
    #ifdef USE_PMI2_API
    *length = PMI2_MAX_VALLEN;
    return UPMI_SUCCESS;
    #else
    return PMI_KVS_Get_value_length_max( length );
    #endif
}

int UPMI_KVS_GET_MY_NAME( char kvsname[], int length ) {
    int pmi_ret_val;
    #ifdef USE_PMI2_API
    pmi_ret_val = PMI2_Job_GetId( kvsname, length );
    if (pmi_ret_val == PMI2_ERR_OTHER && _size == 1) {
        _singleton_mode = 1;
        sprintf(kvsname, "singleton_kvs_%llu", getpid());
        pmi_ret_val = UPMI_SUCCESS;
    }
    #else
    UPMI_lock();
    pmi_ret_val = PMI_KVS_Get_my_name( kvsname, length );
    UPMI_unlock();
    #endif
    return pmi_ret_val;
}

int UPMI_KVS_PUT( const char kvsname[], const char key[], const char value[] ) { 
    int pmi_ret_val;
    #ifdef USE_PMI2_API
    pmi_ret_val = PMI2_KVS_Put( key, value );
    #else
    UPMI_lock();
    pmi_ret_val = PMI_KVS_Put( kvsname, key, value );
    UPMI_unlock();
    #endif
    return pmi_ret_val;
}

int UPMI_KVS_GET( const char kvsname[], const char key[], char value[], int length ) { 
    int pmi_ret_val;
    #ifdef USE_PMI2_API
    int vallen;
    pmi_ret_val = PMI2_KVS_Get( kvsname, PMI2_ID_NULL, key, value, length, &vallen );
    #else
    UPMI_lock();
    pmi_ret_val = PMI_KVS_Get( kvsname, key, value, length );
    UPMI_unlock();
    #endif
    return pmi_ret_val;
}

int UPMI_KVS_COMMIT( const char kvsname[] ) { 
    #ifdef USE_PMI2_API
    //return PMI2_KVS_Fence();
    return UPMI_SUCCESS;
    #else
    return PMI_KVS_Commit( kvsname );
    #endif
}

int UPMI_PUBLISH_NAME( const char service_name[], const char port[], const struct MPID_Info *info_ptr ) { 
    int pmi_ret_val;
    #ifdef USE_PMI2_API
    pmi_ret_val = PMI2_Nameserv_publish( service_name, info_ptr, port );
    #else
    UPMI_lock();
    pmi_ret_val = PMI_Publish_name( service_name, port );
    UPMI_unlock();
    #endif
    return pmi_ret_val;
}

int UPMI_UNPUBLISH_NAME( const char service_name[], const struct MPID_Info *info_ptr ) { 
    int pmi_ret_val;
    #ifdef USE_PMI2_API
    pmi_ret_val = PMI2_Nameserv_unpublish( service_name, info_ptr );
    #else
    UPMI_lock();
    pmi_ret_val = PMI_Unpublish_name( service_name );
    UPMI_unlock();
    #endif
    return pmi_ret_val;
}

int UPMI_LOOKUP_NAME( const char service_name[], char port[], const struct MPID_Info *info_ptr ) { 
    int pmi_ret_val;
    #ifdef USE_PMI2_API
    pmi_ret_val = PMI2_Nameserv_lookup( service_name, info_ptr, port, sizeof port );  
    #else
    UPMI_lock();
    pmi_ret_val = PMI_Lookup_name( service_name, port );
    UPMI_unlock();
    #endif
    return pmi_ret_val;
}

int UPMI_GET_NODE_ATTR( const char name[], char value[], int valuelen, int *found, int waitfor ) {
    #ifdef USE_PMI2_API
    return PMI2_Info_GetNodeAttr( name, value, valuelen, found, waitfor );
    #else
    return UPMI_FAIL;
    #endif
}

int UPMI_GET_NODE_ATTR_INT_ARRAY( const char name[], int array[], int arraylen, int *outlen, int *found ) {
    #ifdef USE_PMI2_API
    return PMI2_Info_GetNodeAttrIntArray( name, array, arraylen, outlen, found );
    #else
    return UPMI_FAIL;
    #endif
}

int UPMI_PUT_NODE_ATTR( const char name[], const char value[] ) {
    #ifdef USE_PMI2_API
    return PMI2_Info_PutNodeAttr( name, value );
    #else
    return UPMI_FAIL;
    #endif
}

int UPMI_GET_JOB_ATTR( const char name[], char value[], int valuelen, int *found ) {
    #ifdef USE_PMI2_API
    return PMI2_Info_GetJobAttr( name, value, valuelen, found );
    #else
    return UPMI_FAIL;
    #endif
}

int UPMI_GET_JOB_ATTR_INT_ARRAY( const char name[], int array[], int arraylen, int *outlen, int *found ) {
    #ifdef USE_PMI2_API
    return PMI2_Info_GetJobAttrIntArray( name, array, arraylen, outlen, found );
    #else
    return UPMI_FAIL;
    #endif
}

int UPMI_JOB_SPAWN(int count,
                   const char * cmds[],
                   int argcs[],
                   const char ** argvs[],
                   const int maxprocs[],
                   const int info_keyval_sizes[],
                   const void *info_keyval_vectors[],
                   int preput_keyval_size,
                   const void *preput_keyval_vector,
                   char jobId[],
                   int jobIdSize,
                   int errors[])
{
    int pmi_ret_val;
    #ifdef USE_PMI2_API
    pmi_ret_val = PMI2_Job_Spawn( count, cmds, argcs, argvs, maxprocs,
                           info_keyval_sizes, (const struct MPID_Info**)info_keyval_vectors,
                           preput_keyval_size, (const struct MPID_Info**)preput_keyval_vector,
                           jobId, jobIdSize, errors );
    #else
    UPMI_lock();
    pmi_ret_val = PMI_Spawn_multiple( count, cmds, argvs, maxprocs,
                               info_keyval_sizes, (const struct PMI_keyval_t**)info_keyval_vectors,
                               preput_keyval_size, (const struct PMI_keyval_t*)preput_keyval_vector,
                               errors );
    UPMI_unlock();
    #endif
    return pmi_ret_val;
}


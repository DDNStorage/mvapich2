!   -*- Mode: Fortran; -*-
!
!   (C) 2014 by Argonne National Laboratory.
!   See COPYRIGHT in top-level directory.
!
subroutine PMPIR_Scatterv_f08ts(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, &
    recvtype, root, comm, ierror)
    use, intrinsic :: iso_c_binding, only : c_int
    use :: mpi_f08, only : MPI_Datatype, MPI_Comm
    use :: mpi_c_interface, only : c_Datatype, c_Comm
    use :: mpi_c_interface, only : MPIR_Scatterv_cdesc, MPIR_Comm_size_c

    implicit none

    type(*), dimension(..), intent(in) :: sendbuf
    type(*), dimension(..) :: recvbuf
    integer, intent(in) :: recvcount
    integer, intent(in) :: root
    integer, intent(in) :: sendcounts(*)
    integer, intent(in) :: displs(*)
    type(MPI_Datatype), intent(in) :: sendtype
    type(MPI_Datatype), intent(in) :: recvtype
    type(MPI_Comm), intent(in) :: comm
    integer, optional, intent(out) :: ierror

    integer(c_int) :: recvcount_c
    integer(c_int) :: root_c
    integer(c_int), allocatable :: sendcounts_c(:)
    integer(c_int), allocatable :: displs_c(:)
    integer(c_Datatype) :: sendtype_c
    integer(c_Datatype) :: recvtype_c
    integer(c_Comm) :: comm_c
    integer(c_int) :: ierror_c
    integer(c_int) :: err, length ! To get length of assumed-size arrays

    if (c_int == kind(0)) then
        ierror_c = MPIR_Scatterv_cdesc(sendbuf, sendcounts, displs, sendtype%MPI_VAL, &
            recvbuf, recvcount, recvtype%MPI_VAL, root, comm%MPI_VAL)
    else
        sendtype_c = sendtype%MPI_VAL
        recvcount_c = recvcount
        recvtype_c = recvtype%MPI_VAL
        root_c = root
        comm_c = comm%MPI_VAL

        err = MPIR_Comm_size_c(comm_c, length)
        sendcounts_c = sendcounts(1:length)
        displs_c = displs(1:length)

        ierror_c = MPIR_Scatterv_cdesc(sendbuf, sendcounts_c, displs_c, sendtype_c, recvbuf, recvcount_c, recvtype_c, &
            root_c, comm_c)
    end if

    if(present(ierror)) ierror = ierror_c

end subroutine PMPIR_Scatterv_f08ts

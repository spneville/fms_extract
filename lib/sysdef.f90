  module sysdef

    implicit none

    save

!-----------------------------------------------------------------------
! Parameters common to all trajectories
!-----------------------------------------------------------------------
    integer                                     :: natm,nsta,nintraj,&
                                                   nstep
    real*8, dimension(:), allocatable           :: alpha,atnum,atmass,&
                                                   r0
    real*8                                      :: dt,tf
    character(len=2), dimension(:), allocatable :: atlbl
    logical(kind=4)                             :: lmomrep

  end module sysdef

  module sysdef

    implicit none

    save

!-----------------------------------------------------------------------
! Parameters common to all trajectories
!-----------------------------------------------------------------------
    integer*8                                   :: natm,nsta,nintraj,&
                                                   nstep
    real*8, dimension(:), allocatable           :: alpha,atnum,atmass
    real*8                                      :: dt,tf
    character(len=2), dimension(:), allocatable :: atlbl

  end module sysdef

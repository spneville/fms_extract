  module trajdef

    implicit none

    save

!#######################################################################
! trajectory: derived type storing all the parameters for all
!             trajectories forming the wavepacket
!#######################################################################

    type trajectory
       integer*8                               :: ntraj
       integer*8, dimension(:), allocatable    :: ista,tkill,ispawn,&
                                                  tspawn
       real*8, dimension(:,:,:), allocatable   :: r,p
       real*8, dimension(:,:), allocatable     :: phase
       complex*16, dimension(:,:), allocatable :: coe
    end type trajectory

    type(trajectory), dimension(:), allocatable :: traj

!#######################################################################

  end module trajdef

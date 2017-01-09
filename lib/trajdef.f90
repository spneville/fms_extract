  module trajdef

    implicit none

    save

!#######################################################################
! trajectory: derived type storing all the parameters for all
!             trajectories forming the wavepacket, and the electronic
!             structure information at the centres of the Gaussian
!             basis functions
!#######################################################################

    type trajectory
       integer                                   :: ntraj
       integer, dimension(:), allocatable        :: ista,tkill,ispawn,&
                                                    tspawn
       real*8, dimension(:,:,:), allocatable     :: r,p
       real*8, dimension(:,:), allocatable       :: phase
       real*8, dimension(:,:,:), allocatable     :: ener
       real*8, dimension(:,:,:,:), allocatable   :: nact
       real*8, dimension(:,:,:,:,:), allocatable :: dipole
       complex*16, dimension(:,:), allocatable   :: coe
    end type trajectory

    type(trajectory), dimension(:), allocatable :: traj

!#######################################################################

  end module trajdef

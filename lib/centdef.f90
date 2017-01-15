  module centdef

    implicit none

    save

    integer                            :: maxcent
    integer, dimension(:), allocatable :: ncent
    real*8, parameter                  :: othrsh=1e-5
    
!#######################################################################
! centroid: derived type storing all the parameters for all centroids
!           formed from the trajectories of a given IFG
!#######################################################################

    type centroid
       integer                               :: indx1,indx2
       integer                               :: tspawn,tkill
       real*8, dimension(:,:), allocatable   :: r,p
       real*8, dimension(:), allocatable     :: alpha
       complex*16, dimension(:), allocatable :: ovrlp
       logical                               :: lcontrib
    end type centroid

    type(centroid), dimension(:,:), allocatable :: cent

  contains

!#######################################################################
! centindx: returns the unique index for a given pair of trajectory
!           numbers n1 and n2 from the IFG indexed ifg.
!
! For the lower-triangle of trajectory indices,
!
! 1,1 
! 2,1 2,2
! 3,1 3,3 3,3
!  .   .   .   .
!  .   .   .       .
!  .   .   .           .
!
! we define the centroid index for n1,n2 as the number of unique,
! non-equal pairs of indices up to and including the current pair,
! with the lower-triangle of trajectory indices being traversed in
! COLUMN-MAJOR order.
!#######################################################################
    
    function centindx(n1,n2,ifg)

      use trajdef
      
      implicit none

      integer :: centindx,n1,n2,ifg,rindx,cindx,ntraj,nprev,i

      ! No. trajectories in the given IFG
      ntraj=traj(ifg)%ntraj

      ! Row index
      rindx=max(n1,n2)

      ! Column index
      cindx=min(n1,n2)

      ! Number of non-equal pairs in the columns before
      ! the current column
      nprev=0
      do i=1,cindx-1
         nprev=nprev+ntraj-i
      enddo
      
      ! Index
      centindx=nprev+rindx-cindx
      
      return
      
    end function centindx
      
!#######################################################################
    
  end module centdef

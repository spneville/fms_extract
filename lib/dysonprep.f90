  module dysonprep

    implicit none

  contains

!#######################################################################

    subroutine mkdysoninp

      use sysdef
      use trajdef
      use dysonmod
      use expec, only: dstep
      use postprep
      
      implicit none

      integer            :: i,k,n
      real*8, parameter  :: coethrsh=0.05d0
      real*8             :: coesq
      complex*16         :: coe
      character(len=120) :: amain,asub

!-----------------------------------------------------------------------
! Create main directories: one for each pair of IFGs and trajectories
!-----------------------------------------------------------------------
      ! Loop over IFGs
      do i=1,nintraj
         ! Loop over trajectories for the current IFG
         do k=1,traj(i)%ntraj
            ! Get main directory name
            call getmaindir(amain,i,k)
            ! Create main directory
            call makedir(amain)
            ! Write the state number of the current trajectory to file
            call wrstatenumber(amain,i,k)
            ! Write the IFG number of the current trajectory to file
            call wrifgnumber(amain,i)
         enddo
      enddo

!-----------------------------------------------------------------------
! For each IFG/trajectory pair, loop over timesteps and create the
! Columbus and superdyson preparation inputs
!
! For a given trajectory, if t>t_kill, then we do not create the
! directory
!-----------------------------------------------------------------------
      write(6,'(/,70a)') ('-',i=1,70)
      write(6,'(25x,a)') 'Creating Input Files'
      write(6,'(70a,/)') ('-',i=1,70)
      ! Loop over IFGs
      do i=1,nintraj
         ! Loop over trajectories for the current IFG
         do k=1,traj(i)%ntraj

            write(6,'(2x,a3,x,i3,x,a12,i3)') 'IFG',i,', Trajectory',k

            ! Loop over timesteps
            do n=1,nstep,dstep
               
               ! Only proceed if the current trajectory has not yet
               ! been killed
               if (n.ge.traj(i)%tkill(k)) cycle

               ! Create subdirectory for the current timestep
               call getsubdir(asub,i,k,n)
               call makedir(asub)

               ! Write the coefficient of current trajectory to file
               call wrcoeff(asub,i,k,n)
               
               ! Copy the Columbus and superdyson preparation
               ! directories to the current subdirectory
               call copyinpdirs(asub)

               ! Write the current Cartesian coordinates to file
               call copygeom(asub,i,k,n)
               
            enddo
         enddo
      enddo

      return

    end subroutine mkdysoninp

!#######################################################################

  end module dysonprep

  module dipoleprep

    implicit none

  contains

!#######################################################################

    subroutine prep_dipole_inp

      implicit none

!-----------------------------------------------------------------------
! (1) Trajectories
!-----------------------------------------------------------------------
      call wrinp_traj

!-----------------------------------------------------------------------
! (2) Centroids
!-----------------------------------------------------------------------
      call wrinp_cent
      
      return
      
    end subroutine prep_dipole_inp

!#######################################################################

    subroutine wrinp_traj

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
! Output what we are doing
!-----------------------------------------------------------------------
      write(6,'(/,70a)') ('-',i=1,70)
      write(6,'(25x,a)') 'Creating Input Files'
      write(6,'(70a,/)') ('-',i=1,70)
      
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
! Columbus input
!
! For a given trajectory, if t>t_kill, then we do not create the
! directory
!-----------------------------------------------------------------------
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
               
               ! Make the columbus directory
               call system('mkdir '//trim(asub)//'/columbus')

               ! Write the current Cartesian coordinates to file
               ! N.B., we use the copygeom_gamess subroutine as
               !       this just copies across the geom file to
               !       the columbus directory. The name of this
               !       subroutine is rather shit/misleading...
               call copygeom_gamess(asub,i,k,n)
               
            enddo
         enddo
      enddo
      
      return
      
    end subroutine wrinp_traj

!#######################################################################

    subroutine wrinp_cent

      use sysdef
      use trajdef
      use centdef
      use dysonmod
      use expec, only: dstep
      use postprep
      
      implicit none

      integer            :: i,k,n,ifg,n1,n2,indx,istep
      real*8, parameter  :: coethrsh=0.05d0
      real*8             :: coesq
      complex*16         :: coe
      character(len=120) :: amain,asub

!-----------------------------------------------------------------------
! Output what we are doing
!-----------------------------------------------------------------------
      write(6,'(/,70a)') ('-',i=1,70)
      write(6,'(25x,a)') 'Creating Input Files'
      write(6,'(70a,/)') ('-',i=1,70)
      
!-----------------------------------------------------------------------
! Create the main directories: one for each IFG/centroid pair
!-----------------------------------------------------------------------
      ! Loop over the IFGs
      do ifg=1,nintraj

         ! Loop over the 1st trajectory index
         do n1=1,traj(ifg)%ntraj-1

            ! Loop over the 2nd trajectory index
            do n2=n1+1,traj(ifg)%ntraj

               ! Get the current centroid index
               indx=centindx(n1,n2,ifg)

               ! Cycle if current centroid is inactive
               if (.not.cent(ifg,indx)%lcontrib) cycle
               
               ! Get the main directory name
               call getmaindir_cent(amain,ifg,n1,n2)

               ! Create main directory
               call makedir(amain)

               ! Write the state numbers of the trajectories to file
               call wrstatenumber_cent(amain,ifg,n1,n2)

               ! Write the IFG number of the current trajectory to file
               call wrifgnumber(amain,ifg)
               
            enddo
         enddo
         
      enddo
         
!-----------------------------------------------------------------------
! For each IFG/centroid pair, loop over timesteps and create the
! Columbus input
!
! For a given centroid, we only create a directory if both associated
! trajectories are alive
!-----------------------------------------------------------------------
      ! Loop over IFGs
      do ifg=1,nintraj

         ! Loop over the 1st trajectory index
         do n1=1,traj(ifg)%ntraj-1

            ! Loop over the 2nd trajectory index
            do n2=n1+1,traj(ifg)%ntraj

               ! Get the current centroid index
               indx=centindx(n1,n2,ifg)

               ! Cycle if current centroid is inactive
               if (.not.cent(ifg,indx)%lcontrib) cycle

               ! Ouput progress
               write(6,'(2x,a3,x,i3,x,a10,i4)') 'IFG',ifg,', Centroid',indx
               
               ! Loop over the timesteps for which both trajectories are
               ! alive
               do istep=1,nstep,dstep

                  ! Cycle if the two trajectories are not both alive
                  if (istep.lt.cent(ifg,indx)%tspawn &
                       .or.istep.ge.cent(ifg,indx)%tkill) cycle
                  
                  ! Cycle if the overlap between the trajectories is below
                  ! threshold
                  if (abs(cent(ifg,indx)%ovrlp(istep)).lt.othrsh) cycle
                  
                  ! Create the subdirectory for the current timestep
                  call getsubdir_cent(asub,ifg,n1,n2,istep)
                  call makedir(asub)

                  ! Write the coefficient of current trajectories to file
                  call wrcoeff_cent(asub,ifg,n1,n2,istep)

                  ! Make the columbus directory
                  call system('mkdir '//trim(asub)//'/columbus')

                  ! Write the current Cartesian coordinates to file
                  call copygeom_cent(asub,ifg,indx,istep)

               enddo

            enddo

         enddo

      enddo
      
      return
      
    end subroutine wrinp_cent
      
!#######################################################################
    
  end module dipoleprep

  module adcprep

  contains

!#######################################################################

    subroutine mkadcinp
      
      use sysdef
      use trajdef
      use expec
      use postprep
      
      implicit none
      
      integer*8          :: i,k,n
      real*8, parameter  :: coethrsh=0.05d0
      real*8             :: coesq,talive
      complex*16         :: coe
      character(len=120) :: amain,asub,path

!-----------------------------------------------------------------------
! Read in the ADC input templates
!-----------------------------------------------------------------------
      call rdtemplates

!-----------------------------------------------------------------------
! Create main directories: one for each pair of IFGs and trajectories
!-----------------------------------------------------------------------
      ! Loop over IFGs
      do i=1,nintraj
         ! Loop over trajectories for the current IFG
         do k=1,traj(i)%ntraj           
            ! Skip the current trajectory if it did not live
            ! long enough
            talive=traj(i)%tkill(k)-traj(i)%tspawn(k)
            if (talive.lt.thrsh_alive) cycle
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
      write(6,'(70a)') ('-',i=1,70)
      
      ! Loop over IFGs
      do i=1,nintraj
         ! Loop over trajectories for the current IFG
         do k=1,traj(i)%ntraj

            ! Skip the current trajectory if it did not live
            ! long enough
            talive=traj(i)%tkill(k)-traj(i)%tspawn(k)
            if (talive.lt.thrsh_alive) cycle

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

               ! Write the state index of the current trajectory to
               ! current subdirectory
               call wrstatenumber(asub,i,k)

               ! Copy the Columbus and superdyson preparation
               ! directories to the current subdirectory
               call copyinpdirs_adc(asub)

               ! Write the current Cartesian coordinates to the
               ! columbus geom file               
               call copygeom_adc(asub,i,k,n)

               ! Write the ADC input files for the current subdirectory
               call wradcinp(asub,i,k,n)
               
            enddo
         enddo
      enddo
      
      return
      
    end subroutine mkadcinp

!#######################################################################

    subroutine rdtemplates

      use iomod
      use expec
        
      implicit none

      integer*8          :: unit,i,n
      character(len=120) :: string

      ! Loop over the different template files
      do n=1,3
      
         ! Open the template file
         call getfreeunit(unit)
         open(unit,file=adcfile(n),form='formatted',status='old')
      
         ! Determine the no.lines in the template file
         nlines(n)=0
5        read(unit,'(a)',end=10) string
         nlines(n)=nlines(n)+1
         goto 5
         
10       continue
         
         ! Read template file
         rewind(unit)
         do i=1,nlines(n)
            read(unit,'(a)') adcinp(n,i)
         enddo
         
         ! Close the template file
         close(unit)
         
      enddo
      
      return
        
    end subroutine rdtemplates
      
!#######################################################################
    
  end module adcprep

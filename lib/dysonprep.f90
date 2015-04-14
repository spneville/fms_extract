  module dysonprep

    implicit none

  contains

!#######################################################################

    subroutine mkdysoninp

      use sysdef
      use trajdef
      use dysonmod
      use expec, only: dstep

      implicit none

      integer*8          :: i,k,n
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

    subroutine getmaindir(amain,i,k)

      implicit none

      integer*8          :: i,k
      character(len=120) :: amain

      amain=''

      if (i.lt.10) then
         write(amain(1:5),'(a4,i1)') 'ifg0',i
      else
         write(amain(1:5),'(a3,i2)') 'ifg',i
      endif

      if (k.lt.10) then
         write(amain(6:12),'(a6,i1)') '_traj0',k
      else
         write(amain(6:12),'(a5,i2)') '_traj',k
      endif

      return

    end subroutine getmaindir

!#######################################################################

    subroutine getsubdir(asub,i,k,n)

      implicit none

      integer*8          :: i,k,n
      character(len=120) :: asub

      asub=''

      if (i.lt.10) then
         write(asub(1:5),'(a4,i1)') 'ifg0',i
      else
         write(asub(1:5),'(a3,i2)') 'ifg',i
      endif

      if (k.lt.10) then
         write(asub(6:12),'(a6,i1)') '_traj0',k
      else
         write(asub(6:12),'(a5,i2)') '_traj',k
      endif

      if (n.lt.10) then
         write(asub(13:23),'(a9,i1)') '/step0000',n
      else if (n.lt.100) then
         write(asub(13:23),'(a8,i2)') '/step000',n
      else if (n.lt.1000) then
         write(asub(13:23),'(a7,i3)') '/step00',n
      else if (n.lt.10000) then
         write(asub(13:23),'(a6,i4)') '/step0',n
      else
         write(asub(13:23),'(a5,i5)') '/step',n
      endif

      return

    end subroutine getsubdir

!#######################################################################

    subroutine makedir(string)

      implicit none

      character(len=*)   :: string
      character(len=120) :: acmnd
      logical(kind=4)    :: ldir
      
      inquire(file=trim(string),exist=ldir)

      if (ldir) then
         acmnd='rm -rf '//trim(string)//'/*'
      else
         acmnd='mkdir '//trim(string)
      endif
      
      call system(acmnd)

      return

    end subroutine makedir

!#######################################################################

    subroutine wrstatenumber(amain,ifg,itraj)

      use trajdef

      implicit none

      integer*8          :: ifg,itraj,unit
      character(len=120) :: amain
      character(len=130) :: aout

!-----------------------------------------------------------------------
! Open file
!-----------------------------------------------------------------------
      unit=20
      aout=trim(amain)//'/state_id'
      open(unit,file=aout,form='formatted',status='unknown')

!-----------------------------------------------------------------------
! Write state number to file
!-----------------------------------------------------------------------
      write(unit,'(i2)') traj(ifg)%ista(itraj)

!-----------------------------------------------------------------------
! Close file
!-----------------------------------------------------------------------
      close(unit)

      return
      
    end subroutine wrstatenumber

!#######################################################################

    subroutine wrifgnumber(amain,ifg)

      use trajdef

      implicit none

      integer*8          :: ifg,itraj,unit
      character(len=120) :: amain
      character(len=130) :: aout

!-----------------------------------------------------------------------
! Open file
!-----------------------------------------------------------------------
      unit=20
      aout=trim(amain)//'/ifg_number'
      open(unit,file=aout,form='formatted',status='unknown')

!-----------------------------------------------------------------------
! Write state number to file
!-----------------------------------------------------------------------
      write(unit,'(i2)') ifg

!-----------------------------------------------------------------------
! Close file
!-----------------------------------------------------------------------
      close(unit)

      return

    end subroutine wrifgnumber

!#######################################################################

    subroutine wrcoeff(asub,ifg,itraj,istep)

      use trajdef

      implicit none

      integer*8          :: ifg,itraj,istep,unit
      real*8             :: coe_r,coe_i
      character(len=120) :: asub
      character(len=130) :: aout

!-----------------------------------------------------------------------
! Open file
!-----------------------------------------------------------------------
      unit=20
      aout=trim(asub)//'/coeff'      
      open(unit,file=aout,form='formatted',status='unknown')

!-----------------------------------------------------------------------
! Write coefficient to file
!-----------------------------------------------------------------------
      coe_r=real(traj(ifg)%coe(itraj,istep))
      coe_i=aimag(traj(ifg)%coe(itraj,istep))
      write(unit,'(a)') 'Re, Im'
      write(unit,'(2(F10.7,2x))') coe_r,coe_i

!-----------------------------------------------------------------------
! Close file
!-----------------------------------------------------------------------
      close(unit)

      return

    end subroutine wrcoeff

!#######################################################################

    subroutine copyinpdirs(asub)
      
      use dysonmod

      implicit none

      character(len=120) :: asub
      character(len=150) :: acmnd 

!-----------------------------------------------------------------------
! Columbus, neutral
!-----------------------------------------------------------------------
      acmnd='cp -r '//trim(acol_n)//' '//trim(asub)
      call system(acmnd)

!-----------------------------------------------------------------------
! Columbus, cation
!-----------------------------------------------------------------------
      acmnd='cp -r '//trim(acol_c)//' '//trim(asub)
      call system(acmnd)

!-----------------------------------------------------------------------
! super dyson prep, neutral
!-----------------------------------------------------------------------
      acmnd='cp -r '//trim(aprep_n)//' '//trim(asub)
      call system(acmnd)

!-----------------------------------------------------------------------
! super dyson prep, cation
!-----------------------------------------------------------------------
      acmnd='cp -r '//trim(aprep_c)//' '//trim(asub)
      call system(acmnd)

      return

    end subroutine copyinpdirs

!#######################################################################

    subroutine copygeom(asub,ifg,itraj,istep)
      
      use sysdef
      use trajdef
      use dysonmod

      implicit none

      integer*8                 :: ifg,itraj,istep,unit,i,m,k,spawnstep
      real*8, dimension(natm*3) :: x
      character(len=120)        :: asub
      character(len=130)        :: ageom

!-----------------------------------------------------------------------
! Determine the current Cartesian coordinates
!
! If t<tspawn, then use the coordinates of the parent trajectory
!-----------------------------------------------------------------------
      spawnstep=traj(ifg)%tspawn(itraj)
      if (istep.lt.spawnstep) then
         k=traj(ifg)%ispawn(itraj)
         x=traj(ifg)%r(k,istep,:)
      else
         x=traj(ifg)%r(itraj,istep,:)
      endif

!-----------------------------------------------------------------------
! Columbus, neutral
!-----------------------------------------------------------------------
      unit=20
      ageom=trim(asub)//'/'//trim(acol_n)//'/geom'      
      open(unit,file=ageom,form='formatted',status='unknown')      
      do i=1,natm
         write(unit,'(1x,a2,4x,F3.1,3(F14.8),3x,F11.8)') atlbl(i),atnum(i),&
              (x(m), m=i*3-2,i*3),atmass(i)
      enddo
      close(unit)

!-----------------------------------------------------------------------
! Columbus, cation
!-----------------------------------------------------------------------
      ageom=trim(asub)//'/'//trim(acol_c)//'/geom'      
      open(unit,file=ageom,form='formatted',status='unknown')      
      do i=1,natm
         write(unit,'(1x,a2,4x,F3.1,3(F14.8),3x,F11.8)') atlbl(i),atnum(i),&
              (x(m), m=i*3-2,i*3),atmass(i)
      enddo
      close(unit)

!-----------------------------------------------------------------------
! super dyson prep, neutral
!-----------------------------------------------------------------------
      ageom=trim(asub)//'/'//trim(aprep_n)//'/geom'      
      open(unit,file=ageom,form='formatted',status='unknown')      
      do i=1,natm
         write(unit,'(1x,a2,4x,F3.1,3(F14.8),3x,F11.8)') atlbl(i),atnum(i),&
              (x(m), m=i*3-2,i*3),atmass(i)
      enddo
      close(unit)

!-----------------------------------------------------------------------
! super dyson prep, cation
!-----------------------------------------------------------------------
      ageom=trim(asub)//'/'//trim(aprep_c)//'/geom'      
      open(unit,file=ageom,form='formatted',status='unknown')      
      do i=1,natm
         write(unit,'(1x,a2,4x,F3.1,3(F14.8),3x,F11.8)') atlbl(i),atnum(i),&
              (x(m), m=i*3-2,i*3),atmass(i)
      enddo
      close(unit)

      return

    end subroutine copygeom

!#######################################################################

  end module dysonprep

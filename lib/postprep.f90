!#######################################################################
! postprep: subroutines for the creation of the input for quantum
!           chemistry calculations used in a post-processing step           
!#######################################################################

  module postprep

  contains

!#######################################################################

    subroutine getmaindir(amain,i,k)

      implicit none
      
      integer            :: i,k
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

      integer            :: ifg,itraj,unit
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

      integer            :: ifg,itraj,unit
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

    subroutine getsubdir(asub,i,k,n)

      implicit none

      integer            :: i,k,n
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

    subroutine wrcoeff(asub,ifg,itraj,istep)

      use trajdef

      implicit none

      integer            :: ifg,itraj,istep,unit
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
      use expec, only: ldummy

      implicit none

      integer                   :: ifg,itraj,istep,unit,i,m,k,spawnstep
      real*8, dimension(natm*3) :: x
      real*8, dimension(3)      :: xcom
      character(len=120)        :: asub
      character(len=130)        :: ageom

!-----------------------------------------------------------------------
! Determine the current Cartesian coordinates
!-----------------------------------------------------------------------
      call currgeom(x,ifg,itraj,istep)

!-----------------------------------------------------------------------
! Determing the centre of mass if required
!-----------------------------------------------------------------------
      if (ldummy) call getcom(x,xcom)

!-----------------------------------------------------------------------
! Columbus, neutral
!-----------------------------------------------------------------------
      unit=20
      ageom=trim(asub)//'/'//trim(acol_n)//'/geom'      
      open(unit,file=ageom,form='formatted',status='unknown')
      if (ldummy) write(unit,'(1x,a2,4x,F3.1,3(F14.8),3x,F11.8)') 'X',&
           0.0d0,(xcom(m), m=1,3),0.0d0
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
      if (ldummy) write(unit,'(1x,a2,4x,F3.1,3(F14.8),3x,F11.8)') 'X',&
           0.0d0,(xcom(m), m=1,3),0.0d0
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
      if (ldummy) write(unit,'(1x,a2,4x,F3.1,3(F14.8),3x,F11.8)') 'X',&
           0.0d0,(xcom(m), m=1,3),0.0d0
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
      if (ldummy) write(unit,'(1x,a2,4x,F3.1,3(F14.8),3x,F11.8)') 'X',&
           0.0d0,(xcom(m), m=1,3),0.0d0
      do i=1,natm
         write(unit,'(1x,a2,4x,F3.1,3(F14.8),3x,F11.8)') atlbl(i),atnum(i),&
              (x(m), m=i*3-2,i*3),atmass(i)
      enddo
      close(unit)

      return

    end subroutine copygeom

!#######################################################################
    
    subroutine copygeom_adc(asub,ifg,itraj,istep)
      
      use sysdef
      use trajdef
      use dysonmod
      use expec, only: ldummy

      implicit none

      integer                   :: ifg,itraj,istep,unit,i,m,k,spawnstep
      real*8, dimension(natm*3) :: x
      real*8, dimension(3)      :: xcom
      character(len=120)        :: asub
      character(len=130)        :: ageom

!-----------------------------------------------------------------------
! Determine the current Cartesian coordinates
!-----------------------------------------------------------------------
      call currgeom(x,ifg,itraj,istep)

!-----------------------------------------------------------------------
! Columbus, neutral
!-----------------------------------------------------------------------
      unit=20
      ageom=trim(asub)//'/geom'
      open(unit,file=ageom,form='formatted',status='unknown')

      if (ldummy) then
         call getcom(x,xcom)
         write(unit,'(1x,a2,4x,F3.1,3(F14.8),3x,F11.8)') 'X ',0.0d0,&
              (xcom(m), m=1,3),0.0d0
      endif

      do i=1,natm
         write(unit,'(1x,a2,4x,F3.1,3(F14.8),3x,F11.8)') atlbl(i),atnum(i),&
              (x(m), m=i*3-2,i*3),atmass(i)
      enddo

      close(unit)

      return

    end subroutine copygeom_adc

!#######################################################################

    subroutine copygeom_gamess(asub,ifg,itraj,istep)
      
      use sysdef
      use trajdef
      use dysonmod
      use expec, only: ldummy

      implicit none

      integer                   :: ifg,itraj,istep,unit,i,m,k,spawnstep
      real*8, dimension(natm*3) :: x
      real*8, dimension(3)      :: xcom
      character(len=120)        :: asub
      character(len=130)        :: ageom

!-----------------------------------------------------------------------
! Determine the current Cartesian coordinates
!-----------------------------------------------------------------------
      call currgeom(x,ifg,itraj,istep)

!-----------------------------------------------------------------------
! Columbus, neutral
!-----------------------------------------------------------------------
      unit=20
      ageom=trim(asub)//'/columbus/geom'
      open(unit,file=ageom,form='formatted',status='unknown')

      if (ldummy) then
         call getcom(x,xcom)
         write(unit,'(1x,a2,4x,F3.1,3(F14.8),3x,F11.8)') 'X ',0.0d0,&
              (xcom(m), m=1,3),0.0d0
      endif

      do i=1,natm
         write(unit,'(1x,a2,4x,F3.1,3(F14.8),3x,F11.8)') atlbl(i),atnum(i),&
              (x(m), m=i*3-2,i*3),atmass(i)
      enddo

      close(unit)

      return

    end subroutine copygeom_gamess

!#######################################################################
    
    subroutine getcom(x,xcom)
      
      use sysdef
      
      implicit none
      
      integer                   :: i,j
      real*8, dimension(natm*3) :: x
      real*8, dimension(3)      :: xcom
      real*8                    :: totmass
      
      xcom=0.0d0
      totmass=0.0d0
      do i=1,natm         
         totmass=totmass+atmass(i)
         do j=1,3
            xcom(j)=xcom(j)+x(i*3-3+j)*atmass(i)
         enddo
      enddo
      xcom=xcom/totmass

      return
      
    end subroutine getcom

!#######################################################################

    subroutine wradcinp(asub,ifg,itraj,istep)

      use sysdef
      use trajdef
      use expec
      use iomod
      
      implicit none

      integer                        :: ifg,itraj,istep,unit,i,j,k,&
                                        spawnstep,n,c,ilbl
      real*8, dimension(natm*3)      :: x
      character(len=120)             :: asub
      character(len=130)             :: filename
      character(len=150)             :: dir
      character(len=1), dimension(3) :: dpl
      character(len=3), dimension(4) :: type
      
!-----------------------------------------------------------------------
! Determine the current Cartesian coordinates
!-----------------------------------------------------------------------
      call currgeom(x,ifg,itraj,istep)

!-----------------------------------------------------------------------
! Initialization of various things
!-----------------------------------------------------------------------
      call getfreeunit(unit)
      
      type(1:4)=(/ 'ip ','dav','si ','dys' /)
      dpl(1:3)=(/ 'x','y','z' /)
      
!-----------------------------------------------------------------------
! 1: IP
!-----------------------------------------------------------------------
      if (adcfile(1).ne.'') then
         n=1
      
         ! Make the directory
         dir=trim(asub)//'/adc_'//trim(type(n))
         call makedir(dir)
         
         ! Open the ADC input file
         call getadcfilename(filename,dir,type(n),ifg,itraj,istep,' ')
         open(unit,file=filename,form='formatted',status='unknown')
         
         ! Write the ADC input file
         do i=1,nlines(n)
            if (adcinp(n,i).eq.'$geom') then
               ! Geometry section
               do j=1,natm
                  write(unit,'(a2,3(2x,F14.8))') atlbl(j),&
                       (x(k)*0.529177249d0,k=j*3-2,j*3)
               enddo
            else
               ! Everything else
               write(unit,'(a)') trim(adcinp(n,i))
            endif
         enddo
         
         ! Close the ADC input file
         close(unit)
      endif

!-----------------------------------------------------------------------
! 2 and 3: Davidson and Stieltjes
!-----------------------------------------------------------------------
      do n=2,3

         if (adcfile(n).eq.'') cycle

         ! Loop over dipole components
         do c=1,3

            ! Make the directory
            dir=trim(asub)//'/adc_'//trim(type(n))//'_'//dpl(c)
            call makedir(dir)
            
            ! Open the ADC input file
            call getadcfilename(filename,dir,type(n),ifg,itraj,istep,dpl(c))
            open(unit,file=filename,form='formatted',status='unknown')

            ! Write the ADC input file
            do i=1,nlines(n)
               if (adcinp(n,i).eq.'$geom') then
                  ! Geometry section
                  do j=1,natm
                     write(unit,'(a2,3(2x,F14.8))') atlbl(j),&
                          (x(k)*0.529177249d0,k=j*3-2,j*3)
                  enddo
               else if (index(adcinp(n,i),'dipole_component').ne.0) then
                  ! Dipole component
                  ilbl=index(adcinp(n,i),'=')
                  write(unit,'(a)') trim(adcinp(n,i)(1:ilbl))//dpl(c)
               else
                  ! Everything else
                  write(unit,'(a)') trim(adcinp(n,i))
               endif
            enddo
            
            ! Close the ADC input file
            close(unit)
            
         enddo
            
      enddo

!-----------------------------------------------------------------------
! 4: Dyson
!-----------------------------------------------------------------------
      if (adcfile(4).ne.'') then
         n=4

         ! Make the directory
         dir=trim(asub)//'/adc_'//trim(type(n))
         call makedir(dir)

         ! Open the ADC input file
         call getadcfilename(filename,dir,type(n),ifg,itraj,istep,' ')
         open(unit,file=filename,form='formatted',status='unknown')

         ! Write the ADC input file
         do i=1,nlines(n)
            if (adcinp(n,i).eq.'$geom') then
               ! Geometry section
               do j=1,natm
                  write(unit,'(a2,3(2x,F14.8))') atlbl(j),&
                       (x(k)*0.529177249d0,k=j*3-2,j*3)
               enddo
            else
               ! Everything else
               write(unit,'(a)') trim(adcinp(n,i))
            endif
         enddo
         
         ! Close the ADC input file
         close(unit)
      endif

      return
      
    end subroutine wradcinp

!#######################################################################
    
    subroutine currgeom(x,ifg,itraj,istep)

      use sysdef
      use trajdef
      
      implicit none
      
      integer                   :: ifg,itraj,istep,spawnstep,k,k1,&
                                   spawncurr,trajlbl
      real*8, dimension(natm*3) :: x
      
!-----------------------------------------------------------------------
! Take the centre of the closest populated ancestor of the current 
! trajectory
!-----------------------------------------------------------------------
      spawncurr=traj(ifg)%tspawn(itraj)
      if (istep.ge.spawncurr) then
         ! The current trajectory has spawned: use its centre
         trajlbl=itraj
      else
         ! The trajectory has not yet spawned: search for the closest
         ! spawned ancestor and use its centre
         k=itraj
         ! Check whether the next ancestor down the line has spawned yet,
         ! and if so, then select it
10       continue
         ! Index of the next parent up the line
         k1=traj(ifg)%ispawn(k)
         ! Timestep at which the next parent was spawned
         spawnstep=traj(ifg)%tspawn(k1)
         ! Accept if the next parent has spawned, else reject
         ! and move onto the next parent up the line
         if (istep.ge.spawnstep) then
            trajlbl=k1
         else
            k=k1
            goto 10
         endif
      endif

      ! Assign the Cartesian coordinates
      x=traj(ifg)%r(trajlbl,istep,:)
      
!      if (sum(x).eq.0.0d0) then
!         print*,itraj,istep
!         print*,"HERE"
!         STOP
!      endif

      return
           
    end subroutine currgeom

!#######################################################################

    subroutine getadcfilename(filename,dir,type,ifg,itraj,istep,dpl)

      implicit none

      integer            :: ifg,itraj,istep,k
      character(len=130) :: filename
      character(len=150) :: dir
      character(len=3)   :: type
      character(len=1)   :: dpl

      filename=trim(dir)//'/adc_'

      k=len_trim(filename)
      if (ifg.lt.10) then
         write(filename(k+1:k+5),'(a4,i1)') 'ifg0',ifg
      else
         write(filename(k+1:k+5),'(a3,i2)') 'ifg',ifg
      endif

      k=len_trim(filename)
      if (itraj.lt.10) then
         write(filename(k+1:k+7),'(a6,i1)') '_traj0',itraj
      else
         write(filename(k+1:k+7),'(a5,i2)') '_traj',itraj
      endif

      k=len_trim(filename)
      if (istep.lt.10) then
         write(filename(k+1:k+10),'(a9,i1)') '_step0000',istep
      else if (istep.lt.100) then
         write(filename(k+1:k+10),'(a8,i2)') '_step000',istep
      else if (istep.lt.1000) then
         write(filename(k+1:k+10),'(a7,i3)') '_step00',istep
      else if (istep.lt.10000) then
         write(filename(k+1:k+10),'(a6,i4)') '_step0',istep
      else
         write(filename(k+1:k+10),'(a5,i6)') '_step',istep
      endif

      if (type.eq.'ip'.or.type.eq.'dys') then
         filename=trim(filename)//'_'//trim(type)//'.inp'
      else
         filename=trim(filename)//'_'//trim(type)//'_'//trim(dpl)//'.inp'
      endif

      return

    end subroutine getadcfilename

!#######################################################################

    subroutine wrgamessinp(asub,ifg,itraj,istep)

      use sysdef
      use trajdef
      use expec
      use iomod
      
      implicit none

      integer                   :: ifg,itraj,istep,unit,n,i,j,k
      real*8, dimension(natm*3) :: x
      character(len=120)        :: asub
      character(len=170)        :: filename

!-----------------------------------------------------------------------
! Determine the current Cartesian coordinates
!-----------------------------------------------------------------------
      call currgeom(x,ifg,itraj,istep)

!-----------------------------------------------------------------------
! Write the GAMESS input files
!-----------------------------------------------------------------------
      call getfreeunit(unit)

      ! Loop over GAMESS input files
      do n=1,ngamessinp

         ! Open the current file
         filename=trim(asub)//'/gamess/'//trim(gamessfile(n))
         open(unit,file=filename,form='formatted',status='unknown')
         
         ! Write to the current file
         do i=1,ngmslines(n)
            if (index(gmsinp(n,i),'$geom').ne.0) then               
               ! Geometry section
               do j=1,natm
                  write(unit,'(1x,a2,F4.1,3(2x,F14.8))') &
                       atlbl(j),atnum(j),&
                       (x(k)*0.529177249d0,k=j*3-2,j*3)
               enddo
            else
               ! Everything else
               write(unit,'(a)') trim(gmsinp(n,i))
            endif
         enddo

         ! Close the current file
         close(unit)

      enddo

      return

    end subroutine wrgamessinp

!#######################################################################

  end module postprep

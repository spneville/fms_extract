!######################################################################
! Calculation of the time-dependent incoherent average of orbital
! energies using the output of the ADC code
!######################################################################

  module adc_orben

    use expec
    use sysdef
    
    implicit none

    save

    integer                                       :: nmaindir,nifg,nmo
    integer, dimension(:), allocatable            :: staindx,ifgindx,&
                                                     step,icontrib,&
                                                     nifg_curr
    real*8, dimension(:,:), allocatable           :: cnorm,avorben
    real*8, dimension(:,:), allocatable           :: emo
    complex*16, dimension(:), allocatable         :: coeff
    character(len=120), dimension(:), allocatable :: amaindir

  contains

!######################################################################

    subroutine calc_adc_orben

      implicit none

      integer                                       :: i,j,k,ifg,itmp,&
                                                       nsubdir&
                                                       &,ncontrib
      real*8                                        :: csq,chk
      character(len=120), dimension(:), allocatable :: asubdir

!-----------------------------------------------------------------------
! Read the names of the main directories
!-----------------------------------------------------------------------
      call rdmaindirfile

!-----------------------------------------------------------------------
! Determine the number of IFGs being considered, which may be less than
! nintraj
!-----------------------------------------------------------------------
      call getnifg

!-----------------------------------------------------------------------
! Determine the no. MOs
!-----------------------------------------------------------------------
      call getnmo

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      call alloc_orben

!-----------------------------------------------------------------------
! Determine the norm of the C-vectors for each IFG at each timestep
! using only the trajectories for which ADC results are available
!-----------------------------------------------------------------------
      call getcnorm

!-----------------------------------------------------------------------
! Determine the state indices for each trajectory (main directory) being
! considered
!-----------------------------------------------------------------------
      call getstaindx

!-----------------------------------------------------------------------
! Determine the number of IFGs that contribute at each timestep
!-----------------------------------------------------------------------
      call getnifg_curr

!-----------------------------------------------------------------------
! Loop over the main trajectories (main directories), reading the
! orbital energies for each and 
!-----------------------------------------------------------------------
      do i=1,nmaindir
         
         ! Ouput our progress
         write(6,'(2a)') 'Processing directory: ',trim(amaindir(i))

         ! Get the list of timesteps/subdirectories
         call getsubdirs(amaindir(i),nsubdir,asubdir,step)

         ! Read the coefficients for each timestep/subdirectory
         call getcoeff(coeff,amaindir(i),asubdir,nsubdir)

         ! Determine which timesteps/subdirectories contribute
         call getcontrib_orben(ncontrib,amaindir(i),asubdir,&
              nsubdir)
         
         ! Cycle if no timesteps/subdirectories contribute to the
         ! spectrum
         if (ncontrib.eq.0) cycle

         ! Read the orbital energies for each timestep/subdirectory
         call rdorben(i,asubdir,nsubdir)

         ! Add the contributions of the current orbital energies
         ! to the average orbital energies
         do j=1,nsubdir

            if (icontrib(j).eq.0) cycle

            csq=real(coeff(j)*conjg(coeff(j)))
            ifg=ifgindx(i)

            avorben(:,j)=avorben(:,j)&
                 +emo(:,j)*(csq/cnorm(ifg,j))/nifg_curr(j)

         enddo

      enddo

!-----------------------------------------------------------------------
! Output the average orbital energies
!-----------------------------------------------------------------------
      call wravorben

      return

    end subroutine calc_adc_orben

!######################################################################

    subroutine alloc_orben

      implicit none

      integer :: itmp

      ! State indices
      allocate(staindx(nmaindir))
      staindx=0

      ! Average orbital energies
      allocate(avorben(nmo,1+nstep/dstep))
      avorben=0.0d0

      ! No. IFGs that contribute at each timestep
      allocate(nifg_curr(1+nstep/dstep))
      nifg_curr=0

      return

    end subroutine alloc_orben

!######################################################################

    subroutine rdmaindirfile

      use expec
      use parsemod

      implicit none

      integer                                       :: unit,n
      character(len=120)                            :: string

      integer                                       :: inkw
      integer, parameter                            :: maxkw=200
      integer, dimension(maxkw)                     :: ilkw
      character(len=120), dimension(maxkw)          :: keyword

!-----------------------------------------------------------------------
! Open file
!-----------------------------------------------------------------------
      unit=20
      open(unit,file=adcdir_file,form='formatted',status='old')

!-----------------------------------------------------------------------
! Determine the number of main directories and allocate the amaindir
! array
!-----------------------------------------------------------------------
      nmaindir=0
5     continue
      call rdinp(unit,keyword,inkw,ilkw)

      if (keyword(1).ne.'end-file') then
         nmaindir=nmaindir+1
         goto 5
      endif

      if (allocated(amaindir)) deallocate(amaindir)
      allocate(amaindir(nmaindir))

!-----------------------------------------------------------------------
! Read the main directory names
!-----------------------------------------------------------------------
      rewind(unit)
      n=0
10    continue
      call rdinp(unit,keyword,inkw,ilkw)
      if (keyword(1).ne.'end-file') then
         n=n+1
         amaindir(n)=keyword(1)         
         goto 10
      endif

!-----------------------------------------------------------------------
! Close file
!-----------------------------------------------------------------------
      close(unit)

      return

    end subroutine rdmaindirfile

!######################################################################

    subroutine getnifg

      use sysdef

      implicit none

      integer                     :: i,k,unit
      integer, dimension(nintraj) :: cnt
      character(len=130)          :: ain

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      if (allocated(ifgindx)) deallocate(ifgindx)
      allocate(ifgindx(nmaindir))
             
!-----------------------------------------------------------------------
! Loop over the main directories, reading the ifg_number file for each
!-----------------------------------------------------------------------
      unit=20
      cnt=0
      do i=1,nmaindir
         ain=trim(amaindir(i))//'/ifg_number'
         open(unit,file=ain,form='formatted',status='old')
         read(unit,*) k
         cnt(k)=cnt(k)+1
         ifgindx(i)=k
         close(unit)
      enddo

      nifg=0
      do i=1,nintraj
         if (cnt(i).gt.0) nifg=nifg+1
      enddo
      
      return

    end subroutine getnifg

!######################################################################

    subroutine getcnorm

      use sysdef 
      use expec
      use iomod
      
      implicit none

      integer                                       :: itmp,n,k,ifg,&
                                                       nsubdir,ncontrib
      integer, dimension(:), allocatable            :: step
      complex*16, dimension(:), allocatable         :: coeff
      character(len=120), dimension(:), allocatable :: asubdir

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      itmp=nstep/dstep+1
      if (allocated(cnorm)) deallocate(cnorm)
      allocate(cnorm(nintraj,itmp))
      cnorm=0.0d0

!-----------------------------------------------------------------------
! Determine the approximate norm for each IFG at each timestep
!-----------------------------------------------------------------------
      ! Loop over main directories/trajectories
      do n=1,nmaindir

         ! IFG no. for the current trajectory
         ifg=ifgindx(n)
         
         ! Get the list of timesteps/subdirectories for the
         ! current trajectory
         call getsubdirs(amaindir(n),nsubdir,asubdir,step)

         ! Read the coefficients for each timestep/subdirectory
         call getcoeff(coeff,amaindir(n),asubdir,nsubdir)

         ! Determine which timesteps/subdirectories contribute
         call getcontrib_orben(ncontrib,amaindir(n),asubdir,&
              nsubdir)

         ! Add the current contribution to the cnorm array
         do k=1,nsubdir
            if (icontrib(k).eq.0) cycle
            cnorm(ifg,k)=cnorm(ifg,k)+real(conjg(coeff(k))*coeff(k))
         enddo

      enddo

      cnorm=sqrt(cnorm)

      return

    end subroutine getcnorm

!######################################################################

    subroutine getnifg_curr

      implicit none

      integer                                       :: n,k,i,ifg,&
                                                       nsubdir,&
                                                       ncontrib
      integer, dimension(:), allocatable            :: step
      integer, dimension(:,:), allocatable          :: ifgcontrib
      character(len=120), dimension(:), allocatable :: asubdir

!----------------------------------------------------------------------
! Allocation and initialisation
!----------------------------------------------------------------------
      nifg_curr=0
      allocate(ifgcontrib(nintraj,1+nstep/dstep))
      ifgcontrib=0

!----------------------------------------------------------------------
! Determine how many IFGs contribute at each timestep (nifg_curr)
!----------------------------------------------------------------------
      ! Loop over main directories/trajectories
      do n=1,nmaindir
         
         ! IFG no. for the current trajectory
         ifg=ifgindx(n)

         ! Get the list of timesteps/subdirectories for the
         ! current trajectory
         call getsubdirs(amaindir(n),nsubdir,asubdir,step)

         ! Determine which timesteps/subdirectories contribute
         call getcontrib_orben(ncontrib,amaindir(n),asubdir,&
              nsubdir)

         ! Add the contribution of the current trajectory to nifg_curr
         do k=1,nsubdir
            ifgcontrib(ifg,k)=ifgcontrib(ifg,k)+icontrib(k)
         enddo

      enddo

      ! Determine nifg_curr
      do k=1,1+nstep/dstep         
         do i=1,nintraj
            if (ifgcontrib(i,k).gt.0) nifg_curr(k)=nifg_curr(k)+1
         enddo
      enddo

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
      deallocate(ifgcontrib)

      return

    end subroutine getnifg_curr

!######################################################################

    subroutine getsubdirs(amaindir,nsubdir,asubdir,step)

      implicit none

      integer                                       :: unit,nsubdir,i
      integer, dimension(:), allocatable            :: step
      character(len=120)                            :: amaindir
      character(len=130)                            :: alist,string
      character(len=120), dimension(:), allocatable :: asubdir
      logical(kind=4)                               :: lexists
      
!-----------------------------------------------------------------------
! Read the sublist file
!-----------------------------------------------------------------------
      alist=trim(amaindir)//'/sublist'

      inquire(file=alist,exist=lexists)
      if (.not.lexists) then
         write(6,'(2a)') 'sublist file not found in the directory ',&
              trim(amaindir)
         STOP
      endif

      unit=20
      open(unit,file=alist,form='formatted',status='old')

      nsubdir=0
10    continue
      read(unit,'(a)',end=100) string
      nsubdir=nsubdir+1
      goto 10

100   continue
      if (allocated(asubdir)) deallocate(asubdir)
      allocate(asubdir(nsubdir))
      
      rewind(unit)
      do i=1,nsubdir
         read(unit,'(a)') asubdir(i)
      enddo

!-----------------------------------------------------------------------
! Determine the timesteps
!-----------------------------------------------------------------------
      if (allocated(step)) deallocate(step)
      allocate(step(nsubdir))

      do i=1,nsubdir
         read(asubdir(i)(5:9),*) step(i)
      enddo

      close(unit)

      return

    end subroutine getsubdirs

!######################################################################

    subroutine getcoeff(coeff,amaindir,asubdir,nsubdir)

      implicit none

      integer                                :: nsubdir,i,unit
      real*8                                 :: cr,ci
      complex*16, dimension(:), allocatable  :: coeff
      character(len=120)                     :: amaindir
      character(len=120), dimension(nsubdir) :: asubdir
      character(len=270)                     :: acoeff

!-----------------------------------------------------------------------
! Allocate the coeff array
!-----------------------------------------------------------------------
      if (allocated(coeff)) deallocate(coeff)
      allocate(coeff(nsubdir))
      coeff=0.0d0

!-----------------------------------------------------------------------
! Loop over timesteps/subdirectories and read the coefficents for each
!-----------------------------------------------------------------------
      unit=20
      do i=1,nsubdir
         acoeff=trim(amaindir)//'/'//trim(asubdir(i))//'coeff'
         open(unit,file=acoeff,form='formatted',status='old')
         read(unit,*)
         read(unit,*) cr,ci
         coeff(i)=dcmplx(cr,ci)
         close(unit)
      enddo

      return

    end subroutine getcoeff

!######################################################################

    subroutine getcontrib_orben(ncontrib,amaindir,asubdir,nsubdir)      

      implicit none

      integer                                :: ncontrib,nsubdir,i
      character(len=120)                     :: amaindir
      character(len=120), dimension(nsubdir) :: asubdir
      character(len=250)                     :: filename
      logical                                :: exists
      
!-----------------------------------------------------------------------
! Allocate the icontrib array
!-----------------------------------------------------------------------
      if (allocated(icontrib)) deallocate(icontrib)
      allocate(icontrib(nsubdir))
      icontrib=1

!-----------------------------------------------------------------------
! Loop over timesteps/subdirectories and determine whether each one
! contributes, i.e., whether the ADC calculation took place
!-----------------------------------------------------------------------
      ncontrib=0
      do i=1,nsubdir
      
         filename=trim(amaindir)//'/'//trim(asubdir(i))//'/adc_dav_x/osc.dat'
         inquire(file=filename,exist=exists)
         if (.not.exists) goto 10
      
         ncontrib=ncontrib+1
         cycle

10       continue
         icontrib(i)=0

      enddo

      return

    end subroutine getcontrib_orben

!######################################################################

    subroutine getstaindx

      implicit none

      integer            :: i,unit
      character(len=130) :: ain

      unit=20

      do i=1,nmaindir
         ain=trim(amaindir(i))//'/state_id'
         open(unit,file=ain,form='formatted',status='old')
         read(unit,*) staindx(i)
         close(unit)
      enddo

      return

    end subroutine getstaindx

!######################################################################

    subroutine getnmo
    
      implicit none

      integer                                       :: nsubdir,k1,k2,&
                                                       k3,unit
      character(len=120), dimension(:), allocatable :: asubdir
      character(len=250)                            :: filename,string

!-----------------------------------------------------------------------
! Read the no. MOs from the ADC output for the first trajectory for
! the first timestep
!-----------------------------------------------------------------------
      ! Get the list of subdirectories for the 1st trajectory
      call getsubdirs(amaindir(1),nsubdir,asubdir,step)

      ! Open one of the ADC log files (any will do)
      k1=index(amaindir(1),'/ifg')+1
      k2=len_trim(amaindir(1))
      if (index(asubdir(1),'/').eq.0) then
         k3=len_trim(asubdir(1))
      else
         k3=len_trim(asubdir(1))-1
      endif
      
      filename=trim(amaindir(1))//'/'//trim(asubdir(1)) &
           //'/adc_dav_x/adc_'//amaindir(1)(k1:k2)//'_' &
           //asubdir(1)(1:k3)//'_dav_x.log'
      
      unit=77
      open(unit,file=filename,form='formatted',status='old')

      ! Read the no. MOs from ADC log file
10    read(unit,'(a)') string
      if (index(string,'vectors in the set').eq.0) goto 10
      read(string,'(25x,i4)') nmo

      ! Close the ADC log file
      close(unit)

      return

    end subroutine getnmo

!######################################################################

    subroutine rdorben(itraj,asubdir,nsubdir)

      implicit none

      integer                                :: itraj,nsubdir,i
      character(len=120), dimension(nsubdir) :: asubdir

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
      if (allocated(emo)) deallocate(emo)
      allocate(emo(nmo,nsubdir))
      emo=0.0d0

!----------------------------------------------------------------------
! For the current trajectory/main directory, read the orbital energies
! for each timestep/subdirectory
!----------------------------------------------------------------------
      do i=1,nsubdir
         
         ! Cycle if the current timestep doesn't contribute
         if (icontrib(i).eq.0) cycle

         ! Read the orbital energies for the current timestep
         call rdemo(itraj,i,asubdir(i))        

      enddo

      return

    end subroutine rdorben

!######################################################################

    subroutine rdemo(itraj,istep,asubdir)

      implicit none

      integer                :: itraj,istep,k1,k2,k3,unit,i,n
      character(len=120)     :: asubdir,string
      character(len=250)     :: filename


!-----------------------------------------------------------------------
! Open one of the ADC log files (any will do)
!-----------------------------------------------------------------------
      k1=index(amaindir(itraj),'/ifg')+1
      k2=len_trim(amaindir(itraj))
      if (index(asubdir,'/').eq.0) then
         k3=len_trim(asubdir)
      else
         k3=len_trim(asubdir)-1
      endif
      
      filename=trim(amaindir(itraj))//'/'//trim(asubdir) &
           //'/adc_dav_x/adc_'//amaindir(itraj)(k1:k2)//'_' &
           //asubdir(1:k3)//'_dav_x.log'

      unit=77
      open(unit,file=filename,form='formatted',status='old')

!-----------------------------------------------------------------------
! Read the MO energies from the ADC log file
!-----------------------------------------------------------------------
      ! Read to the MO section
10    read(unit,'(a)') string
      if (index(string,'orbital energy').eq.0) goto 10
      read(unit,*)

      ! Read the MO energies
      n=0
15    continue
      do i=1,2
         read(unit,'(a)') string
      enddo
      if (index(string,'***********************************').eq.0) then
         n=n+1
         read(string,'(43x,F15.10)') emo(n,istep)
         goto 15
      endif

!-----------------------------------------------------------------------
! Close the ADC log file
!-----------------------------------------------------------------------
      close(unit)

      return

    end subroutine rdemo

!######################################################################
    
    subroutine wravorben

      implicit none

      integer :: i,j,unit

!----------------------------------------------------------------------
! Open the average orbital energy file
!----------------------------------------------------------------------
      unit=77
      open(unit,file='avorben.dat',form='formatted',status='unknown')

!----------------------------------------------------------------------
! Write the average orbital energy file
!----------------------------------------------------------------------
      do i=1,(nstep/dstep)+1
         !write(unit,'(11(x,F10.7))') (i-1)*dt*dstep/41.341375d0,&
         !     (avorben(j,i),j=1,10)
         write(unit,'(3(4x,F12.7))'),(i-1)*dt*dstep/41.341375d0,&
              (avorben(2,i)-avorben(1,i))*27.2114d0
      enddo
      
!----------------------------------------------------------------------
! Close the average obital energy file
!----------------------------------------------------------------------
      close(unit)

      return

    end subroutine wravorben

!######################################################################

  end module adc_orben

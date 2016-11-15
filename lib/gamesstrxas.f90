  module gamesstrxas

    implicit none

    save

    integer                                       :: nmaindir,nifg,&
                                                     ne,nt,maxfunc,nfunc
    integer, dimension(:), allocatable            :: staindx,ifgindx
    integer                                       :: failunit,okunit,cifunit
    real*8, dimension(:,:), allocatable           :: spec,par,cnorm
    real*8, parameter                             :: eh2ev=27.2113845d0
    real*8, parameter                             :: c_au=137.03604d0
    real*8                                        :: pi=3.14159265358979d0
    character(len=120), dimension(:), allocatable :: amaindir

  contains

!#######################################################################
! gamess_trxas: Calculates the pre-edge time resolved X-ray absorption
!               spectrum using Columbus-MRCI/GAMESS-ORMAS
!               cross-sections
!#######################################################################

    subroutine gamess_trxas

      implicit none

!-----------------------------------------------------------------------
! Read the names of the main directories
!-----------------------------------------------------------------------
      call rdmaindirfile

!-----------------------------------------------------------------------
! Calculate the spectrum
!-----------------------------------------------------------------------
      call calc_spectrum

      return
      
    end subroutine gamess_trxas

!#######################################################################

    subroutine rdmaindirfile
       
      use expec
      use parsemod

      implicit none
       
      integer                              :: unit,n
      character(len=120)                   :: string

      integer                              :: inkw
      integer, parameter                   :: maxkw=200
      integer, dimension(maxkw)            :: ilkw
      character(len=120), dimension(maxkw) :: keyword
      
!-----------------------------------------------------------------------
! Open file
!-----------------------------------------------------------------------
      unit=20
      open(unit,file=gmsdir_file,form='formatted',status='old')
       
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
     
!#######################################################################

    subroutine calc_spectrum

      use expec
      use trajdef
      use sysdef
      use projmod

      implicit none

      integer                                       :: i,nsubdir,&
                                                       ncontrib
      integer, dimension(:), allocatable            :: step,icontrib
      real*8, dimension(:), allocatable             :: ip,einit
      complex*16, dimension(:), allocatable         :: coeff
      character(len=120), dimension(:), allocatable :: asubdir

!-----------------------------------------------------------------------
! Convert the energy bounds to a.u.
!
! Note that, quite confusingly, we work with energies in a.u. and time
! in fs.
!-----------------------------------------------------------------------
      egrid(1:2)=egrid(1:2)/eh2ev

!-----------------------------------------------------------------------
! Determine the number of IFGs being considered, which may be less than
! nintraj
!-----------------------------------------------------------------------
      call getnifg

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      ! Electronic state indices
      allocate(staindx(nmaindir))

      ! TRXAS arrays
      ne=int(egrid(3))
      nt=int(tgrid(3))
      allocate(spec(ne,nt))
      spec=0.0d0

!-----------------------------------------------------------------------
! Determine the state indices for each trajectory (main directory) being
! considered
!-----------------------------------------------------------------------
      call getstaindx

!-----------------------------------------------------------------------
! Loop over the main trajectories (main directories), reading the
! TDMs, energy differences and coefficients for each timestep
!-----------------------------------------------------------------------
      nfunc=0
      do i=1,nmaindir

         write(6,'(2a)') 'Processing directory: ',trim(amaindir(i))
         
         ! Get the list of timesteps/subdirectories
         call getsubdirs(amaindir(i),nsubdir,asubdir,step)

         ! Read the coefficients for each timestep/subdirectory
         call getcoeff(coeff,amaindir(i),asubdir,nsubdir)

         ! Determine which timesteps/subdirectories contribute to the
         ! spectrum
         call getcontrib_gmsspec(ncontrib,icontrib,coeff,nsubdir)

         ! Cycle if no timesteps/subdirectories contribute to the
         ! spectrum
         if (ncontrib.eq.0) cycle

         ! Read the ionisation potential for each timestep/subdirectory
         call getip_gms(ip,amaindir(i),asubdir,nsubdir,icontrib)


      enddo

      STOP


      return
       
    end subroutine calc_spectrum

!#######################################################################

        subroutine getnifg

      use sysdef

      implicit none

      integer                     :: i,k,unit
      integer, dimension(nintraj) :: cnt
      character(len=130)          :: ain

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
     
!#######################################################################

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

!#######################################################################

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

!#######################################################################

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

!#######################################################################

    subroutine getcontrib_gmsspec(ncontrib,icontrib,coeff,nsubdir)

      implicit none

      integer                            :: ncontrib,nsubdir,i
      integer, dimension(:), allocatable :: icontrib
      real*8, parameter                  :: thrsh=0.001d0
      complex*16, dimension(nsubdir)     :: coeff

!-----------------------------------------------------------------------
! Allocate the icontrib array
!-----------------------------------------------------------------------
      if (allocated(icontrib)) deallocate(icontrib)
      allocate(icontrib(nsubdir))
      icontrib=0

!-----------------------------------------------------------------------
! Determine which timesteps/subdirectories contribute to the spectrum
!-----------------------------------------------------------------------
      ncontrib=0
      do i=1,nsubdir
         if (abs(coeff(i)).gt.thrsh) then
            ncontrib=ncontrib+1
            icontrib(i)=1
         endif
      enddo

      return

    end subroutine getcontrib_gmsspec

!#######################################################################

    subroutine getip_gms(ip,amaindir,asubdir,nsubdir,icontrib)

      implicit none

      integer                                :: nsubdir,unit,i
      integer, dimension(nsubdir)            :: icontrib
      real*8, dimension(:), allocatable      :: ip
      real*8                                 :: e0
      complex*16, dimension(nsubdir)         :: coeff
      character(len=120)                     :: amaindir,string
      character(len=120), dimension(nsubdir) :: asubdir
      character(len=250)                     :: filename
      
!-----------------------------------------------------------------------
! Allocate the IP and E0 arrays
!-----------------------------------------------------------------------
      if (allocated(ip)) deallocate(ip)
      allocate(ip(nsubdir))
      ip=0.0d0

!-----------------------------------------------------------------------
! Loop over timesteps/subdirectories and read the IP for each
!-----------------------------------------------------------------------
      unit=20

      do i=1,nsubdir

         if (icontrib(i).eq.0) cycle

         filename=trim(amaindir)//'/'//trim(asubdir(i))//&
              '/gamess/core_ion_exci.dat'

         open(unit,file=filename,form='formatted',status='old')

         read(unit,'(23x,F16.10)') ip(i)
         print*,trim(filename),ip(i)

         close(unit)

      enddo

      return

    end subroutine getip_gms

!#######################################################################

  end module gamesstrxas

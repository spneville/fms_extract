  module gamesstrxas

    implicit none

    save

    integer                                       :: nmaindir,nifg,&
                                                     ne,nt,maxfunc,nfunc
    integer, dimension(:), allocatable            :: staindx,ifgindx
    integer                                       :: failunit,okunit,cifunit
    integer, parameter                            :: maxsta=100
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
                                                       ncontrib,itmp
      integer, dimension(:), allocatable            :: step,icontrib,&
                                                       ncore
      real*8, dimension(:), allocatable             :: ip,einit
      real*8, dimension(:,:), allocatable           :: ecore,tdmsq
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
! Determine the norm of the C-vectors for each IFG at each timestep
! using only the trajectories for which GAMESS/ORMAS results are
! available
!-----------------------------------------------------------------------
      call getcnorm_gms

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

         ! Read the core-excitation energies for each
         ! timesetep/subdirectory
         ! N.B., this needs to be done before getcontrib_gmsspec
         ! is called so that we know whether the ORMAS calculation
         ! failed
         call getcoreexci_gms(amaindir(i),asubdir,nsubdir,ecore,&
              ncore,icontrib,coeff)

         ! Determine which timesteps/subdirectories contribute to the
         ! spectrum
         call getcontrib_gmsspec(ncontrib,icontrib,coeff,nsubdir,&
              ecore,ncore,amaindir(i),asubdir)

         ! Cycle if no timesteps/subdirectories contribute to the
         ! spectrum
         if (ncontrib.eq.0) cycle

         ! Read the ionisation potential for each timestep/subdirectory
         call getip_gms(ip,amaindir(i),asubdir,nsubdir,icontrib)

         ! Read the energy of the initial state
         call geteinit(amaindir(i),asubdir,nsubdir,icontrib,&
              einit,ncore,staindx(i))

         ! Read the transition dipole moments for each
         ! timesetep/subdirectory
         call gettdmsq(amaindir(i),asubdir,nsubdir,icontrib,&
              ncore,tdmsq)

         ! Calculate the contribution of the current trajectory to
         ! the TRXAS
         call trxas_currifg(nsubdir,icontrib,ncore,ecore,tdmsq,&
              coeff,step,ip,ifgindx(i),einit)

      enddo

      ! Write the TRXAS to file
      call wrtrxas

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

    subroutine getcnorm_gms

      use sysdef
      use expec
      use iomod

      implicit none
      
      integer                                       :: n,k,ifg,itmp,&
                                                       nsubdir,&
                                                       ncontrib
      integer, dimension(:), allocatable            :: step,icontrib,&
                                                       ncore
      real*8                                        :: ftmp
      real*8, dimension(:,:), allocatable           :: ecore
      complex*16, dimension(:), allocatable         :: coeff
      character(len=120), dimension(:), allocatable :: asubdir

! THIS IS QUITE DODGY AS WE ASSUME NON-OVERLAPING BASIS FUNCTIONS

!----------------------------------------------------------------------
! Allocate the cnorm array
!----------------------------------------------------------------------
      itmp=nstep/dstep+1
      allocate(cnorm(nifg,itmp))
      cnorm=0.0d0

!----------------------------------------------------------------------
! Fill in the cnorm array
!----------------------------------------------------------------------
      ! Loop over main directories/trajectories
      do n=1,nmaindir

         ! IFG no. for the current trajectory
         ifg=ifgindx(n)

         ! Get the list of timesteps/subdirectories for the
         ! current trajectory
         call getsubdirs(amaindir(n),nsubdir,asubdir,step)
         
         ! Read the coefficients for each timestep/subdirectory
         call getcoeff(coeff,amaindir(n),asubdir,nsubdir)

         ! Read the core-excitation energies for each
         ! timesetep/subdirectory
         ! N.B., this needs to be done before getcontrib_gmsspec
         ! is called so that we know whether the ORMAS calculation
         ! failed
         call getcoreexci_gms(amaindir(n),asubdir,nsubdir,ecore,&
              ncore,icontrib,coeff)

         ! Determine which timesteps/subdirectories contribute to the
         ! spectrum
         call getcontrib_gmsspec(ncontrib,icontrib,coeff,nsubdir,&
              ecore,ncore,amaindir(n),asubdir)

         ! Add the current contribution to the cnorm array
         do k=1,nsubdir
            if (icontrib(k).eq.0) cycle
            itmp=step(k)/dstep+1
            cnorm(ifg,itmp)=cnorm(ifg,itmp)+conjg(coeff(k))*coeff(k)
         enddo
         
      enddo

      cnorm=sqrt(cnorm)

      return

    end subroutine getcnorm_gms
     
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

    subroutine getcontrib_gmsspec(ncontrib,icontrib,coeff,nsubdir,&
         ecore,ncore,amaindir,asubdir)

      implicit none

      integer                                :: ncontrib,nsubdir,i,j,&
                                                unit
      integer, dimension(:), allocatable     :: icontrib
      integer, dimension(nsubdir)            :: ncore
      real*8, dimension(maxsta,nsubdir)      :: ecore
      real*8, parameter                      :: thrsh=0.01d0
      real*8                                 :: diff
      complex*16, dimension(nsubdir)         :: coeff
      character(len=120)                     :: amaindir
      character(len=120), dimension(nsubdir) :: asubdir
      character(len=250)                     :: filename,string
      logical                                :: found,empty

!-----------------------------------------------------------------------
! Allocate the icontrib array
!-----------------------------------------------------------------------
      if (allocated(icontrib)) deallocate(icontrib)
      allocate(icontrib(nsubdir))
      icontrib=0

!-----------------------------------------------------------------------
! Determine which timesteps/subdirectories contribute to the spectrum
!-----------------------------------------------------------------------
      ! Criterion 1: magnitude of the expansion coefficient
      ncontrib=0
      do i=1,nsubdir
         if (abs(coeff(i)).gt.thrsh) then
            ncontrib=ncontrib+1
            icontrib(i)=1
         endif
      enddo

      ! Criterion 2: success of the core-excited ORMAS calculation.
      ! The ORMAS eigensolver has a tendency to diverge at times,
      ! and this is what we check for here.
      ! Note that here we assume that the ORMAS calculation at the
      ! 1st timestep always succeeds, which isn't really ideal...
      do i=1,nsubdir
         do j=1,ncore(i)
            diff=(ecore(j,i)-ecore(1,1))*eh2ev
            if (diff.gt.100.0d0) then
               icontrib(i)=0
               ncontrib=ncontrib-1
               exit
            endif
         enddo
      enddo

      ! Criterion 3: success of the core-ionised ORMAS calculation.
      ! This seems less suscepible to going wrong than the core-
      ! excited calculation, so we will simply check to see if the
      ! core_ion_exci.dat file exists
      unit=20
      do i=1,nsubdir
         ! (i) Does the file exist?
         filename=trim(amaindir)//'/'//trim(asubdir(i))//&
              '/gamess/core_ion_exci.dat'
         inquire(file=filename,exist=found)
         if (.not.found.and.icontrib(i).eq.1) then
            icontrib(i)=0
            ncontrib=ncontrib-1
         endif
         ! (ii) Is the file non-empty
         if (found) then
            empty=.true.
            open(unit,file=filename,form='formatted',status='old')
10          read(unit,'(a)',end=999) string
            if (index(string,'ENERGY=').ne.0) then
               empty=.false.
            else
               goto 10
            endif
999         continue
            if (empty.and.icontrib(i).eq.1) then
               icontrib(i)=0
               ncontrib=ncontrib-1
            endif            
            close(unit)
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
! Allocate the IP array
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

         close(unit)

      enddo

      return

    end subroutine getip_gms

!#######################################################################

    subroutine getcoreexci_gms(amaindir,asubdir,nsubdir,ecore,ncore,&
         icontrib,coeff)

      implicit none

      integer                                :: nsubdir,i,n,&
                                               unit,ierr
      integer, dimension(:), allocatable     :: ncore
      integer, dimension(nsubdir)            :: icontrib
      real*8, dimension(:,:), allocatable    :: ecore
      real*8, parameter                      :: thrsh=0.001d0
      complex*16, dimension(nsubdir)         :: coeff
      character(len=120)                     :: amaindir
      character(len=120), dimension(nsubdir) :: asubdir
      character(len=250)                     :: filename

!-----------------------------------------------------------------------
! Allocate the ecore and ncore arrays
!-----------------------------------------------------------------------
      if (allocated(ecore)) deallocate(ecore)
      allocate(ecore(maxsta,nsubdir))
      ecore=0.0d0

      if (allocated(ncore)) deallocate(ncore)
      allocate(ncore(nsubdir))
      ncore=0

!-----------------------------------------------------------------------
! Loop over timesteps/subdirectories and read the core-excitation
! energies for each
!-----------------------------------------------------------------------
      unit=20

      do i=1,nsubdir

         ! Read the core-excitation energies if |C_j(t)| > eps
         if (abs(coeff(i)).gt.thrsh) then

            filename=trim(amaindir)//'/'//trim(asubdir(i))//&
                 '/gamess/core_exci.dat'

            open(unit,file=filename,form='formatted',status='old',&
                 iostat=ierr)
            if (ierr.ne.0) goto 999

            ncore(i)=0
10          read(unit,'(25x,F15.10)',end=100) ecore(ncore(i)+1,i)

            ncore(i)=ncore(i)+1
            goto 10
            
100         continue
            
            close(unit)

            cycle
            
            ! If the ORMAS output file was not found, then set the
            ! core excitation energy ludicrously high s.t. this is
            ! picked up in getcontrib_gmsspec
999         continue
            ecore(1,i)=1000.0d0
            
         endif

      enddo

      return
      
    end subroutine getcoreexci_gms

!#######################################################################

    subroutine geteinit(amaindir,asubdir,nsubdir,icontrib,&
         einit,ncore,ista)

      integer                                :: nsubdir,ista,i,n,&
                                                unit
      integer, dimension(nsubdir)            :: ncore
      integer, dimension(nsubdir)            :: icontrib
      real*8, dimension(:), allocatable      :: einit
      character(len=120)                     :: amaindir
      character(len=120), dimension(nsubdir) :: asubdir
      character(len=250)                     :: filename,string
      logical                                :: found

!-----------------------------------------------------------------------
! Allocate the einit array
!-----------------------------------------------------------------------
      if (allocated(einit)) deallocate(einit)
      allocate(einit(nsubdir))
      einit=0.0d0

!-----------------------------------------------------------------------
! Loop over timesteps/subdirectories and read the initial state energy
! for each
!-----------------------------------------------------------------------
      unit=20

      do i=1,nsubdir

         if (icontrib(i).eq.0) cycle

         filename=trim(amaindir)//'/'//trim(asubdir(i))//&
              '/columbus/columbus.out'

         open(unit,file=filename,form='formatted',status='old')
         
         found=.false.
10       read(unit,'(a)',end=100) string
         if (index(string,'- ciudg.x -').eq.0) goto 10
         found=.true.

100      continue

         if (.not.found) then
            write(6,'(/,2x,a,/)') 'MRCI energies not found in the &
                 file '//trim(filename)
            STOP
         endif
         
         n=0
20       read(unit,'(a)') string
         if (index(string,'total mr-sdci energy').eq.0) goto 20
         n=n+1
         if (n.ne.ista) goto 20
         
         read(string,'(34x,F15.10)') einit(i)
         
         close(unit)

      enddo
      
      return

    end subroutine geteinit

!#######################################################################

    subroutine gettdmsq(amaindir,asubdir,nsubdir,icontrib,ncore,tdmsq)

      implicit none

      integer                                :: nsubdir,i,f,k,unit
      integer, dimension(nsubdir)            :: ncore
      integer, dimension(nsubdir)            :: icontrib
      real*8, dimension(:,:), allocatable    :: tdmsq
      real*8, dimension(3)                   :: tdm
      character(len=120)                     :: amaindir
      character(len=120), dimension(nsubdir) :: asubdir
      character(len=250)                     :: filename,string

!-----------------------------------------------------------------------
! Allocate the tdmsq array
!-----------------------------------------------------------------------
      if (allocated(tdmsq)) deallocate(tdmsq)
      allocate(tdmsq(maxsta,nsubdir))
      tdmsq=0.0d0

!-----------------------------------------------------------------------
! Loop over timesteps/subdirectories and read the transition dipole
! moments for each
!-----------------------------------------------------------------------
      unit=20

      do i=1,nsubdir
         
         if (icontrib(i).eq.0) cycle
         
         do f=1,ncore(i)
             
            filename=trim(amaindir)//'/'//trim(asubdir(i))//&
              '/gamess/core'
             k=len_trim(filename)
             if (f-1.lt.10) then
                write(filename(k+1:k+1),'(i1)') f-1
             else
                write(filename(k+1:k+2),'(i2)') f-1
             endif
             filename=trim(filename)//'_tdm.dat'
             
             open(unit,file=filename,form='formatted',status='old')
             
             read(unit,'(20x,3(6x,F15.12))') (tdm(k), k=1,3)

             tdmsq(f,i)=dot_product(tdm,tdm)

             close(unit)

         enddo

      enddo

      return

    end subroutine gettdmsq

!#######################################################################

    subroutine trxas_currifg(nsubdir,icontrib,ncore,ecore,tdmsq,&
         coeff,step,ip,ifg,einit)

      use expec
      use sysdef

      implicit none

      integer                           :: nsubdir,n,k,ifg,i,j
      integer, dimension(nsubdir)       :: icontrib,ncore,step
      real*8, dimension(maxsta,nsubdir) :: ecore,tdmsq
      real*8, dimension(nsubdir)        :: ip,einit
      real*8                            :: tcurr,t,e,dele,delt,csq,&
                                           lfunc,prefac,vee,tcent,&
                                           tsig,musq,shape
      complex*16, dimension(nsubdir)    :: coeff

!-----------------------------------------------------------------------
! Set the grid spacings
!-----------------------------------------------------------------------
      dele=(egrid(2)-egrid(1))/egrid(3)
      delt=(tgrid(2)-tgrid(1))/tgrid(3)

!-----------------------------------------------------------------------
! Calculate the contribution of the current trajectory to the TRXAS
!-----------------------------------------------------------------------
      tsig=fwhm_t/2.35482d0
      
      ! Loop over timesteps/subdirectories
      do n=1,nsubdir
         
         ! Cycle if the current timestep doesn't contribute
         if (icontrib(n).eq.0) cycle

         ! Current time in fs
         tcent=(step(n)-1)*dt/41.341375d0

         ! Loop over final states
         do k=1,ncore(n)

            ! Cycle if the current state lies above the IP
            if (ecore(k,n).ge.ip(n)) cycle

            ! Prefactor
            csq=conjg(coeff(n))*coeff(n)
            prefac=4.0d0*pi*csq/c_au
            ! 'Renormalise'
            prefac=prefac/cnorm(ifg,n)

            ! VEE
            vee=ecore(k,n)-einit(n)            

            ! TDM^2
            musq=tdmsq(k,n)

            ! Loop over grid points
            do i=1,int(egrid(3))
               e=egrid(1)+(i-1)*dele
               call lineshape(shape,gamma,e,vee)
               do j=1,int(tgrid(3))
                  t=tgrid(1)+(j-1)*delt
                  spec(i,j)=spec(i,j)+(1.0d0/real(nifg))&
                       *prefac*vee*musq*shape &
                       * exp(-((t-tcent)**2)/(2.0d0*tsig**2))
               enddo
            enddo

         enddo

      enddo

      return

    end subroutine trxas_currifg

!#######################################################################

    subroutine lineshape(lfunc,gamma,e,deltae)

      implicit none

      real*8 :: lfunc,gamma,e,deltae,numer,denom

      numer=0.5d0*gamma

      denom=(0.25d0*gamma**2) + (deltae-e)**2

      lfunc=numer/denom

      return

    end subroutine lineshape

!#######################################################################

    subroutine wrtrxas
      
      use expec
      use sysdef

      implicit none
      
      integer :: iout,i,j,k,count,negrid,n
      real*8  :: dele,delt,e,t,func,tsig,prefac,tcent

!-----------------------------------------------------------------------
! Set the grid spacings
!-----------------------------------------------------------------------
      dele=(egrid(2)-egrid(1))/egrid(3)
      delt=(tgrid(2)-tgrid(1))/tgrid(3) 

!-----------------------------------------------------------------------
! Open the TR-TXAS output file
!-----------------------------------------------------------------------
      iout=20
      open(iout,file='trxas.dat',form='formatted',status='unknown')

!-----------------------------------------------------------------------
! Calculate and output the TR-TXAS
!-----------------------------------------------------------------------
      do i=1,int(egrid(3))
         write(iout,*)
         do j=1,int(tgrid(3))
            e=egrid(1)+(i-1)*dele
            t=tgrid(1)+(j-1)*delt
            write(iout,*) e*eh2ev,t,spec(i,j)            
         enddo
      enddo
      
      return

    end subroutine wrtrxas
    
!#######################################################################

  end module gamesstrxas

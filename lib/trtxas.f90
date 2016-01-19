  module adctas

    save
    
    integer*8                                     :: nmaindir,nifg,&
                                                     ne,nt,maxfunc,nfunc
    integer*8, dimension(:), allocatable          :: staindx,ifgindx
    real*8, dimension(:,:), allocatable           :: spec,par,cnorm
    real*8, parameter                             :: eh2ev=27.2113845d0
    real*8, parameter                             :: c_au=137.03604d0
    real*8                                        :: pi=3.14159265358979d0
    character(len=120), dimension(:), allocatable :: amaindir

  contains

!#######################################################################
! adc_trtxas: calculates the time-resolved transient X-ray absorption
!             spectrum using ADC cross-section
!#######################################################################

    subroutine adc_trtxas

      use expec

      implicit none

      integer*8 :: i
      
      write(6,'(/,70a)') ('-',i=1,70)
      write(6,'(11x,a)') 'Calculating the TR-TXAS using ADC &
           cross-sections'
      write(6,'(70a)') ('-',i=1,70)

!-----------------------------------------------------------------------
! Read the names of the main directories
!-----------------------------------------------------------------------
      call rdmaindirfile

!-----------------------------------------------------------------------
! Calculate the TR-TXAS
!-----------------------------------------------------------------------
      call calc_trtxas

      return

    end subroutine adc_trtxas

!#######################################################################

    subroutine rdmaindirfile

      use expec
      use parsemod

      implicit none

      integer*8                            :: unit,n
      character(len=120)                   :: string

      integer*8                            :: inkw
      integer*8, parameter                 :: maxkw=60
      integer*8, dimension(maxkw)          :: ilkw
      character(len=120), dimension(maxkw) :: keyword

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

    subroutine calc_trtxas

      use expec
      use trajdef
      use sysdef

      implicit none

      integer*8                                     :: i,nsubdir,itmp,&
                                                       ncontrib,nstates_f
      integer*8, dimension(:), allocatable          :: step,icontrib
      real*8, dimension(:), allocatable             :: ip,einit
      real*8, dimension(:,:), allocatable           :: tdmsq,deltae
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
      call trtxas_alloc

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
! Determine the no. final space Davidson states
!-----------------------------------------------------------------------
      call getsubdirs(amaindir(1),nsubdir,asubdir,step)
      call getcontrib(ncontrib,icontrib,amaindir(1),asubdir,nsubdir)
      call get_nstates_f(nstates_f,amaindir(1),asubdir,nsubdir,icontrib)

!-----------------------------------------------------------------------
! Allocate the TR-TXAS parameter array
!-----------------------------------------------------------------------      
      itmp=(nstep/dstep)+1
      maxfunc=nmaindir*itmp*nstates_f
      allocate(par(maxfunc,4))
      par=0.0d0

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
         call getcontrib(ncontrib,icontrib,amaindir(i),asubdir,nsubdir)

         ! Cycle if no timesteps/subdirectories contribute to the
         ! spectrum
         if (ncontrib.eq.0) cycle

         ! Read the ionisation potential for each timestep/subdirectory
         call getip_adc(ip,amaindir(i),asubdir,nsubdir,icontrib)

         ! Read the excitation energies and TDMs for each 
         ! timestep/subdirectory
         call rddavidson(tdmsq,amaindir(i),asubdir,nsubdir,icontrib,&
              staindx(i),einit,nstates_f,deltae)
         
         ! Calcuate the contribution to the TR-TXAS from the current
         ! trajectory
         call trtxas_currtraj(nsubdir,icontrib,nstates_f,ip,deltae,&
              tdmsq,coeff,step,ifgindx(i))

      enddo

!-----------------------------------------------------------------------
! Calculate and output the total TRPES
!-----------------------------------------------------------------------
      call trtxas_total(nmaindir)

      return

    end subroutine calc_trtxas

!#######################################################################

    subroutine getnifg

      use sysdef

      implicit none

      integer*8                     :: i,k,unit
      integer*8, dimension(nintraj) :: cnt
      character(len=130)            :: ain

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

      integer*8          :: i,unit
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

    subroutine get_nstates_f(nstates_f,amaindir,asubdir,nsubdir,&
         icontrib)

      implicit none
      
      integer*8                              :: nstates_f,nsubdir,i,&
                                                k1,k2,k3,unit,ilbl
      integer*8, dimension(nsubdir)          :: icontrib
      character(len=120)                     :: amaindir,string
      character(len=120), dimension(nsubdir) :: asubdir
      character(len=250)                     :: filename

      unit=20

      do i=1,nsubdir
         
         if (icontrib(i).eq.0) cycle

         k1=index(amaindir,'/')+1
         k2=len_trim(amaindir)
         if (index(asubdir(i),'/').eq.0) then
            k3=len_trim(asubdir(i))
         else
            k3=len_trim(asubdir(i))-1
         endif

         filename=trim(amaindir)//'/'//trim(asubdir(i)) &
              //'/adc_dav_x/adc_'//amaindir(k1:k2)//'_' &
              //asubdir(i)(1:k3)//'_dav_x.inp'

         open(unit,file=filename,form='formatted',status='old')

10       read(unit,'(a)') string
         if (index(string,'diag_final_section').eq.0) goto 10

15       read(unit,'(a)') string
         if (index(string,'nstates').eq.0) goto 15
         
         k1=index(string,'=')+1
         k2=len_trim(string)
         read(string(k1:k2),*) nstates_f
         
         close(unit)
         
         exit

      enddo

      return

    end subroutine get_nstates_f

!#######################################################################

    subroutine trtxas_alloc

      use expec
      
      implicit none
      
      integer*8 :: ne,nt
      
      ! Electronic state indices
      allocate(staindx(nmaindir))

      ! TR-TXAS
      ne=int(egrid(3))
      nt=int(tgrid(3))
      allocate(spec(ne,nt))

      return

    end subroutine trtxas_alloc

!#######################################################################

    subroutine getcnorm

      use sysdef 
      use expec

      implicit none

      integer*8                                     :: n,k,ifg,itmp,&
                                                       nsubdir,&
                                                       ncontrib
      integer*8, dimension(:), allocatable          :: step,icontrib
      real*8                                        :: ftmp
      complex*16, dimension(:), allocatable         :: coeff
      character(len=120), dimension(:), allocatable :: asubdir

      itmp=nstep/dstep+1
      allocate(cnorm(nifg,itmp))
      cnorm=0.0d0

      ! Loop over main directories/trajectories
      do n=1,nmaindir
         
         ! IFG no. for the current trajectory
         ifg=ifgindx(n)
         
         ! Get the list of timesteps/subdirectories for the
         ! current trajectory
         call getsubdirs(amaindir(n),nsubdir,asubdir,step)
         
         ! Read the coefficients for each timestep/subdirectory
         call getcoeff(coeff,amaindir(n),asubdir,nsubdir)

         ! Determine which timesteps/subdirectories contribute to the
         ! spectrum
         call getcontrib(ncontrib,icontrib,amaindir(n),asubdir,nsubdir)
         
         ! Add the current contribution to the cnorm array
         do k=1,nsubdir            
            if (icontrib(k).eq.0) cycle

            itmp=step(k)/dstep+1
            cnorm(ifg,itmp)=cnorm(ifg,itmp)+conjg(coeff(k))*coeff(k)
         enddo


      enddo

      cnorm=sqrt(cnorm)

!      do n=1,nifg
!         do k=1,50
!            print*,cnorm(n,k)
!         enddo
!      enddo

      return

    end subroutine getcnorm

!#######################################################################

    subroutine getsubdirs(amaindir,nsubdir,asubdir,step)

      implicit none

      integer*8                                     :: unit,nsubdir,i
      integer*8, dimension(:), allocatable          :: step
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

      integer*8                              :: nsubdir,i,unit
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

    subroutine getcontrib(ncontrib,icontrib,amaindir,asubdir,nsubdir)

      implicit none

      integer*8                              :: ncontrib,nsubdir,i
      integer*8, dimension(:), allocatable   :: icontrib
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
! contributes to the TR-TXAS
!-----------------------------------------------------------------------
      do i=1,nsubdir

         ! IP-ADC calculation
         filename=trim(amaindir)//'/'//trim(asubdir(i))//'/adc_ip/davstates.dat'
         inquire(file=filename,exist=exists)
         if (.not.exists) goto 10
         
         ! ADC, Davidson, x
         filename=trim(amaindir)//'/'//trim(asubdir(i))//'/adc_dav_x/osc.dat'
         inquire(file=filename,exist=exists)
         if (.not.exists) goto 10

         ! ADC, Davidson, y
         filename=trim(amaindir)//'/'//trim(asubdir(i))//'/adc_dav_y/osc.dat'
         inquire(file=filename,exist=exists)
         if (.not.exists) goto 10
         
         ! ADC, Davidson, z
         filename=trim(amaindir)//'/'//trim(asubdir(i))//'/adc_dav_z/osc.dat'
         inquire(file=filename,exist=exists)
         if (.not.exists) goto 10

         ncontrib=ncontrib+1

         cycle

10       continue
         icontrib(i)=0

      enddo

      return

    end subroutine getcontrib

!#######################################################################

    subroutine getip_adc(ip,amaindir,asubdir,nsubdir,icontrib)

      implicit none

      integer*8                              :: nsubdir,unit,i
      integer*8, dimension(nsubdir)          :: icontrib
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

         filename=trim(amaindir)//'/'//trim(asubdir(i)) &
              //'/adc_ip/davstates.dat'
         open(unit,file=filename,form='formatted',status='old')

5        read(unit,'(a)') string
         if (index(string,'MP2 energy:').eq.0) goto 5
         read(string,'(28x,F14.8)') e0

10       read(unit,'(a)') string
         if (index(string,'Energy:').eq.0) goto 10
         read(string,'(27x,F14.8)') ip(i)
         ip(i)=ip(i)-e0

         close(unit)
      enddo

      return

    end subroutine getip_adc

!#######################################################################

    subroutine rddavidson(tdmsq,amaindir,asubdir,nsubdir,icontrib,&
         staindx,einit,nstates_f,deltae)
      
      implicit none

      integer*8                              :: nsubdir,staindx,unit,&
                                                i,c,count,k1,k2,k3,&
                                                adcstate,nstates_f,&
                                                indx
      integer*8, dimension(nsubdir)          :: icontrib
      real*8, dimension(:), allocatable      :: einit
      real*8, dimension(:,:), allocatable    :: tdmsq,deltae
      real*8                                 :: ener,osc
      complex*16, dimension(nsubdir)         :: coeff
      character(len=120)                     :: amaindir,string
      character(len=120), dimension(nsubdir) :: asubdir
      character(len=250)                     :: filename
      character(len=1), dimension(3)         :: dpllbl

!-----------------------------------------------------------------------
! Initialisation of various things
!-----------------------------------------------------------------------
      unit=20
      dpllbl=(/ 'x','y','z' /)

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      if (allocated(tdmsq)) deallocate(tdmsq)
      allocate(tdmsq(nsubdir,nstates_f))
      tdmsq=0.0d0

      if (allocated(deltae)) deallocate(deltae)
      allocate(deltae(nsubdir,nstates_f))
      deltae=0.0d0

      if (allocated(einit)) deallocate(einit)
      allocate(einit(nsubdir))
      einit=0.0d0

!-----------------------------------------------------------------------
! Read in the energy of the initial state
!-----------------------------------------------------------------------      
      ! Loop over timesteps/subdirectories
      do i=1,nsubdir
         
         if (icontrib(i).eq.0) cycle

         ! (1) Read the index of the ADC sate chosen by the target
         !     matching routine
         k1=index(amaindir,'/')+1
         k2=len_trim(amaindir)
         if (index(asubdir(i),'/').eq.0) then
            k3=len_trim(asubdir(i))
         else
            k3=len_trim(asubdir(i))-1
         endif

         filename=trim(amaindir)//'/'//trim(asubdir(i)) &
              //'/adc_dav_x/adc_'//amaindir(k1:k2)//'_' &
              //asubdir(i)(1:k3)//'_dav_x.log'

         open(unit,file=filename,form='formatted',status='old')

5        read(unit,'(a)') string
         if (index(string,'selected').eq.0) goto 5
         
         read(string,'(7x,i2)') adcstate
         
         close(unit)

         ! (2) Read the energy of the initial state
         filename=trim(amaindir)//'/'//trim(asubdir(i)) &
              //'/adc_dav_x/davstates.dat'

         open(unit,file=filename,form='formatted',status='old')

         if (adcstate.eq.0) then
            ! Ground state
10          read(unit,'(a)') string
            if (index(string,'Ground state MP2 energy:').eq.0) goto 10
            read(string,'(28x,F14.8)') einit(i)
         else
            ! Excited states
            count=0
15          read(unit,'(a)') string
            if (index(string,'Energy:').eq.0) then
               goto 15
            else
               count=count+1
               if (count.eq.adcstate) then
                  read(string,'(27x,F14.8)') einit(i)
               else
                  goto 15
               endif
            endif
         endif
         
         close(unit)

      enddo

!-----------------------------------------------------------------------
! Read in the TDM^2 values
!-----------------------------------------------------------------------
      ! Loop over timesteps/subdirectories
      do i=1,nsubdir

         if (icontrib(i).eq.0) cycle

         ! (1) Excitation energies
         k1=index(amaindir,'/')+1
         k2=len_trim(amaindir)
         if (index(asubdir(i),'/').eq.0) then
            k3=len_trim(asubdir(i))
         else
            k3=len_trim(asubdir(i))-1
         endif
         
         filename=trim(amaindir)//'/'//trim(asubdir(i)) &
              //'/adc_dav_x/adc_'//amaindir(k1:k2)//'_' &
              //asubdir(i)(1:k3)//'_dav_x.log'

         open(unit,file=filename,form='formatted',status='old')

         count=0
20       read(unit,'(a)') string
         if (index(string,'Final space CVS-ADC(2)').eq.0) goto 20

25       read(unit,'(a)') string
         if (index(string,'Energy:').eq.0) then
            goto 25
         else
            count=count+1
            read(string,'(27x,F14.8)') ener
            deltae(i,count)=ener-einit(i)
            if (count.lt.nstates_f) goto 25
         endif

         close(unit)         

         ! (2) TDM^2 values
         do c=1,3
                        
            filename=trim(amaindir)//'/'//trim(asubdir(i)) &
              //'/adc_dav_'//dpllbl(c)//'/osc.dat'
            
            open(unit,file=filename,form='formatted',status='old')

30          read(unit,*,end=35) ener,osc

            call get_state(indx,ener,deltae(i,:),nstates_f)

            tdmsq(i,indx)=tdmsq(i,indx)+(3.0d0/(2.0d0*deltae(i,indx)))*osc

            goto 30

35          continue

            close(unit)
            
         enddo

      enddo

      return

    end subroutine rddavidson

!#######################################################################

    subroutine get_state(indx,ener,deltae,nstates_f)
      
      implicit none

      integer*8                    :: indx,nstates_f,i
      real*8                       :: ener,ftmp
      real*8, dimension(nstates_f) :: deltae
      real*8, parameter            :: tol=0.0001d0

      indx=-1

      do i=1,nstates_f
         if (abs(ener/eh2ev-deltae(i)).le.tol) then
            indx=i
            exit
         endif
      enddo

      if (indx.eq.-1) then
         write(6,'(/,2x,a,/)') 'Error in subroutine get_state in &
              trtxas.f90'
         stop
      endif

      return

    end subroutine get_state

!#######################################################################

    subroutine trtxas_currtraj(nsubdir,icontrib,nstates_f,ip,deltae,&
         tdmsq,coeff,step,ifg)

      use expec
      use sysdef

      implicit none

      integer*8                            :: nsubdir,nstates_f,n,k,&
                                              ifg
      integer*8, dimension(nsubdir)        :: icontrib,step
      real*8, dimension(nsubdir)           :: ip
      real*8, dimension(nsubdir,nstates_f) :: deltae,tdmsq
      real*8                               :: tcurr,t,e,dele,delt,csq,&
                                              lfunc,prefac
      complex*16, dimension(nsubdir)       :: coeff
      
      ! Temporary: Gamma = 0.5 ev
      real*8, parameter :: gamma=13.60569225d0

      ! Loop over timesteps
      do n=1,nsubdir
        
         ! Cycle if the current timestep doesn't contribute
         if (icontrib(n).eq.0) cycle
 
         ! Current time in fs
         t=(step(n)-1)*dt/41.341375d0
        
         ! Loop over final states
         do k=1,nstates_f
            
            ! Cycle if the current state lies above the IP
            if (deltae(n,k).ge.ip(n)) cycle

            ! Keep track of the no. functions that will form
            ! the final TR-TXAS
            nfunc=nfunc+1

            ! Fill in the parameter array for the current 
            ! function:
            !
            ! (1) Prefactor
            csq=conjg(coeff(n))*coeff(n)
            par(nfunc,1)=4.0d0*pi*csq/c_au
            ! 'Renormalise'
            par(nfunc,1)=par(nfunc,1)/cnorm(ifg,n)

            ! (2) TDM^2
            par(nfunc,2)=tdmsq(n,k)

            ! (3) Delta E_IF
            par(nfunc,3)=deltae(n,k)

            ! (4) Centre wrt t
            par(nfunc,4)=t

         enddo

      enddo

      return

    end subroutine trtxas_currtraj

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

    subroutine trtxas_total(nmaindir)

      use expec
      use sysdef

      implicit none

      integer*8 :: nmaindir,iout,i,j,k
      real*8    :: dele,delt,e,t,func,prefac,musq,shape,deif,tcent

      ! Temporary: (i)  Gamma = 0.5 ev
      !            (ii) Pump-probe cross-corr. = 50 fs
      real*8, parameter :: gamma=0.5d0/eh2ev
      real*8, parameter :: tsig=50.0d0/2.35482d0

      write(6,'(/,2x,a,/)') 'Constructing the TR-TXAS...'

!-----------------------------------------------------------------------
! Set the grid spacings
!-----------------------------------------------------------------------
      dele=(egrid(2)-egrid(1))/egrid(3)
      delt=(tgrid(2)-tgrid(1))/tgrid(3) 

!-----------------------------------------------------------------------
! Open the TR-TXAS output file
!-----------------------------------------------------------------------
      iout=20
      open(iout,file='trtxas.dat',form='formatted',status='unknown')

!-----------------------------------------------------------------------
! Calculate and output the TR-TXAS
!-----------------------------------------------------------------------
      ! Loop over grid points
      do i=1,int(egrid(3))+1
         write(iout,*)
         do j=1,int(tgrid(3))+1

            e=egrid(1)+(i-1)*dele
            t=tgrid(1)+(j-1)*delt

            ! Loop over the individual functions, and sum the 
            ! contributions from each
            func=0.0d0
            do k=1,nfunc
               
               ! Prefactor
               prefac=par(k,1)

               ! TDM^2
               musq=par(k,2)

               ! Delta E_IF
               deif=par(k,3)

               ! Centre wrt t
               tcent=par(k,4)

               ! Lorentzian lineshape
               call lineshape(shape,gamma,e,deif)

               ! Function value
               func=func+(1.0d0/real(nifg))*prefac*e*musq*shape &
                    * exp(-((t-tcent)**2)/(2.0d0*tsig**2))

            enddo

            ! Ouput the current photon energy, time delay and
            ! function value
            write(iout,*) e*eh2ev,t,func

         enddo
      enddo

!-----------------------------------------------------------------------
! Close the TR-TXAS output file
!-----------------------------------------------------------------------
      close(iout)

      return

    end subroutine trtxas_total

!#######################################################################

  end module adctas

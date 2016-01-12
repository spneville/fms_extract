  module trpes
    
    implicit none

  contains

!#######################################################################
! trpes_dnorm: calculates the TRPES using Dyson orbital norms
!#######################################################################

    subroutine trpes_dnorm

      use expec

      implicit none

      integer*8                                     :: i,nmaindir
      character(len=120), dimension(:), allocatable :: amaindir

      write(6,'(/,70a)') ('-',i=1,70)
      write(6,'(11x,a)') 'Calculating the TRPES using Dyson orbital norms'
      write(6,'(70a)') ('-',i=1,70)

!-----------------------------------------------------------------------
! Read the names of the main directories
!-----------------------------------------------------------------------
      call rdmaindirfile(nmaindir,amaindir)

!-----------------------------------------------------------------------
! Calculate the TRPES
!-----------------------------------------------------------------------
      call calc_trpes(nmaindir,amaindir)
      
      return

    end subroutine trpes_dnorm

!#######################################################################

    subroutine rdmaindirfile(nmaindir,amaindir)

      use expec
      use parsemod

      implicit none

      integer*8                                     :: nmaindir,unit,n
      character(len=120), dimension(:), allocatable :: amaindir
      character(len=120)                            :: string

      integer*8                                     :: inkw
      integer*8, parameter                          :: maxkw=60
      integer*8, dimension(maxkw)                   :: ilkw
      character(len=120), dimension(maxkw)          :: keyword

!-----------------------------------------------------------------------
! Open file
!-----------------------------------------------------------------------
      unit=20
      open(unit,file=dnormfile,form='formatted',status='old')

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

    subroutine calc_trpes(nmaindir,amaindir)

      use expec
      use trajdef
      use sysdef

      implicit none

      integer*8                                     :: nmaindir,i,nifg,&
                                                       nsubdir,maxgauss
      integer*8, dimension(nmaindir)                :: staindx
      integer*8, dimension(:), allocatable          :: step
      real*8, dimension(:,:), allocatable           :: ip,dnorm
      complex*16, dimension(:), allocatable         :: coeff
      character(len=120), dimension(nmaindir)       :: amaindir
      character(len=120), dimension(:), allocatable :: asubdir

!-----------------------------------------------------------------------
! Determine the number of IFGs being considered, which may be less than
! nintraj
!-----------------------------------------------------------------------
      call getnifg(nifg,nmaindir,amaindir)

!-----------------------------------------------------------------------
! Allocate the Gaussian parameter array gausspar
!-----------------------------------------------------------------------
      maxgauss=nmaindir*nstep*nionize
      allocate(gausspar(maxgauss,4))
      gausspar=0.0d0

!-----------------------------------------------------------------------
! Determine the state indices for each trajectory (main directory) being
! considered
!-----------------------------------------------------------------------
      call getstaindx(staindx,nifg,nmaindir,amaindir)

!-----------------------------------------------------------------------
! Loop over the main trajectories (main directories), reading the
! state energies and dyson orbital norms for each timestep
!-----------------------------------------------------------------------
      ngauss=0
      do i=1,nmaindir

         write(6,'(2a)') 'Processing directory: ',trim(amaindir(i))

         ! Get the list of timesteps/subdirectories
         call getsubdirs(amaindir(i),nsubdir,asubdir,step)

         ! Read the coefficients for each timestep/subdirectory
         call getcoeff(coeff,amaindir(i),asubdir,nsubdir)

         ! Determine the ionization energies for each
         ! timestep/subdirectory
         call getip(ip,amaindir(i),asubdir,nsubdir,coeff)

         ! Read the Dyson orbital norms for each ionization
         ! channel for each timestep/subdirectory
         call getdnorm(dnorm,amaindir(i),asubdir,nsubdir,step,coeff)

         ! Calcuate the contribution to the TRPES from the current
         ! trajectory
         call trpes_currtraj(ip,dnorm,nsubdir,step,coeff,staindx(i),&
              nifg)

      enddo

!-----------------------------------------------------------------------
! Calculate and output the total TRPES
!-----------------------------------------------------------------------
      call trpes_total(nmaindir)

      return

    end subroutine calc_trpes

!#######################################################################

    subroutine getnifg(nifg,nmaindir,amaindir)

      use sysdef

      implicit none

      integer*8                               :: nifg,nmaindir,i,k,unit
      integer*8, dimension(nintraj)           :: cnt
      character(len=120), dimension(nmaindir) :: amaindir
      character(len=130)                      :: ain

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
         close(unit)
      enddo

      nifg=0
      do i=1,nintraj
         if (cnt(i).gt.0) nifg=nifg+1
      enddo
      
    end subroutine getnifg

!#######################################################################

    subroutine getstaindx(staindx,nifg,nmaindir,amaindir)

      implicit none

      integer*8, dimension(nmaindir)          :: staindx
      integer*8                               :: nifg,nmaindir,i,unit
      character(len=120), dimension(nmaindir) :: amaindir
      character(len=130)                      :: ain

      do i=1,nmaindir         
         ain=trim(amaindir(i))//'/state_id'
         open(unit,file=ain,form='formatted',status='old')
         read(unit,*) staindx(i)
      enddo

    end subroutine getstaindx

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

    end subroutine getsubdirs

!#######################################################################

    subroutine getip(ip,amaindir,asubdir,nsubdir,coeff)

      use expec

      implicit none

      integer*8                              :: nsubdir,i,j,k,unit,nsta,&
                                                sn,sc
      real*8, dimension(:,:), allocatable    :: ip
      real*8, dimension(30)                  :: en_n,en_c
      real*8                                 :: csq
      complex*16, dimension(nsubdir)         :: coeff
      character(len=120)                     :: amaindir
      character(len=120), dimension(nsubdir) :: asubdir
      character(len=270)                     :: acolout
      character(len=120)                     :: string

!-----------------------------------------------------------------------
! Allocate the ip array
!-----------------------------------------------------------------------
      if (allocated(ip)) deallocate(ip)
      allocate(ip(nsubdir,nionize))
      ip=0.0d0

!-----------------------------------------------------------------------
! Read the energies of the neutral and cationic states from the
! Columbus output files
!-----------------------------------------------------------------------
      en_n=0.0d0
      en_c=0.0d0
      unit=20

      ! Loop over timesteps/subdirectories
      do i=1,nsubdir

         ! Read the energies of the neutral states
         acolout=trim(amaindir)//'/'//trim(asubdir(i))//&
              '/mrci_neutral/columbus.out'
         open(unit,file=acolout,form='formatted',status='old')
10       continue
         read(unit,'(a)',end=999) string
         if (string(21:31).ne.'- ciudg.x -') goto 10
         nsta=0
15       continue
         read(unit,'(a)') string
         if (string(2:4).ne.'---') then
            if (string(2:29).eq.'mr-sdci+(1-c0**2)*(eci-eref)') then
               nsta=nsta+1
               read(string,'(46x,F18.12)') en_n(nsta)
            endif
            goto 15
         endif
         close(unit)

         ! Read the energies of the cationic states
         acolout=trim(amaindir)//'/'//trim(asubdir(i))//&
              '/mrci_cation/columbus.out'
         open(unit,file=acolout,form='formatted',status='old')
20       continue
         read(unit,'(a)',end=999) string
         if (string(21:31).ne.'- ciudg.x -') goto 20
         nsta=0
25       continue
         read(unit,'(a)') string
         if (string(2:4).ne.'---') then
            if (string(2:29).eq.'mr-sdci+(1-c0**2)*(eci-eref)') then
               nsta=nsta+1
               read(string,'(46x,F18.12)') en_c(nsta)
            endif
            goto 25
         endif
         close(unit)

         ! Fill in the ionization potential array ip s.t. its elements
         ! correspond to the user requested ionization channels
         do j=1,nionize            
            sn=iionize(j,1)
            sc=iionize(j,2)
            ip(i,j)=(en_c(sc)-en_n(sn))*27.211d0+ipshift(j)            
         enddo
   
         cycle
      
         ! Die here if we have not found a set of mrci energies
999      continue
         csq=conjg(coeff(i))*coeff(i)
         if (csq.gt.1d-4) then
            write(6,'(/,2x,2(a,x),/)') &
                 'mrci energies not found in the file:',trim(acolout)
            STOP
         endif

      enddo

      return

!-----------------------------------------------------------------------
! Die here if we have not found a set of mrci energies
!-----------------------------------------------------------------------
!999   continue
!      csq=conjg(coeff(i))*coeff(i)
!      if (csq.gt.1d-4) then
!         write(6,'(/,2x,2(a,x),/)') &
!              'mrci energies not found in the file:',trim(acolout)
!         STOP
!      endif

    end subroutine getip

!#######################################################################

    subroutine getdnorm(dnorm,amaindir,asubdir,nsubdir,step,coeff)

      use expec

      implicit none

      integer*8                              :: nsubdir,i,j,k,unit,si,&
                                                dj,n
      integer*8, dimension(nsubdir)          :: step
      real*8, dimension(:,:), allocatable    :: dnorm
      real*8                                 :: csq
      complex*16, dimension(nsubdir)         :: coeff
      character(len=120)                     :: amaindir
      character(len=120), dimension(nsubdir) :: asubdir
      character(len=270)                     :: asdout
      character(len=120)                     :: string
      character(len=2)                       :: asi,adj
      character(len=9)                       :: astep

!-----------------------------------------------------------------------
! Allocate the dnorm array
!-----------------------------------------------------------------------
      if (allocated(dnorm)) deallocate(dnorm)
      allocate(dnorm(nsubdir,nionize))
      dnorm=0.0d0

!-----------------------------------------------------------------------
! Read the Dyson orbital norms for each ionization channel for each
! timestep
!-----------------------------------------------------------------------
      unit=20

      ! Loop over timesteps/subdirectories
      do i=1,nsubdir

         ! Loop over ionization channels
         do j=1,nionize

            ! Open the current superdyson output file
            si=iionize(j,1)-1
            dj=iionize(j,2)-1            
            n=step(i)
            if (n.lt.10) then
               write(astep,'(a8,i1)') 'step0000',n
            else if (n.lt.100) then
               write(astep,'(a7,i2)') 'step000',n
            else if (n.lt.1000) then
               write(astep,'(a6,i3)') 'step00',n
            else if (n.lt.10000) then
               write(astep,'(a5,i4)') 'step0',n
            else
               write(astep,'(a4,i5)') 'step',n
            endif
            write(asi,'(a1,i1)') 'S',si
            write(adj,'(a1,i1)') 'D',dj
            asdout=trim(amaindir)//'/'//trim(astep)//&
                 '_'//asi//'_'//adj//'.log'
            open(unit,file=asdout,form='formatted',status='old')

            ! Read the Dyson orbital norm from the current superdyson
            ! output file
10          read(unit,'(a)',end=999) string
            if (string(2:5).ne.'psid') goto 10
            read(string,'(18x,F10.8)') dnorm(i,j)
            dnorm(i,j)=sqrt(dnorm(i,j))

            ! Close the current superdyson output file
            close(unit)

         enddo

         cycle

         ! Die here if we have not found a Dyson orbital norm 
999      continue
         csq=conjg(coeff(i))*coeff(i)
         if (csq.gt.1d-4) then
            write(6,'(/,2x,2(a,x),/)') &
                 'Dyson orbital norm not found in the file:',trim(asdout)
            STOP
         endif

      enddo

      return

!-----------------------------------------------------------------------
! Die here if we have not found a Dyson orbital norm
!-----------------------------------------------------------------------
!999   continue
!      write(6,'(/,2x,2(a,x),/)') &
!           'Dyson orbital norm not found in the file:',trim(asdout)
!      STOP

    end subroutine getdnorm

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

    subroutine trpes_currtraj(ip,dnorm,nsubdir,step,coeff,ista,nifg)

      use expec
      use sysdef

      implicit none

      integer*8                          :: nsubdir,n,k,ista,nifg
      integer*8, dimension(nsubdir)      :: step
      real*8, dimension(nsubdir,nionize) :: ip,dnorm
      real*8                             :: t,csq
      real*8, parameter                  :: cthresh=0.01d0
      complex*16, dimension(nsubdir)     :: coeff
      logical(kind=4)                    :: linitsta

      ! Loop over timesteps
      do n=1,nsubdir

         ! Current time in fs
         t=(step(n)-1)*dt/41.341375d0

         ! Loop over ionization channels
         do k=1,nionize

            ! Check whether the current trajectories associated
            ! electronic state corresponds to the initial state
            ! of the current ionization channel
            if (iionize(k,1).ne.ista) cycle

            ! Ignore the current contribution if the vertical
            ! (shifted) IP is greater than the probe energy
            if (ip(n,k).gt.eprobe) cycle
            
            ! Ignore the current contribution if the square modulus
            ! of the coefficient is below threshold
            csq=conjg(coeff(n))*coeff(n)
            if (csq.lt.cthresh) cycle

            ! Keep track of the no. two-dimensional Gaussians
            ! that will form our TRPES
            ngauss=ngauss+1

            ! Fill in the Gaussian parameter array for the
            ! current Gaussian:
            !
            ! (1) Coefficient
            csq=conjg(coeff(n))*coeff(n)
            gausspar(ngauss,1)=csq*dnorm(n,k)/nifg

            ! (2) Width parameter, energy domain
            gausspar(ngauss,2)=fwhm_e/2.35482d0
            
            ! (3) Width parameter, time domain
            gausspar(ngauss,3)=fwhm_t/2.35482d0

            ! (4) Centre wrt E
            gausspar(ngauss,4)=eprobe-ip(n,k)

            ! (5) Centre wrt t
            gausspar(ngauss,5)=t

         enddo

      enddo

      return

    end subroutine trpes_currtraj

!#######################################################################

    subroutine trpes_total(nifg)

      use expec
      use sysdef

      implicit none

      integer*8 :: i,j,k,iout,nifg
      real*8    :: a,esig,tsig,tcent,ecent,e,t,dele,delt,func

      write(6,'(/,2x,a,/)') 'Constructing the TRPES...'

!-----------------------------------------------------------------------
! Set the grid spacings
!-----------------------------------------------------------------------
      dele=(egrid(2)-egrid(1))/egrid(3)
      delt=(tgrid(2)-tgrid(1))/tgrid(3)

!-----------------------------------------------------------------------
! Open the TRPES output file
!-----------------------------------------------------------------------
      iout=20
      open(iout,file='trpes.dat',form='formatted',status='unknown')

!-----------------------------------------------------------------------
! Calculate and output the TRPES
!-----------------------------------------------------------------------
      ! Loop over grid points
      do i=1,int(egrid(3))+1
         write(iout,*)
         do j=1,int(tgrid(3))+1      
            ! Set the current energy and time
!            e=(i-1)*dele
!            t=(j-1)*delt

            e=egrid(1)+(i-1)*dele
            t=tgrid(1)+(j-1)*delt

            ! Loop over Gaussians, summing the contributions from each
            func=0.0d0
            do k=1,ngauss
               ! Set the current parameters
               a=gausspar(k,1)
               esig=gausspar(k,2)
               tsig=gausspar(k,3)
               ecent=gausspar(k,4)
               tcent=gausspar(k,5)
               ! Calculate the contribution of the current Gaussian
               func=func+a*exp(-((e-ecent)**2)/(2.0d0*esig**2))&
                    *exp(-((t-tcent)**2)/(2.0d0*tsig**2))
            enddo
            ! Ouput the current photoelectron energy, time delay and
            ! function value
            write(iout,*) e,t,func
         enddo
      enddo

!-----------------------------------------------------------------------
! Close the TRPES output file
!-----------------------------------------------------------------------
      close(iout)

!-----------------------------------------------------------------------
! Write the gnuplot file
!-----------------------------------------------------------------------
      open(iout,file='trpes_dnorm.gnu',form='formatted',status='unknown')
      write(iout,'(a)') '# ~/.gnuplot'
      write(iout,'(a,/)') 'set palette @MATLAB'
      write(iout,'(a)') 'set pm3d map interpolate 0,0'
      write(iout,'(a)') 'set xlabel ''Time (fs)'''
      write(iout,'(a)') 'set ylabel ''E (eV)'''
      write(iout,'(a)') 'splot ''trpes.dat'''
      write(iout,'(a)') 'pause -1'
      close(iout)

      return

    end subroutine trpes_total

!#######################################################################

  end module trpes

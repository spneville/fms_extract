  module adctas

    save
    
    integer                                       :: nmaindir,nifg,&
                                                     ne,nt,maxfunc,nfunc
    integer, dimension(:), allocatable            :: staindx,ifgindx
    integer                                       :: failunit,okunit,cifunit
    integer                                       :: maxion
    real*8, dimension(:,:), allocatable           :: spec,par,cnorm
    real*8, parameter                             :: eh2ev=27.2113845d0
    real*8, parameter                             :: c_au=137.03604d0
    real*8                                        :: pi=3.14159265358979d0
    character(len=120), dimension(:), allocatable :: amaindir
    logical                                       :: lcont,lbound,ldyson

  contains

!#######################################################################
! adc_trtxas: calculates the time-resolved transient X-ray absorption
!             spectrum using ADC cross-sections
!#######################################################################

    subroutine adc_trtxas

      use expec

      implicit none

      integer :: i
      
      write(6,'(/,70a)') ('-',i=1,70)
      write(6,'(11x,a)') 'Calculating the TR-TXAS using ADC &
           cross-sections'
      write(6,'(70a)') ('-',i=1,70)

!-----------------------------------------------------------------------
! Set job type flags
!-----------------------------------------------------------------------
      lbound=.false.
      lcont=.false.
      ldyson=.false.
      if (ijob.eq.8) then
         lbound=.true.
      else if (ijob.eq.9) then
         lcont=.true.
      else if (ijob.eq.11) then
         ldyson=.true.
      endif

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
      use projmod

      implicit none

      integer                                       :: i,nsubdir,itmp,&
                                                       ncontrib,nstates_f
      integer, dimension(:), allocatable            :: step,icontrib
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
! CI filtering
!-----------------------------------------------------------------------
      if (lcifilter) then
         call rdseamfiles2
         call calc_projector2
      endif

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
      if (lbound) then
         call getsubdirs(amaindir(1),nsubdir,asubdir,step)
         call getcontrib(ncontrib,icontrib,amaindir(1),asubdir,&
              nsubdir,0,staindx(1))
         call get_nstates_f(nstates_f,amaindir(1),asubdir,nsubdir,&
              icontrib)
      endif

!-----------------------------------------------------------------------
! Allocate the TR-TXAS parameter array
!-----------------------------------------------------------------------      
      if (lbound) then
         itmp=(nstep/dstep)+1
         maxfunc=nmaindir*itmp*nstates_f
         allocate(par(maxfunc,4))
         par=0.0d0
      else if (lcont) then
         ! Do nothing as we directly fill in the spectrum array
         ! for this case
      else if (ldyson) then
         ! TEMPORARY BODGE: the maximum number of ionization channels
         ! needs to be determined from, e.g., the input vile. For now 
         ! we set this to 200.
         maxion=200
         itmp=(nstep/dstep)+1
         maxfunc=nmaindir*itmp*maxion
         allocate(par(maxfunc,5))
         par=0.0d0
      endif

!-----------------------------------------------------------------------
! Loop over the main trajectories (main directories), reading the
! TDMs, energy differences and coefficients for each timestep
!-----------------------------------------------------------------------
      nfunc=0
      do i=1,nmaindir

!         ! BODGE
!         if (staindx(i).eq.1) cycle
!         ! BODGE

         write(6,'(2a)') 'Processing directory: ',trim(amaindir(i))

         ! Get the list of timesteps/subdirectories
         call getsubdirs(amaindir(i),nsubdir,asubdir,step)

         ! Read the coefficients for each timestep/subdirectory
         call getcoeff(coeff,amaindir(i),asubdir,nsubdir)

         ! Determine which timesteps/subdirectories contribute to the
         ! spectrum
         call getcontrib(ncontrib,icontrib,amaindir(i),asubdir,&
              nsubdir,0,staindx(i))

         ! Cycle if no timesteps/subdirectories contribute to the
         ! spectrum
         if (ncontrib.eq.0) cycle

         ! Read the ionisation potential for each timestep/subdirectory
         call getip_adc(ip,amaindir(i),asubdir,nsubdir,icontrib)

         ! Determine the cross-sections for each timestep/subdirectory
         call getxsec(tdmsq,amaindir(i),asubdir,nsubdir,icontrib,&
              staindx(i),einit,nstates_f,deltae,ip,coeff,ifgindx(i),&
              step)
         
         ! Calcuate the contribution to the bound part of the TR-TXAS 
         ! from the current trajectory
         if (lbound) call trtxas_currtraj(nsubdir,icontrib,nstates_f,&
              ip,deltae,tdmsq,coeff,step,ifgindx(i))

      enddo

!-----------------------------------------------------------------------
! Calculate and output the total spectrum
!-----------------------------------------------------------------------
      call trtxas_total

      return

    end subroutine calc_trtxas

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

    subroutine get_nstates_f(nstates_f,amaindir,asubdir,nsubdir,&
         icontrib)

      implicit none
      
      integer                                :: nstates_f,nsubdir,i,j,&
                                                k1,k2,k3,unit,ilbl,istart,&
                                                iend
      integer, dimension(nsubdir)            :: icontrib
      character(len=120)                     :: amaindir,string
      character(len=120), dimension(nsubdir) :: asubdir
      character(len=250)                     :: filename

      unit=20

      do i=1,nsubdir
         
         if (icontrib(i).eq.0) cycle

         iend=len_trim(amaindir)
         istart=0
         do j=3,iend
            if (amaindir(j-2:j).eq.'../') istart=j+1
         enddo
         if (istart.eq.0) istart=1
         k1=index(amaindir(istart:iend),'/')+istart 
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
      
      integer :: ne,nt
      
      ! Electronic state indices
      allocate(staindx(nmaindir))

      ! TR-TXAS
      ne=int(egrid(3))
      nt=int(tgrid(3))
      allocate(spec(ne,nt))
      spec=0.0d0

      return

    end subroutine trtxas_alloc

!#######################################################################

    subroutine getcnorm

      use sysdef 
      use expec
      use iomod

      implicit none

      integer                                       :: n,k,ifg,itmp,&
                                                       nsubdir,&
                                                       ncontrib
      integer, dimension(:), allocatable            :: step,icontrib
      real*8                                        :: ftmp
      complex*16, dimension(:), allocatable         :: coeff
      character(len=120), dimension(:), allocatable :: asubdir

      ! Open the files containing the geometries at which the
      ! target matching code failed/succeeded, and, if required,
      ! the geometries filtered out due to their proximity to
      ! an interesection seam
      failunit=55
      open(failunit,file='targfail.xyz',form='formatted',status='unknown')
      okunit=56
      open(okunit,file='targok.xyz',form='formatted',status='unknown')
      if (lcifilter) then
         cifunit=57
         open(cifunit,file='cifilter.xyz',form='formatted',status='unknown')
      endif

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
         call getcontrib(ncontrib,icontrib,amaindir(n),asubdir,nsubdir,&
              1,staindx(n))
         
         ! Add the current contribution to the cnorm array
         do k=1,nsubdir
            if (icontrib(k).eq.0) cycle
            itmp=step(k)/dstep+1
            cnorm(ifg,itmp)=cnorm(ifg,itmp)+conjg(coeff(k))*coeff(k)
         enddo

      enddo

      cnorm=sqrt(cnorm)

      close(failunit)
      close(okunit)

      return

    end subroutine getcnorm

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

    subroutine getcontrib(ncontrib,icontrib,amaindir,asubdir,nsubdir,&
         ifailchk,ista)

      use expec, only: lcifilter,cifdthrsh,cifstate

      implicit none

      integer                                :: ncontrib,nsubdir,i,ista
      integer                                :: ifailchk
      integer, dimension(:), allocatable     :: icontrib
      real*8                                 :: dist
      character(len=120)                     :: amaindir
      character(len=120), dimension(nsubdir) :: asubdir
      character(len=250)                     :: filename
      character(len=2)                       :: afail
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
      if (lbound) then
         do i=1,nsubdir
            ! IP-ADC calculation
            filename=trim(amaindir)//'/'//trim(asubdir(i))//'/adc_ip/davstates.dat'
            inquire(file=filename,exist=exists)
            if (.not.exists) then
               afail='ip'
               goto 10
            endif

            ! ADC, Davidson, x
            filename=trim(amaindir)//'/'//trim(asubdir(i))//'/adc_dav_x/osc.dat'
            inquire(file=filename,exist=exists)
            if (.not.exists) then
               afail='x'
               goto 10
            endif
            
            ! ADC, Davidson, y
            filename=trim(amaindir)//'/'//trim(asubdir(i))//'/adc_dav_y/osc.dat'
            inquire(file=filename,exist=exists)
            if (.not.exists) then
               afail='y'
               goto 10
            endif
            
            ! ADC, Davidson, z
            filename=trim(amaindir)//'/'//trim(asubdir(i))//'/adc_dav_z/osc.dat'
            inquire(file=filename,exist=exists)
            if (.not.exists) then
               afail='z'
               goto 10
            endif

            ! CI filtering
            if (lcifilter) then
               call getseamdistance(amaindir,asubdir(i),dist,ifailchk)
               if (dist.lt.cifdthrsh) then
                  if (cifstate.eq.0) then
                     goto 10
                  else
                     if (ista.eq.cifstate) goto 10
                  endif
               endif
            endif

            ncontrib=ncontrib+1
            
            ! Write the successful geometries to file
            if (ifailchk.eq.1) call okgeom(amaindir,asubdir(i))

            cycle
            
10          continue
            icontrib(i)=0

            ! If the calculation didn't complete due to the failure of
            ! the target matching code, then write the geometry to file
            if (ifailchk.eq.1.and.afail.ne.'ip') &
                 call failgeom(amaindir,asubdir(i),afail)
            
         enddo

      else if (lcont) then
         do i=1,nsubdir
            ! IP-ADC calculation
            filename=trim(amaindir)//'/'//trim(asubdir(i))//'/adc_ip/davstates.dat'
            inquire(file=filename,exist=exists)
            if (.not.exists) goto 20
            
            ! ADC, SI, x
            filename=trim(amaindir)//'/'//trim(asubdir(i))//'/adc_si_x/osc.dat'
            inquire(file=filename,exist=exists)
            if (.not.exists) goto 20

            ! ADC, SI, y
            filename=trim(amaindir)//'/'//trim(asubdir(i))//'/adc_si_y/osc.dat'
            inquire(file=filename,exist=exists)
            if (.not.exists) goto 20

            ! ADC, SI, z
            filename=trim(amaindir)//'/'//trim(asubdir(i))//'/adc_si_z/osc.dat'
            inquire(file=filename,exist=exists)
            if (.not.exists) goto 20

            ncontrib=ncontrib+1
            
            cycle
            
20          continue
            icontrib(i)=0
         enddo
         
      else if (ldyson) then
         do i=1,nsubdir
            ! ADC Dyson orbital calculation
            filename=trim(amaindir)//'/'//trim(asubdir(i))&
                 //'/adc_dys/dyson_norms'
            inquire(file=filename,exist=exists)
            if (.not.exists) then
               icontrib(i)=0
            else
               ncontrib=ncontrib+1
            endif
         enddo

      endif

      return

    end subroutine getcontrib

!#######################################################################

    subroutine getseamdistance(amaindir,asubdir,dist,iwrgeom)

      use sysdef
      use expec
      use intcoo

      implicit none

      integer                           :: k1,k2,k3,iadc,i,j,istart,&
                                           iend,iwrgeom,unit
      real*8, dimension(natm*3)         :: xcoo
      real*8                            :: dist
      character(len=2), dimension(natm) :: aatm
      character(len=120)                :: amaindir
      character(len=120)                :: asubdir,string
      character(len=250)                :: filename
      logical                           :: lopen

      iend=len_trim(amaindir)
      istart=0
      do i=3,iend
         if (amaindir(i-2:i).eq.'../') istart=i+1
      enddo
      if (istart.eq.0) istart=1
      k1=index(amaindir(istart:iend),'/')+istart
      k2=len_trim(amaindir)
      if (index(asubdir,'/').eq.0) then
         k3=len_trim(asubdir)
      else
         k3=len_trim(asubdir)-1
      endif

      filename=trim(amaindir)//'/'//trim(asubdir) &
           //'/adc_dav_z'//'/adc_'//amaindir(k1:k2)//'_' &
           //asubdir(1:k3)//'_dav_z.log'

      iadc=323
      open(iadc,file=filename,form='formatted',status='old')

      ! Read the Cartesian coordinates
5     read(iadc,'(a)') string
      if (index(string,'Angstrom').eq.0) goto 5
      do i=1,3
         read(iadc,*)
      enddo
      do i=1,natm
         read(iadc,*) aatm(i),(xcoo(j),j=i*3-2,i*3)
      enddo

      close(iadc)

      dist=x2int(xcoo,1)

      ! Write the Cartesian coordinates to file if we are closer than
      ! threshold to the seam of interest
      if (iwrgeom.eq.1.and.dist.lt.cifdthrsh) then         
         write(cifunit,'(i2)') natm
         write(cifunit,'(a)') trim(filename)
         do i=1,natm
            write(cifunit,'(a2,3(2x,F10.7))') atlbl(i),(xcoo(j),j=i*3-2,i*3)
         enddo
      endif

      return
      
    end subroutine getseamdistance

!#######################################################################

    subroutine failgeom(amaindir,asubdir,atype)

      use sysdef
      use iomod

      implicit none

      integer                           :: k1,k2,k3,iadc,i,j,istart,&
                                           iend
      real*8, dimension(natm*3)         :: xcoo
      character(len=2), dimension(natm) :: aatm
      character(len=2)                  :: atype
      character(len=120)                :: amaindir
      character(len=120)                :: asubdir,string
      character(len=250)                :: filename
      logical(kind=4)                   :: exists,failed

      iend=len_trim(amaindir)
      istart=0
      do i=3,iend
         if (amaindir(i-2:i).eq.'../') istart=i+1
      enddo
      if (istart.eq.0) istart=1
      k1=index(amaindir(istart:iend),'/')+istart
      k2=len_trim(amaindir)
      if (index(asubdir,'/').eq.0) then
         k3=len_trim(asubdir)
      else
         k3=len_trim(asubdir)-1
      endif
           
      filename=trim(amaindir)//'/'//trim(asubdir) &
           //'/adc_dav_'//trim(atype)//'/adc_'//amaindir(k1:k2)//'_' &
           //asubdir(1:k3)//'_dav_'//trim(atype)//'.log'

      inquire(file=filename,exist=exists)
      
      if (exists) then
         iadc=323
         open(iadc,file=filename,form='formatted',status='old')
         
         ! Read the Cartesian coordinates
5        read(iadc,'(a)',end=999) string
         if (index(string,'Angstrom').eq.0) goto 5
         do i=1,3
            read(iadc,*)
         enddo
         do i=1,natm
            read(iadc,*) aatm(i),(xcoo(j),j=i*3-2,i*3)
         enddo

         ! Check whether the calculation failed due to the target
         ! matching going wrong, and if so write the geometry to file
         failed=.true.
10       read(iadc,'(a)',end=888) string
         if (index(string,'State').ne.0.and.index(string,'(selected)').ne.0) then
            failed=.false.
         else
            goto 10
         endif
888      continue
         if (failed) then
            write(failunit,'(i2)') natm
            write(failunit,'(a)') trim(filename)
            do i=1,natm
               write(failunit,'(a2,3(2x,F10.7))') aatm(i),(xcoo(j),j=i*3-2,i*3)
            enddo
         endif

999      continue
         close(iadc)

      endif

      return

    end subroutine failgeom

!#######################################################################

    subroutine okgeom(amaindir,asubdir)

      use sysdef
      use iomod

      implicit none
      
      integer                           :: k1,k2,k3,iadc,i,j,istart,&
                                           iend
      real*8, dimension(natm*3)         :: xcoo
      character(len=2), dimension(natm) :: aatm
      character(len=120)                :: amaindir
      character(len=120)                :: asubdir,string
      character(len=250)                :: filename

      iend=len_trim(amaindir)
      istart=0
      do i=3,iend
         if (amaindir(i-2:i).eq.'../') istart=i+1
      enddo
      if (istart.eq.0) istart=1
      k1=index(amaindir(istart:iend),'/')+istart
      k2=len_trim(amaindir)
      if (index(asubdir,'/').eq.0) then
         k3=len_trim(asubdir)
      else
         k3=len_trim(asubdir)-1
      endif      

      filename=trim(amaindir)//'/'//trim(asubdir) &
           //'/adc_dav_x'//'/adc_'//amaindir(k1:k2)//'_' &
           //asubdir(1:k3)//'_dav_x'//'.log'
      
      iadc=323
      open(iadc,file=filename,form='formatted',status='old')
      
5     read(iadc,'(a)') string
      if (index(string,'Angstrom').eq.0) goto 5
      do i=1,3
         read(iadc,*)
      enddo
      do i=1,natm
         read(iadc,*) aatm(i),(xcoo(j),j=i*3-2,i*3)
      enddo

      write(okunit,'(i2)') natm
      write(okunit,'(a)') trim(filename)
      do i=1,natm
         write(okunit,'(a2,3(2x,F10.7))') aatm(i),(xcoo(j),j=i*3-2,i*3)
      enddo

      close(iadc)
      
      return

    end subroutine okgeom

!#######################################################################

    subroutine getip_adc(ip,amaindir,asubdir,nsubdir,icontrib)

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
! Allocate the IP and E0 array
!-----------------------------------------------------------------------
      if (allocated(ip)) deallocate(ip)
      allocate(ip(nsubdir))
      ip=0.0d0

!-----------------------------------------------------------------------
! Loop over timesteps/subdirectories and read the IP for each
!-----------------------------------------------------------------------
      if (.not.ldyson) then

         unit=20
         do i=1,nsubdir
            
            if (icontrib(i).eq.0) cycle

            filename=trim(amaindir)//'/'//trim(asubdir(i)) &
                 //'/adc_ip/davstates.dat'
            open(unit,file=filename,form='formatted',status='old')

5           read(unit,'(a)') string
            if (index(string,'MP2 energy:').eq.0) goto 5
            read(string,'(28x,F14.8)') e0

10          read(unit,'(a)') string
            if (index(string,'Energy:').eq.0) goto 10
            read(string,'(27x,F14.8)') ip(i)
            ip(i)=ip(i)-e0

            close(unit)
         enddo
      
      endif

      return

    end subroutine getip_adc

!#######################################################################

    subroutine getxsec(tdmsq,amaindir,asubdir,nsubdir,icontrib,&
         staindx,einit,nstates_f,deltae,ip,coeff,ifg,step)

      implicit none

      integer                                :: nsubdir,staindx,&
                                                nstates_f,ifg
      integer, dimension(nsubdir)            :: icontrib,step
      real*8, dimension(:), allocatable      :: einit
      real*8, dimension(:,:), allocatable    :: tdmsq,deltae
      real*8, dimension(nsubdir)             :: ip
      complex*16, dimension(nsubdir)         :: coeff
      character(len=120)                     :: amaindir
      character(len=120), dimension(nsubdir) :: asubdir
      
      if (lbound) then
         call rddavidson(tdmsq,amaindir,asubdir,nsubdir,icontrib,&
              staindx,einit,nstates_f,deltae)
      else if (lcont) then
         call rdstieltjes(amaindir,asubdir,nsubdir,icontrib,ip,&
              coeff,ifg,step,einit)

      else if (ldyson) then
         call rddysnorm(amaindir,asubdir,nsubdir,icontrib,coeff,step)
      endif

      return

    end subroutine getxsec

!#######################################################################

    subroutine rddavidson(tdmsq,amaindir,asubdir,nsubdir,icontrib,&
         staindx,einit,nstates_f,deltae)
      
      implicit none

      integer                                :: nsubdir,staindx,unit,&
                                                i,j,c,count,k1,k2,k3,&
                                                adcstate,nstates_f,&
                                                indx,istart,iend
      integer, dimension(nsubdir)            :: icontrib
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
      call get_einit(nsubdir,icontrib,amaindir,asubdir,einit)

!-----------------------------------------------------------------------
! Read in the TDM^2 values
!-----------------------------------------------------------------------
      ! Loop over timesteps/subdirectories
      do i=1,nsubdir

         if (icontrib(i).eq.0) cycle

         ! (1) Excitation energies
         iend=len_trim(amaindir)
         istart=0
         do j=3,iend
            if (amaindir(j-2:j).eq.'../') istart=j+1
         enddo
         if (istart.eq.0) istart=1
         k1=index(amaindir(istart:iend),'/')+istart
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

            call get_state(indx,ener,deltae(i,:),nstates_f,filename)

            tdmsq(i,indx)=tdmsq(i,indx)+(3.0d0/(2.0d0*deltae(i,indx)))*osc

            goto 30

35          continue

            close(unit)
            
         enddo

      enddo

      return

    end subroutine rddavidson

!#######################################################################

    subroutine get_einit(nsubdir,icontrib,amaindir,asubdir,einit)

      implicit none

      integer                                :: nsubdir,i,j,k1,k2,k3,unit,&
                                                adcstate,count,istart,iend
      integer, dimension(nsubdir)            :: icontrib
      real*8, dimension(nsubdir)             :: einit
      character(len=120)                     :: amaindir,string
      character(len=120), dimension(nsubdir) :: asubdir
      character(len=250)                     :: filename

      unit=30

      ! Loop over timesteps/subdirectories
      do i=1,nsubdir
         
         if (icontrib(i).eq.0) cycle

         ! (1) Read the index of the ADC sate chosen by the target
         !     matching routine
         iend=len_trim(amaindir)
         istart=0
         do j=3,iend
            if (amaindir(j-2:j).eq.'../') istart=j+1
         enddo
         if (istart.eq.0) istart=1
         k1=index(amaindir(istart:iend),'/')+istart
         k2=len_trim(amaindir)
         if (index(asubdir(i),'/').eq.0) then
            k3=len_trim(asubdir(i))
         else
            k3=len_trim(asubdir(i))-1
         endif

         if (lbound) then
            filename=trim(amaindir)//'/'//trim(asubdir(i)) &
                 //'/adc_dav_x/adc_'//amaindir(k1:k2)//'_' &
                 //asubdir(i)(1:k3)//'_dav_x.log'
         else if (lcont) then
            filename=trim(amaindir)//'/'//trim(asubdir(i)) &
                 //'/adc_si_x/adc_'//amaindir(k1:k2)//'_' &
                 //asubdir(i)(1:k3)//'_si_x.log'
         endif

         open(unit,file=filename,form='formatted',status='old')

5        read(unit,'(a)') string
         if (index(string,'selected').eq.0) goto 5
         
         read(string,'(7x,i2)') adcstate
         
         close(unit)

         ! (2) Read the energy of the initial state
         if (lbound) then
            filename=trim(amaindir)//'/'//trim(asubdir(i)) &
                 //'/adc_dav_x/davstates.dat'
         else if (lcont) then
            filename=trim(amaindir)//'/'//trim(asubdir(i)) &
                 //'/adc_si_x/davstates.dat'
         endif

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

      return

    end subroutine get_einit

!#######################################################################

    subroutine get_state(indx,ener,deltae,nstates_f,filename)
      
      implicit none

      integer                      :: indx,nstates_f,i
      real*8                       :: ener,ftmp
      real*8, dimension(nstates_f) :: deltae
      real*8, parameter            :: tol=0.0001d0
      character(len=250)           :: filename

      indx=-1

      do i=1,nstates_f
         if (abs(ener/eh2ev-deltae(i)).le.tol) then
            indx=i
            exit
         endif
      enddo

      if (indx.eq.-1) then
         write(6,'(2(/,2x,a),/)') 'Error in subroutine get_state &
              in trtxas.f90','File: '//trim(filename)
         stop
      endif

      return

    end subroutine get_state

!#######################################################################

    subroutine rdstieltjes(amaindir,asubdir,nsubdir,icontrib,ip,coeff,&
         ifg,step,einit)
      
      use expec
      use sysdef
      use mcsplinemod

      implicit none
      
      integer                                :: nsubdir,unit,i,c,k,j,&
                                                negrid,k1,k2,ifg,m,n
      integer, dimension(nsubdir)            :: icontrib,step
      real*8, dimension(3,siord-1)           :: si_e,si_f
      real*8, dimension(:), allocatable      :: xsec
      real*8, dimension(nsubdir)             :: ip,e0
      real*8                                 :: csq,t,ipactual
      real*8, dimension(:), allocatable      :: einit
      complex*16, dimension(nsubdir)         :: coeff
      character(len=120), dimension(nsubdir) :: asubdir
      character(len=120)                     :: amaindir
      character(len=1), dimension(3)         :: dpllbl
      character(len=250)                     :: filename

      real*8 :: tsig,tcurr

!-----------------------------------------------------------------------
! Initialisation of various things
!-----------------------------------------------------------------------
      unit=20
      dpllbl=(/ 'x','y','z' /)

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      negrid=int(egrid(3))
      allocate(xsec(negrid))

      if (allocated(einit)) deallocate(einit)
      allocate(einit(nsubdir))
      einit=0.0d0

!-----------------------------------------------------------------------
! Read the initial state energies for each timestep/subdirectory
!-----------------------------------------------------------------------
      call get_einit(nsubdir,icontrib,amaindir,asubdir,einit)

!-----------------------------------------------------------------------
! Read the ground state energies for each timestep/subdirectory
!-----------------------------------------------------------------------
      call get_e0(nsubdir,icontrib,amaindir,asubdir,e0)

!-----------------------------------------------------------------------
! Read in the discretised Stieltjes cross-sections
!-----------------------------------------------------------------------
      ! Loop over timesteps/subdirectories
      do i=1,nsubdir
         
         ! Cycle if the current timestep doesn't contribute
         if (icontrib(i).eq.0) cycle

         ! Current time in fs
         tcurr=(step(i)-1)*dt/41.341375d0

         ! Read the Stieltjes imaging cross-sections
         do c=1,3
            filename=trim(amaindir)//'/'//trim(asubdir(i)) &
                 //'/adc_si_'//dpllbl(c)//'/xsec_order'            
            k=len_trim(filename)
            if (siord.lt.10) then
               write(filename(k+1:k+3),'(a2,i1)') '00',siord
            else 
               write(filename(k+1:k+3),'(a1,i2)') '0',siord
            endif
            open(unit,file=filename,form='formatted',status='old')
            do j=1,siord-1
               read(unit,*) si_e(c,j),si_f(c,j)
            enddo
            close(unit)
         enddo

         ! Interpolate to determine the cross-sections at the energy 
         ! grid points
         xsec=0.0d0
         ipactual=ip(i)-(einit(i)-e0(i))
         call interpolate_stieltjes(si_e,si_f,xsec,negrid,siord,&
              egrid,ipactual)

         xsec=xsec*(3.0d0/2.0d0)*4.0d0*pi/c_au

         ! Fill in the spectrum
         nfunc=nfunc+1         
         k1=(nfunc-1)*negrid+1
         k2=nfunc*negrid

         csq=conjg(coeff(i))*coeff(i)
         
         ! (1) Amplitude
!         par(k1:k2,1)=csq*xsec(1:negrid)/cnorm(ifg,i)

         ! (2) Centre wrt time
!         par(k1:k2,2)=t

         ! TEST
         tsig=fwhm_t/2.35482d0
         do m=1,int(egrid(3))
            do n=1,int(tgrid(3))
               t=tgrid(1)+(n-1)*((tgrid(2)-tgrid(1))/tgrid(3))
               spec(m,n)=spec(m,n)+(1.0d0/real(nifg)) &
                    * csq*xsec(m)/cnorm(ifg,i) &
                    * exp(-((t-tcurr)**2)/(2.0d0*tsig**2))
            enddo
         enddo
         ! TEST

      enddo

      return

    end subroutine rdstieltjes

!#######################################################################

    subroutine get_e0(nsubdir,icontrib,amaindir,asubdir,e0)

      implicit none

      integer                                :: nsubdir,i,j,k1,k2,k3,unit,&
                                                adcstate,count,istart,iend
      integer, dimension(nsubdir)            :: icontrib
      real*8, dimension(nsubdir)             :: e0
      character(len=120)                     :: amaindir,string
      character(len=120), dimension(nsubdir) :: asubdir
      character(len=250)                     :: filename

      unit=30

      ! Loop over timesteps/subdirectories
      do i=1,nsubdir
         
         if (icontrib(i).eq.0) cycle

         ! (1) Read the index of the ADC sate chosen by the target
         !     matching routine
         iend=len_trim(amaindir)
         istart=0
         do j=3,iend
            if (amaindir(j-2:j).eq.'../') istart=j+1
         enddo
         if (istart.eq.0) istart=1
         k1=index(amaindir(istart:iend),'/')+istart
         k2=len_trim(amaindir)
         if (index(asubdir(i),'/').eq.0) then
            k3=len_trim(asubdir(i))
         else
            k3=len_trim(asubdir(i))-1
         endif

         if (lbound) then
            filename=trim(amaindir)//'/'//trim(asubdir(i)) &
                 //'/adc_dav_x/adc_'//amaindir(k1:k2)//'_' &
                 //asubdir(i)(1:k3)//'_dav_x.log'
         else if (lcont) then
            filename=trim(amaindir)//'/'//trim(asubdir(i)) &
                 //'/adc_si_x/adc_'//amaindir(k1:k2)//'_' &
                 //asubdir(i)(1:k3)//'_si_x.log'
         endif

         open(unit,file=filename,form='formatted',status='old')

5        read(unit,'(a)') string
         if (index(string,'MP2 energy:').eq.0) goto 5
         
         read(string,'(32x,F15.10)') e0(i)
         
         close(unit)

      enddo

      return

    end subroutine get_e0

!#######################################################################

    subroutine rddysnorm(amaindir,asubdir,nsubdir,icontrib,coeff,step)

      use sysdef
      use expec

      implicit none

      integer                                :: nsubdir,i,unit
      integer, dimension(nsubdir)            :: icontrib,step
      real*8                                 :: e,de,norm,t,csq
      complex*16, dimension(nsubdir)         :: coeff
      character(len=120)                     :: amaindir
      character(len=120), dimension(nsubdir) :: asubdir
      character(len=250)                     :: filename

!-----------------------------------------------------------------------
! Read in the Dyson norms and vertical ionization energies
!-----------------------------------------------------------------------
      unit=20

      ! Loop over timesteps/subdirectories
      do i=1,nsubdir

         ! Cycle if the current timestep doesn't contribute
         if (icontrib(i).eq.0) cycle

         ! Current time in fs
         t=(step(i)-1)*dt/41.341375d0

         ! Read the Dyson norm file
         filename=trim(amaindir)//'/'//trim(asubdir(i)) &
              //'/adc_dys/dyson_norms'
         open(unit,file=filename,form='formatted',status='old')
         read(unit,*)
5        read(unit,'(3(2x,F14.8))',end=10) e,de,norm

         if (de.gt.0.0d0.and.de.le.eprobe) then
            ! Fill in the next row in the parameter array
            nfunc=nfunc+1
            
            ! (1) Coefficient
            csq=conjg(coeff(i))*coeff(i)
            par(nfunc,1)=norm*csq

            ! (2) Width parameter, energy domain
            par(nfunc,2)=fwhm_e/2.35482d0
            
            ! (3) Width parameter, time domain
            par(nfunc,3)=fwhm_t/2.35482d0

            ! (4) Centre wrt E
            par(nfunc,4)=eprobe-de
            
            ! (5) Centre wrt t
            par(nfunc,5)=t
         endif

         goto 5

10       continue

         close(unit)

      enddo

      return

    end subroutine rddysnorm

!#######################################################################

    subroutine trtxas_currtraj(nsubdir,icontrib,nstates_f,ip,deltae,&
         tdmsq,coeff,step,ifg)

      use expec
      use sysdef

      implicit none

      integer                              :: nsubdir,nstates_f,n,k,&
                                              ifg
      integer, dimension(nsubdir)          :: icontrib,step
      real*8, dimension(nsubdir)           :: ip
      real*8, dimension(nsubdir,nstates_f) :: deltae,tdmsq
      real*8                               :: tcurr,t,e,dele,delt,csq,&
                                              lfunc,prefac
      complex*16, dimension(nsubdir)       :: coeff

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

    subroutine trtxas_total      

      implicit none

      write(6,'(/,2x,a,/)') 'Constructing the TR-TXAS...'

      if (lbound) then
         call trtxas_total_bound
      else if (lcont) then
         call trtxas_total_cont
      else if (ldyson) then         
         call trpes_adc_total
      endif

      return

    end subroutine trtxas_total

!#######################################################################

    subroutine trtxas_total_bound

      use expec
      use sysdef

      implicit none
      
      integer :: iout,i,j,k
      real*8  :: dele,delt,e,t,func,prefac,musq,shape,deif,tcent,&
                 tsig

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
      do i=1,int(egrid(3))
         write(iout,*)
         do j=1,int(tgrid(3))

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

               ! FWHM in the time domain
               tsig=fwhm_t/2.35482d0

               ! Function value
               func=func+(1.0d0/real(nifg))*prefac*e*musq*shape &
                    * exp(-((t-tcent)**2)/(2.0d0*tsig**2))

            enddo

            ! Ouput the current photon energy, time delay and
            ! function value
            write(iout,*) e*eh2ev,t,func

         enddo
      enddo

      return

    end subroutine trtxas_total_bound

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

    subroutine trtxas_total_cont

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
      open(iout,file='trtxas.dat',form='formatted',status='unknown')

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

    end subroutine trtxas_total_cont

!#######################################################################

    subroutine trpes_adc_total

      use expec
      use sysdef

      implicit none

      integer :: i,j,k,iout
      real*8  :: a,esig,tsig,tcent,ecent,e,t,dele,delt,func

!-----------------------------------------------------------------------
! Convert the energy grid back to units of eV
!-----------------------------------------------------------------------
      egrid(1:2)=egrid(1:2)*eh2ev
      
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
      do i=1,int(egrid(3))+1
         write(iout,*)
         do j=1,int(tgrid(3))+1

            e=egrid(1)+(i-1)*dele
            t=tgrid(1)+(j-1)*delt

            ! Loop over Gaussians, summing the contributions from each
            func=0.0d0
            do k=1,nfunc
               ! Set the current parameters
               a=par(k,1)
               esig=par(k,2)
               tsig=par(k,3)
               ecent=par(k,4)
               tcent=par(k,5)
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

    end subroutine trpes_adc_total

!#######################################################################

  end module adctas

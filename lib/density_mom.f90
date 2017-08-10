  module density_mom
    
    implicit none

  contains

!#######################################################################
! reddens_mom: calculates the reduced density in the momentum
!              representation wrt a given internal coordinate using a 
!              Monte Carlo procedure as described in
!              J. Phys. Chem. A, 111, 11305 (2007)
!
! As it stands, this routine is very confusing as we end up with the
! momenta being held in the arrays traj%r and the negatives of the
! positions held in the arrays traj%p
!
! Hence, on the surface things may look erroneous but, hopefully, they
! are actually OK...
!#######################################################################

    subroutine reddens_mom

      use trajdef
      use sysdef
      use intcoo
      use expec
      use gausstools

      implicit none

      integer                          :: n,m,i,ibas,ibin,itmp,nalive,&
                                          iout,count
      real*8, dimension(3*natm)        :: xcoo,pcoo
      real*8                           :: icoo,dens,impfunc
      real*8, dimension(int(dgrid(4))) :: cent
      logical(kind=4)                  :: lpop,lbound

!-----------------------------------------------------------------------
! Exit if the permutation of identical nuclei has been requested - I
! currently cannot be arsed to implement this
!-----------------------------------------------------------------------
      if (npermute.gt.0) then
         print*,
         print*,"THE PERMUTATION OF IDENTICAL NUCLEI IS NOT YET &
              SUPPORTED IN THE MOMENTUM REDUCED DENSITY CODE"
         print*,
         STOP
      endif
      
!-----------------------------------------------------------------------
! Transform the FMS wavefunction to the momentum representation
!-----------------------------------------------------------------------
      call pos2mom

!-----------------------------------------------------------------------
! Initialise arrays
!-----------------------------------------------------------------------
      call densinit

!-----------------------------------------------------------------------
! Open output file
!-----------------------------------------------------------------------
      iout=30
      open(iout,file='dens.dat',form='formatted',status='unknown')

!-----------------------------------------------------------------------
! Get the centres of the partitions of the internal coordinate value
! interval
!-----------------------------------------------------------------------
      call getcent(cent)

!-----------------------------------------------------------------------
! If we are dealing with a rectilinear Cartesian vector, then read the
! vector file
!-----------------------------------------------------------------------
      if (ityp.lt.0) call rdvec

!-----------------------------------------------------------------------
! Calculate the reduced density integrated over each partition of the
! internal coordinate interval via Monte Carlo integration
!-----------------------------------------------------------------------
      write(6,'(/,2x,a,/)') 'Calculating reduced densities...'
      
      ! Loop over timesteps
      do n=1,nstep,dstep

         write(6,'(2x,a3,F12.4,1x,a4)') 't =',dt*(n-1),'a.u.'

         pfunc=0.0d0

         ! Loop over IFGs
         do i=1,nintraj

            ! Loop over Monte Carlo steps
            do m=1,nmc

               ! Make sure that we have no contributions from 
               ! IFGs that are dead
               lpop=ispop(i,n)
               if (.not.lpop) cycle

               ! Select a live trajectory (indexed by ibas)
               call select_traj(ibas,i,n)
               
               ! Using the selected trajectory, sample the Cartesian
               ! momenta using the corresponding Gaussian
               ! distribution
               count=0
10             continue
               count=count+1
               call sample_cart(ibas,i,n,pcoo,xcoo)
              
               ! Calculate the internal coordinate of interest at the
               ! chosen current geometry
               icoo=x2int_mom(pcoo,xcoo)
               
               ! If the internal coordinate value is not contained within
               ! the user specified interval, then sample a different
               ! Cartesian geometry
               lbound=isbound(icoo)
               if (.not.lbound) then
                  if (count.gt.8000) then
                     goto 999
                  else
                     goto 10
                  endif
               endif

               if (.not.lbound) goto 10
               
               ! Determine the bin that contains the internal coordinate 
               ! value 
               ibin=getbin(icoo)
               
               ! Calculate |Psi|^2 at the current geometry
               call psixpsi(dens,pcoo,i,n)

               ! Calculate the value of the importance sampling function at
               ! the current geometry
               call importance(impfunc,pcoo,i,n)

               ! Get the no. IFGs alive at the current timestep
               call naliveifg(nalive,n)

               ! Add the importance-function-weighted density to the current
               ! bin
               pfunc(ibin)=pfunc(ibin)+dens/(impfunc*nmc*nalive)

            enddo
            
         enddo
         
!         ! Check that the probabilities sum to unity
!         sum=0.0d0
!         do k=1,int(dgrid(4))
!            sum=sum+pfunc(k)
!         enddo
!         print*,"SUM:",sum

         ! Output the integrated densities at the current timestep
         call outintdens(iout,n,pfunc,cent)

      enddo
      
!-----------------------------------------------------------------------
! Close the output file
!-----------------------------------------------------------------------
      close(iout)

      return
      
999   continue
      write (6,'(/,2x,a,/)') 'It has not been possible to sample a &
           geometry within the given coordinate range. Please increase &
           this range'
      STOP

    end subroutine reddens_mom

!#######################################################################
! pos2mom: constructs the FMS wavefunction in the momentum
!#######################################################################

    subroutine pos2mom

      use sysdef
      use trajdef

      implicit none

      integer                               :: i,ntraj,istep,itraj
      real*8, dimension(natm*3)             :: alpha_tmp
      real*8, dimension(:,:,:), allocatable :: r_tmp,p_tmp

!-----------------------------------------------------------------------
! Gaussian widths
!-----------------------------------------------------------------------
      do i=1,natm*3
         alpha_tmp(i)=1.0d0/(4.0d0*alpha(i))
      enddo
      
      alpha=alpha_tmp

!-----------------------------------------------------------------------
! Gaussian centres in phase space
!
! Confusingly, we load the momenta into the position array and vice
! versa so that we can use the pre-existing Gaussian overlap routines
! as they are
!
! P0 -> -R0
! R0 -> +P0
!-----------------------------------------------------------------------
      do i=1,nintraj
         ntraj=traj(i)%ntraj
         allocate(r_tmp(ntraj,nstep,natm*3),p_tmp(ntraj,nstep,natm*3))
         do itraj=1,ntraj
            do istep=1,nstep
               r_tmp(itraj,istep,:)=traj(i)%p(itraj,istep,:)
               p_tmp(itraj,istep,:)=-traj(i)%r(itraj,istep,:)               
            enddo
         enddo
         traj(i)%p=p_tmp
         traj(i)%r=r_tmp
         deallocate(r_tmp,p_tmp)
      enddo

      return

    end subroutine pos2mom

!#######################################################################
! densinit: allocates and initialises allocatable arrays used in the
!           calculation of the reduced density
!#######################################################################

    subroutine densinit

      use expec, only: pfunc,dgrid,dstep
      use sysdef, only: nstep

      implicit none
      
      integer :: itmp

!-----------------------------------------------------------------------
! pfunc: probabilities of the internal coordinate being in each
!        partition of the coordinate interval at each timestep
!-----------------------------------------------------------------------
      itmp=int(dgrid(4))
      allocate(pfunc(itmp))

      return

    end subroutine densinit

!#######################################################################
! getcent: determines the centres of the partitions of the internal
!          coordinate value interval
!#######################################################################

    subroutine getcent(cent)

      use expec, only: dgrid

      implicit none

      integer                          :: npart,i
      real*8, dimension(int(dgrid(4))) :: cent
      real*8                           :: a,b

      npart=int(dgrid(4))

      do i=1,npart
         a=dgrid(1)+(i-1)*dgrid(3)
         b=dgrid(1)+i*dgrid(3)
         cent(i)=a+0.5*(b-a)
      enddo

      return

    end subroutine getcent

!#######################################################################
! rdvec: reads the Cartesian vector from file
!#######################################################################

    subroutine rdvec

      use sysdef
      use expec
      use parsemod

      implicit none

      integer                              :: unit,i,k
      integer                              :: inkw
      integer, parameter                   :: maxkw=200
      integer, dimension(maxkw)            :: ilkw
      real*8                               :: dp
      character(len=120), dimension(maxkw) :: keyword

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(cartvec(natm*3))
      allocate(refgeom(natm*3))
      cartvec=0.0d0
      refgeom=0.0d0

!-----------------------------------------------------------------------
! Open the vector file
!-----------------------------------------------------------------------
      unit=20
      open(unit,file=vecfile,form='formatted',status='old')

!-----------------------------------------------------------------------
! Read the reference geometry and vector
!-----------------------------------------------------------------------
      read(unit,*)
      read(unit,*)
      do i=1,natm
         call rdinp(unit,keyword,inkw,ilkw)
         do k=1,3
            read(keyword(1+k),*) refgeom(i*3-3+k)
            read(keyword(4+k),*) cartvec(i*3-3+k)
         enddo         
      enddo

!-----------------------------------------------------------------------
! Close the vector file
!-----------------------------------------------------------------------
      close(unit)

!-----------------------------------------------------------------------
! Normalise the Cartesian vector
!-----------------------------------------------------------------------
      dp=dot_product(cartvec,cartvec)
      cartvec=cartvec/sqrt(dp)

      return

    end subroutine rdvec

!#######################################################################
! select_traj: selects a trajectory index j according to the 
!              distribution f(j)=|C_j|^2
!#######################################################################

    subroutine select_traj(ibas,itraj,istep)

      use trajdef

      implicit none

      integer    :: ibas,itraj,istep,j,ntraj
      real*8     :: randm,pop
      complex*16 :: coe

      ntraj=traj(itraj)%ntraj

      ibas=0
10    continue

      ! Randomly select a trajectory
      call random_number(randm)
      j=ceiling(ntraj*randm)

      ! Choose the selected trajectory if
      ! random no. .le. population of the trajectory
      call random_number(randm)
      coe=traj(itraj)%coe(j,istep)
      pop=real(conjg(coe)*coe)
      if (randm.le.real(pop)) then
         ibas=j
         goto 20
      endif
      
      if (ibas.eq.0) goto 10

20    continue

      return

    end subroutine select_traj

!#######################################################################
! sample_cart: For a given trajectory at a given timestep, samples the
!              cartesian coordinates xcoo from the corresponding
!              Gaussian distributions
!#######################################################################

    subroutine sample_cart(ibas,itraj,istep,pcoo,xcoo)

      use trajdef
      use sysdef

      implicit none

      integer                   :: ibas,itraj,istep,i,j
      real*8, dimension(3*natm) :: xcoo,pcoo
      real*8                    :: pcent,sigma,dx1,dx2,rsq

      ! Loop over Cartesian coordinates
      do i=1,3*natm

         ! Centre of the selected trajectory
         pcent=traj(itraj)%r(ibas,istep,i)

         ! Generate random Cartesian coordinates according to
         ! the Gaussian distribution associated with the
         ! selected trajectory
         sigma = sqrt(1./(4.*alpha(i)))
         rsq = 2.
         do
            if(rsq.lt.1.d0.and.rsq.ne.0.d0)exit
            call random_number(dx1)
            call random_number(dx2)
            dx1=2.d0*dx1-1.d0
            dx2=2.d0*dx2-1.d0
            rsq=dx1*dx1+dx2*dx2
         enddo

         pcoo(i)=pcent+sigma*dx1*sqrt(-2.d0*log(rsq)/rsq)

      enddo

      ! Note that we have the mapping R0 <-> -P0
      xcoo=-traj(itraj)%p(ibas,istep,:)

      return

    end subroutine sample_cart

!#######################################################################

    function isbound(icoo)
      
      use expec

      implicit none
      
      real*8          :: icoo,diff1,diff2
      logical(kind=4) :: isbound

      isbound=.true.
      
      diff1=icoo-dgrid(1)
      diff2=dgrid(2)-icoo
      
      if (diff1.lt.0.or.diff2.lt.0) isbound=.false.

      return

    end function isbound

!#######################################################################
! get_bin: determines the number of the partition (bin) into which
!          the current value of the internal coordinate (icoo) falls
!#######################################################################

    function getbin(icoo) result(ibin)

      use expec

      implicit none

      integer :: ibin
      real*8  :: icoo

      ibin=ceiling((icoo-dgrid(1))/dgrid(3))

      return

    end function getbin

!#######################################################################

    subroutine importance(func,xcoo,itraj,istep)

      use sysdef
      use trajdef
      use gausstools, only: gxg1d

      implicit none

      integer                   :: itraj,istep,ntraj,j,m
      real*8                    :: func,gxg
      real*8, dimension(3*natm) :: xcoo
      complex*16                :: cj

      ntraj=traj(itraj)%ntraj

      func=0.0d0
      do j=1,ntraj
         cj=traj(itraj)%coe(j,istep)
         gxg=1.0d0
         do m=1,3*natm
           gxg=gxg*real(gxg1d(m,j,j,itraj,istep,xcoo(m)))
         enddo
         func=func+gxg*conjg(cj)*cj
      enddo

      return

    end subroutine importance

!#######################################################################
! outintdens: writes to file the reduced densities integerated over
!             each partition at the istep'th timestep
!#######################################################################

    subroutine outintdens(iout,istep,pfunc,cent)

      use expec, only: dgrid
      use sysdef

      implicit none

      integer                          :: istep,iout,npart,i
      real*8, dimension(int(dgrid(4))) :: pfunc,cent
      real*8                           :: t

      ! Current time in fs
      t=dt*(istep-1)/41.341375d0

      ! No. partitions
      npart=int(dgrid(4))

      ! Ouput the integrated densities
      write(iout,*)
      do i=1,npart
         ! Temporary bodge
         if (pfunc(i).lt.0.5) write(iout,*) t,cent(i),pfunc(i)
      enddo

      return

    end subroutine outintdens

!#######################################################################

  end module density_mom

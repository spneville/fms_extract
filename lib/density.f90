  module density
    
    implicit none

  contains

!#######################################################################
! reddens: calculates the reduced density wrt a given internal
!          coordinate using a Monte Carlo procedure as described in
!          J. Phys. Chem. A, 111, 11305 (2007)
!#######################################################################

    subroutine reddens

      use trajdef
      use sysdef
      use intcoo
      use expec
      use gausstools

      implicit none

      integer*8                        :: n,m,i,ibas,ibin,itmp,nalive,&
                                          iout
      real*8, dimension(3*natm)        :: xcoo
      real*8                           :: icoo,dens,impfunc
      real*8, dimension(int(dgrid(4))) :: cent
      logical(kind=4)                  :: lpop,lbound

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
               
               ! Using the selected trajectory, sample Cartesian
               ! coordinates using the corresponding Gaussian
               ! distribution
10             continue
               call sample_cart(ibas,i,n,xcoo)
              
               ! Calculate the internal coordinate of interest at the
               ! chosen current geometry
               icoo=x2int(xcoo)
               
               ! If the internal coordinate value is not contained within
               ! the user specified interval, then sample a different
               ! Cartesian geometry
               lbound=isbound(icoo)
               if (.not.lbound) goto 10
               
               ! Determine the bin that contains the internal coordinate 
               ! value 
               ibin=getbin(icoo)
               
               ! Calculate |Psi|^2 at the current geometry
               call psixpsi(dens,xcoo,i,n)

               ! Calculate the value of the importance sampling function at
               ! the current geometry
               call importance(impfunc,xcoo,i,n)

               ! Get the no. IFGs alive at the current timestep
               call naliveifg(nalive,n)

               ! Add the importance-function-weighted density to the current
               ! bin
               pfunc(ibin)=pfunc(ibin)+dens/(impfunc*nmc*nalive)

            enddo
            
         enddo
         
         ! Check that the probabilities sum to unity
!         sum=0.0d0
!         do k=1,int(dgrid(4))
!            sum=sum+pfunc(n,k)
!         enddo
!         print*,"SUM:",sum,"NALIVE",nalive

         ! Output the integrated densities at the current timestep
         call outintdens(iout,n,pfunc,cent)

      enddo
      
!-----------------------------------------------------------------------
! Close the output file
!-----------------------------------------------------------------------
      close(iout)


      return
  
    end subroutine reddens

!#######################################################################
! densinit: allocates and initialises allocatable arrays used in the
!           calculation of the reduced density
!#######################################################################

    subroutine densinit

      use expec, only: pfunc,dgrid,dstep
      use sysdef, only: nstep

      implicit none
      
      integer*8 :: itmp

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

      integer*8                        :: npart,i
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
! select_traj: selects a trajectory index j according to the 
!              distribution f(j)=|C_j|^2
!#######################################################################

    subroutine select_traj(ibas,itraj,istep)

      use trajdef

      implicit none

      integer*8  :: ibas,itraj,istep,j,ntraj
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

    subroutine sample_cart(ibas,itraj,istep,xcoo)

      use trajdef
      use sysdef

      implicit none

      integer*8                 :: ibas,itraj,istep,i,j
      real*8, dimension(3*natm) :: xcoo
      real*8                    :: rcent,sigma,dx1,dx2,rsq

      ! Loop over Cartesian coordinates
      do i=1,3*natm

         ! Centre of the selected trajectory
         rcent=traj(itraj)%r(ibas,istep,i)

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

         xcoo(i)=rcent+sigma*dx1*sqrt(-2.d0*log(rsq)/rsq)

      enddo

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

      integer*8 :: ibin
      real*8    :: icoo

      ibin=ceiling((icoo-dgrid(1))/dgrid(3))

      return

    end function getbin

!#######################################################################

    subroutine importance(func,xcoo,itraj,istep)

      use sysdef
      use trajdef
      use gausstools, only: gxg1d

      implicit none

      integer*8                 :: itraj,istep,ntraj,j,m
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

      integer*8                        :: istep,iout,npart,i
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
         if (pfunc(i).lt.0.5) write(iout,*),t,cent(i),pfunc(i)
      enddo

      return

    end subroutine outintdens

!#######################################################################

  end module density

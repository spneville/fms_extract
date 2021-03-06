  module density
    
    implicit none

  contains

!#######################################################################
! reddens: calculates the 1D reduced density wrt a given internal
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

      integer                          :: n,m,i,ibas,ibin,itmp,nalive,&
                                          iout,count
      real*8, dimension(3*natm)        :: xcoo,xicoo
      real*8                           :: icoo,dens,impfunc
      real*8, dimension(int(dgrid(4))) :: cent
      logical(kind=4)                  :: lpop,lbound

!-----------------------------------------------------------------------
! Initialise arrays
!-----------------------------------------------------------------------
      call densinit

!-----------------------------------------------------------------------
! If the coordinate type is the distance from a CI seam, then read
! the g-vector, h-vector and CI geometry files and set up the projector
! onto the branching space
!-----------------------------------------------------------------------
      if (ityp.eq.-2) then
         call rdseamfiles
         call calc_projector
      endif
         
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
      write(6,'(/,2x,a,/)') 'Calculating the 1D reduced densities...'
      
      ! Loop over timesteps
      do n=1,nstep,dstep

         ! Exit if we have gone beyond tfinal
         if (dt*(n-1).gt.tfinal) exit

         ! Output our progress
         write(6,'(2x,a3,F12.4,1x,a4)') 't =',dt*(n-1),'a.u.'

         ! Initialise the density for the current timestep
         pfunc=0.0d0

         ! Loop over IFGs
         do i=1,nintraj

            ! Loop over Monte Carlo steps
            do m=1,nmc

               ! Make sure that we have no contributions from 
               ! IFGs that are dead
               if (dstate.eq.0) then
                  lpop=ispop(i,n)
               else
                  lpop=ispop_staproj(i,n,dstate)
               endif                  
               if (.not.lpop) cycle
               
               ! Select a live trajectory (indexed by ibas)
               if (dstate.eq.0) then
                  call select_traj(ibas,i,n)
               else
                  call select_traj_staproj(ibas,i,n,dstate,lpop)
                  if (.not.lpop) cycle
               endif
               
               ! Using the selected trajectory, sample Cartesian
               ! coordinates using the corresponding Gaussian
               ! distribution
               count=0
10             continue
               count=count+1
               call sample_cart(ibas,i,n,xcoo,xicoo)
              
               ! Calculate the internal coordinate of interest at the
               ! chosen current geometry
               !icoo=x2int(xcoo,1)
               icoo=x2int(xicoo,1)
               
               ! If the internal coordinate value is not contained within
               ! the user specified interval, then sample a different
               ! Cartesian geometry
               lbound=isbound(icoo)
               if (.not.lbound) then
                  if (count.gt.1000) then
!                     goto 999
                     cycle
                  else
                     goto 10
                  endif
               endif

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
!               pfunc(ibin)=pfunc(ibin)+dens/(impfunc*nmc*nalive)
               pfunc(ibin)=pfunc(ibin)+dens/(impfunc*nmc*nintraj)
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

999   continue
      write (6,'(/,2x,a,/)') 'It has not been possible to sample a &
           geometry within the given coordinate range. Please increase &
           this range'
      STOP
  
    end subroutine reddens

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
! rdseamfiles: reading of the NACT vector file, the gradient vector
!              files and the CI geometry files    
!#######################################################################

    subroutine rdseamfiles

      use sysdef
      use expec
      
      implicit none

      integer                     :: unit,i,j,k
      real*8                      :: norm,dp
      real*8, dimension(2,natm*3) :: grad
      character(len=2)            :: atmp
      
!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(cigeom(3*natm))
      allocate(branchvec(2,3*natm))

!-----------------------------------------------------------------------
! CI geometry
!-----------------------------------------------------------------------
      unit=335
      open(unit,file=cifile,form='formatted',status='old')
      read(unit,*)
      read(unit,*)
      do i=1,natm
         read(unit,*) atmp,(cigeom(j),j=i*3-2,i*3)
      enddo
      close(unit)

      ! Convert to a.u.
      cigeom=cigeom/0.529177249d0

!-----------------------------------------------------------------------
! NACT vector
!-----------------------------------------------------------------------
      open(unit,file=hfile,form='formatted',status='old')
      do i=1,natm
         read(unit,*) (branchvec(1,j),j=i*3-2,i*3)
      enddo
      close(unit)      

!-----------------------------------------------------------------------
! Gradient difference vector
!-----------------------------------------------------------------------
      do k=1,2
         open(unit,file=gfile(k),form='formatted',status='old')
         do i=1,natm
            read(unit,*) (grad(k,j),j=i*3-2,i*3)
         enddo
         close(unit)
         norm=sqrt(dot_product(grad(k,:),grad(k,:)))
         grad(k,:)=grad(k,:)/norm         
      enddo
      branchvec(2,:)=grad(1,:)-grad(2,:)

!-----------------------------------------------------------------------
! Orthonormalisation
!-----------------------------------------------------------------------
      dp=dot_product(branchvec(2,:),branchvec(1,:))
      branchvec(2,:)=branchvec(2,:)-dp*branchvec(1,:)/dot_product(branchvec(1,:),branchvec(1,:))
      dp=dot_product(branchvec(2,:),branchvec(1,:))
      branchvec(2,:)=branchvec(2,:)-dp*branchvec(1,:)/dot_product(branchvec(1,:),branchvec(1,:))
      
      do k=1,2
         norm=sqrt(dot_product(branchvec(k,:),branchvec(k,:)))
         branchvec(k,:)=branchvec(k,:)/norm
      enddo

      return
      
    end subroutine rdseamfiles

!#######################################################################

    subroutine calc_projector

      use sysdef
      use expec
      
      implicit none

      integer :: i,j,k
      
!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(branchproj(3*natm,3*natm))

!-----------------------------------------------------------------------
! Set up the projector onto the branching space
!-----------------------------------------------------------------------
      branchproj=0.0d0
      do k=1,2
         do i=1,3*natm
            do j=1,3*natm
               branchproj(i,j)=branchproj(i,j)+&
                    branchvec(k,i)*branchvec(k,j)
            enddo
         enddo
      enddo
 
      return
      
    end subroutine calc_projector
      
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
! select_traj_staproj: selects a trajectory index j according to the 
!                      distribution f(j)=|C_j|^2 subject to the
!                      constraint that the trajectory has a state index
!                      ista
!#######################################################################

    subroutine select_traj_staproj(ibas,itraj,istep,ista,lpop)

      use trajdef

      implicit none

      integer         :: ibas,itraj,istep,j,ntraj,ista,k,s
      real*8          :: randm,pop
      complex*16      :: coe
      logical(kind=4) :: lpop

      ntraj=traj(itraj)%ntraj

      ! Check to see whether the current IFG has population in
      ! the ista'th electronic state
      lpop=.false.
      do k=1,ntraj
         s=traj(itraj)%ista(k)
         if (s.eq.ista) then
            lpop=.true.
            exit
         endif
      enddo
      if (.not.lpop) return

      ibas=0
10    continue

      ! Randomly select a trajectory in the electronic state of
      ! interest
      call random_number(randm)
      j=ceiling(ntraj*randm)
      s=traj(itraj)%ista(j)
      if (s.ne.ista) goto 10

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

    end subroutine select_traj_staproj

!#######################################################################
! sample_cart: For a given trajectory at a given timestep, samples the
!              cartesian coordinates xcoo from the corresponding
!              Gaussian distributions
!#######################################################################

    subroutine sample_cart(ibas,itraj,istep,xcoo,xicoo)

      use trajdef
      use sysdef
      use expec
      use geomtools
      
      implicit none

      integer                   :: ibas,itraj,istep,i,j
      real*8, dimension(3*natm) :: xcoo,xicoo
      real*8                    :: rcent,sigma,dx1,dx2,rsq


!-----------------------------------------------------------------------
! Sample the Cartesian coordinates
!-----------------------------------------------------------------------
      ! Loop over Cartesian coordinates
      do i=1,3*natm

         ! Centre of the selected trajectory
         rcent=traj(itraj)%r(ibas,istep,i)
                  
         ! Generate random Cartesian coordinates according to
         ! the Gaussian distribution associated with the
         ! selected trajectory
         sigma=sqrt(1.0d0/(4.0d0*alpha(i)))
         rsq=2.0d0
         do
            if (rsq.lt.1.0d0.and.rsq.ne.0.0d0) exit
            call random_number(dx1)
            call random_number(dx2)
            dx1=2.d0*dx1-1.d0
            dx2=2.d0*dx2-1.d0
            rsq=dx1*dx1+dx2*dx2
         enddo

         xcoo(i)=rcent+sigma*dx1*sqrt(-2.d0*log(rsq)/rsq)

      enddo

!-----------------------------------------------------------------------
! If requested, put the sampled geometry into maximum coincidence with
! the reference geometry subject to the permutation of a set of
! identical nuclei
!
! Note the xicoo is only used in the calculation of the internal
! coordinate value, whilst xcoo is used in the calculation of the
! importance sampling function value
!-----------------------------------------------------------------------
      if (npermute.gt.0) then
         xicoo=maxcoinc(r0,xcoo)
      else
         xicoo=xcoo
      endif
      
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

  end module density

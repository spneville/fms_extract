  module gausstools

    implicit none

    contains

!#######################################################################
! ispop: returns .true. if the population of a bundle of trajectories
!        at a given timestep is above a threshold value (10^-5)
!#######################################################################

    function ispop(itraj,istep)

      use trajdef

      implicit none

      integer           :: itraj,istep,j,ntraj
      real*8, parameter :: tol=10d-5
      complex*16        :: coe,coe2
      logical(kind=4)   :: ispop

      ntraj=traj(itraj)%ntraj

      ispop=.false.
      do j=1,ntraj
         coe=traj(itraj)%coe(j,istep)
         coe2=conjg(coe)*coe
         if (real(coe2).gt.tol) ispop=.true.
      enddo

      return

    end function ispop

!#######################################################################
! ispop_staproj: returns .true. if the population of a bundle of
!                trajectories projected onto a given state at a given
!                timestep is above a threshold value (10^-5)
!#######################################################################
    
    function ispop_staproj(itraj,istep,ista) result(ispop)

      use trajdef

      implicit none

      integer           :: itraj,istep,ista,j,ntraj,s
      real*8, parameter :: tol=10d-5
      complex*16        :: coe,coe2
      logical(kind=4)   :: ispop

      ntraj=traj(itraj)%ntraj

      ispop=.false.
      do j=1,ntraj
         s=traj(itraj)%ista(j)
         if (s.ne.ista) cycle
         coe=traj(itraj)%coe(j,istep)
         coe2=conjg(coe)*coe
         if (real(coe2).gt.tol) ispop=.true.
      enddo

      return

    end function ispop_staproj

!#######################################################################
! ispop_1bas: returns .true. if a given Gaussian basis function is
!             populated at a given timestep
!#######################################################################

    function ispop_1bas(ifg,ibas,istep) result(ispop)

      use trajdef
      
      implicit none

      integer           :: ifg,ibas,istep
      real*8, parameter :: tol=10d-5
      logical           :: ispop
      complex*16        :: coe,coe2

      coe=traj(ifg)%coe(ibas,istep)

      coe2=conjg(coe)*coe

      if (real(coe2).gt.tol) then
         ispop=.true.
      else
         ispop=.false.
      endif
         
      return
      
    end function ispop_1bas
    
!#######################################################################
! psixpsi: calculates the value of |Psi(x)|^2 for the Cartesian
!          coordinates x held in xcoo 
!#######################################################################

    subroutine psixpsi(dens,xcoo,itraj,istep)

      use sysdef
      use trajdef
      
      implicit none

      integer                   :: m,j,k,itraj,istep,ntraj,s1,s2
      real*8, dimension(3*natm) :: xcoo
      real*8                    :: dens,gammaj,gammak
      complex*16                :: im,czero,cj,ck,gjxgk

      im=(0.0d0,1.0d0)
      czero=(0.0d0,0.0d0)

      ntraj=traj(itraj)%ntraj

      dens=0.0d0

      do j=1,ntraj
         do k=1,ntraj

            s1=traj(itraj)%ista(j)
            s2=traj(itraj)%ista(k)
            if (s1.ne.s2) cycle

            cj=traj(itraj)%coe(j,istep)
            ck=traj(itraj)%coe(k,istep)
            if (cj.eq.czero.or.ck.eq.czero) cycle

            gammaj=traj(itraj)%phase(j,istep)
            gammak=traj(itraj)%phase(k,istep)

            gjxgk=exp(im*(gammak-gammaj))

            do m=1,natm*3
               gjxgk=gjxgk*gxg1d(m,j,k,itraj,istep,xcoo(m))
            enddo
            
            dens=dens+real(conjg(cj)*ck*gjxgk)

         enddo
      enddo

      return

    end subroutine psixpsi

!#######################################################################

    function gxg1d(m,j,k,itraj,istep,x)

      use sysdef
      use trajdef

      implicit none

      integer    :: m,j,k,itraj,istep
      real*8     :: x,rj,rk,pj,pk,a,pre,pi,eta,xi
      complex*16 :: gxg1d,im

      pi=4.0d0*datan(1.0d0)
      im=(0.0d0,1.0d0)

!-----------------------------------------------------------------------
! Gaussian parameters
!-----------------------------------------------------------------------
      rj=traj(itraj)%r(j,istep,m)
      rk=traj(itraj)%r(k,istep,m)
      pj=traj(itraj)%p(j,istep,m)
      pk=traj(itraj)%p(k,istep,m)
      a=alpha(m)

!-----------------------------------------------------------------------
! Prefactor
!-----------------------------------------------------------------------
      pre=sqrt(2.0d0*a/pi)

!-----------------------------------------------------------------------
! Real part of the exponent
!-----------------------------------------------------------------------
      eta=a*( (x-rj)**2 + (x-rk)**2 )

!-----------------------------------------------------------------------
! Imaginary part of the exponent
!-----------------------------------------------------------------------
      xi=pk*(x-rk)-pj*(x-rj)

!-----------------------------------------------------------------------
! gj* x gk
!-----------------------------------------------------------------------
      gxg1d=pre*exp(-eta+im*xi)

      return

    end function gxg1d

!#######################################################################

    subroutine naliveifg(nalive,istep)

      use sysdef

      implicit none

      integer         :: nalive,istep,i
      logical(kind=4) :: lalive

      nalive=0
      do i=1,nintraj
         lalive=ispop(i,istep)
         if (lalive) nalive=nalive+1
      enddo

      return

    end subroutine naliveifg

!#######################################################################

    function overlap_general(ifg1,ifg2,traj1,traj2,step1,step2) &
         result(func)

      use sysdef
      use trajdef
      
      implicit none

      integer    :: ifg1,ifg2,traj1,traj2,step1,step2,state1,state2,m
      real*8     :: gamma1,gamma2
      complex*16 :: func,im,czero,c1,c2

      im=(0.0d0,1.0d0)
      czero=(0.0d0,0.0d0)
      
      state1=traj(ifg1)%ista(traj1)
      state2=traj(ifg2)%ista(traj2)

      if (state1.ne.state2) then
         func=czero
      else
         gamma1=traj(ifg1)%phase(traj1,step1)
         gamma2=traj(ifg2)%phase(traj2,step2)
         func=exp(im*(gamma1-gamma2))
         do m=1,natm*3
            func=func*overlap_general1d(m,ifg1,ifg2,traj1,traj2,step1,step2)
         enddo         
      endif
         
      return
      
    end function overlap_general

!#######################################################################

    function overlap_general_nuc_only(ifg1,ifg2,traj1,traj2,step1,step2) &
         result(func)

      use sysdef
      use trajdef
      
      implicit none

      integer    :: ifg1,ifg2,traj1,traj2,step1,step2,m
      real*8     :: gamma1,gamma2
      complex*16 :: func,im,czero,c1,c2

      im=(0.0d0,1.0d0)
      czero=(0.0d0,0.0d0)
      
      gamma1=traj(ifg1)%phase(traj1,step1)
      gamma2=traj(ifg2)%phase(traj2,step2)
      func=exp(im*(gamma1-gamma2))
      do m=1,natm*3
         func=func*overlap_general1d(m,ifg1,ifg2,traj1,traj2,step1,step2)
      enddo

      return
      
    end function overlap_general_nuc_only
    
!#######################################################################

    function overlap_general1d(m,ifg1,ifg2,traj1,traj2,step1,&
         step2) result(func)

      use sysdef
      use trajdef
      
      implicit none

      integer    :: m,ifg1,ifg2,traj1,traj2,step1,step2
      real*8     :: pre,eta,xi,a,r1,r2,p1,p2,rdiff,pdiff,rcent
      complex*16 :: func

!-----------------------------------------------------------------------
! Widths: frozen Gaussian basis, so equal for both gaussians
!-----------------------------------------------------------------------
      a=alpha(m)

!-----------------------------------------------------------------------
! Positions and momenta      
!-----------------------------------------------------------------------
      r1=traj(ifg1)%r(traj1,step1,m)
      p1=traj(ifg1)%p(traj1,step1,m)
      r2=traj(ifg2)%r(traj2,step2,m)
      p2=traj(ifg2)%p(traj2,step2,m)

!-----------------------------------------------------------------------      
! Position and momentum differences, and position centroid
!-----------------------------------------------------------------------      
      rdiff=r1-r2
      pdiff=p1-p2
      rcent=(a*r1+a*r2)/(2.0d0*a)

!-----------------------------------------------------------------------
! Prefactor
!-----------------------------------------------------------------------
      pre=sqrt(2.0d0*sqrt(a*a)/(2.0d0*a))

!-----------------------------------------------------------------------
! Real part of the exponent
!-----------------------------------------------------------------------
      eta=(a*a*rdiff**2+0.25d0*pdiff**2)/(2.0d0*a)

!-----------------------------------------------------------------------
! Imaginary part of the exponent
!-----------------------------------------------------------------------
      xi=(p1*r1-p2*r2)-rcent*pdiff

!-----------------------------------------------------------------------
! One-dimensional overlap
!-----------------------------------------------------------------------
      func=pre*exp(-eta+(0.0d0,1.0d0)*xi)
      
      return
      
    end function overlap_general1d

!#######################################################################
! rcent: for a pair of trajectories indexed n1 and n2 from the IFG
!        indexed ifg, returns the position vector of the correspoding
!        centroid at the istep'th timestep.
!        Here, we make use of the fact that frozen Gaussian basis
!        functions are used, i.e, the widths for the two trajectories
!        are equal.
!#######################################################################
    
    function rcent(n1,n2,ifg,istep)

      use sysdef
      use trajdef

      integer                   :: n1,n2,ifg,istep,i
      real*8, dimension(natm*3) :: rcent
      real*8                    :: a
      
      rcent=0.0d0

      do i=1,natm*3
         a=alpha(i)
         rcent(i)=(a*traj(ifg)%r(n1,istep,i)&
              +a*traj(ifg)%r(n2,istep,i))/(2.0d0*a)
      enddo
      
      return
      
    end function rcent

!#######################################################################
! pcent: for a pair of trajectories indexed n1 and n2 from the IFG
!        indexed ifg, returns the momentum vector of the correspoding
!        centroid at the istep'th timestep.
!        Here, we make use of the fact that frozen Gaussian basis
!        functions are used, i.e, the widths for the two trajectories
!        are equal.
!#######################################################################
    
    function pcent(n1,n2,ifg,istep)

      use sysdef
      use trajdef

      integer                   :: n1,n2,ifg,istep,i
      real*8, dimension(natm*3) :: pcent
      real*8                    :: a
      
      pcent=0.0d0
      
      do i=1,natm*3
         a=alpha(i)
         pcent(i)=(a*traj(ifg)%p(n1,istep,i)&
              +a*traj(ifg)%p(n2,istep,i))/(2.0d0*a)
      enddo
      
      return
      
    end function pcent

!#######################################################################
! nucdip_integral: calculates the nuclear dipole matrix element
!                  < g1 | mu_N | g2 > for a given pair of Gaussian
!                  basis functions at a given timestep
!#######################################################################
    
    function nucdip_integral(ifg,n1,n2,istep) result(func)

      use sysdef
      use trajdef
      
      integer                  :: ifg,n1,n2,istep,i,c,indx
      real*8                   :: a12,z,r1,r2,p1,p2,a1,a2
      complex*16, dimension(3) :: func
      complex*16               :: ovrlp,ci,czero,b12

!-----------------------------------------------------------------------
! Initialisation
!-----------------------------------------------------------------------
      ci=(0.0d0,1.0d0)
      czero=(0.0d0,0.0d0)
      func=czero

!-----------------------------------------------------------------------
! Calculation of < g1 | mu_N | g2 >
!-----------------------------------------------------------------------
      ! Overlap <g1|g2>
      ovrlp=overlap_general_nuc_only(ifg,ifg,n1,n2,istep,istep)

      ! Loop over atoms
      do i=1,natm

         ! Nuclear charge
         z=atnum(i)

         ! Loop over x, y and z components
         do c=1,3

            ! Cartesian coordinate index
            indx=i*3-3+c

            ! Calculation of < g1 | R_indx | g2 >
            r1=traj(ifg)%r(n1,istep,indx)
            r2=traj(ifg)%r(n2,istep,indx)
            p1=traj(ifg)%p(n1,istep,indx)
            p2=traj(ifg)%p(n2,istep,indx)
            a1=alpha(indx)
            a2=alpha(indx)
            a12=2.0d0*alpha(indx)
            b12=-2.0d0*(a1*r1+a2*r2) + ci*(p1-p2)

            func(c)=func(c)+z*(-0.5*b12/a12)

         enddo
         
      enddo

      func=func*ovrlp

      return

    end function nucdip_integral

!#######################################################################

    end module gausstools

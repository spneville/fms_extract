  module gausstools

    implicit none

    contains

!#######################################################################
! ispop: returns .true. if the population of a trajectory at a given
!        timestep is above a threshold value (10^-5)
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
! ispop_staproj: returns .true. if the population of a trajectory
!                projected onto a given state at a given timestep is
!                above a threshold value (10^-5)
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

    end module gausstools

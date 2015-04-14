  module adpop
    
    implicit none

  contains

!#######################################################################
! calcadpop: calculation of adiabatic state populations
!#######################################################################

    subroutine calcadpop
          
      use sysdef
      use trajdef
      use expec

      implicit none
      
      integer*8       :: i,j,k,n,itmp,staj,stak
      complex*16      :: sjk,cj,ck
      logical(kind=4) :: lkill

      write(6,'(/,2x,a,/)') 'Calculating adiabatic populations...'

      ! Loop over timesteps
      do n=1,nstep

         ! Loop over IFGs 
         do i=1,nintraj
         
            ! Loop over all pairs of trajectories for the current IFG
            do j=1,traj(i)%ntraj
               do k=1,traj(i)%ntraj

                  ! Bodges to get around the killing of trajectories:
                  ! (i)  take coefficients at times t>t_kill to be equal 
                  !      to those at t=t_kill;
                  ! (ii) if t>t_kill, then set all off-diagonal elements
                  !      off the overlap matrix for the corresponding 
                  !      trajectory to zero, i.e., assume that the killed
                  !      trajectory has negligible overlap with all other
                  !      trajectories
                  lkill=.false.
                  if (n.gt.traj(i)%tkill(j)) then
                     lkill=.true.
                     itmp=traj(i)%tkill(j)
                     cj=traj(i)%coe(j,itmp)
                  else
                     cj=traj(i)%coe(j,n)
                  endif
                  if (n.gt.traj(i)%tkill(k)) then
                     lkill=.true.
                     itmp=traj(i)%tkill(k)
                     ck=traj(i)%coe(k,itmp)
                  else
                     ck=traj(i)%coe(k,n)
                  endif
                  if (j.ne.k.and.lkill) cycle

                  ! Contributions arise only if both coefficients are
                  ! non-zero
                  if (cj.ne.(0.0d0,0.0d0).and.cj.ne.(0.0d0,0.0d0)) then

                     ! If the trajectory state indices are equal, then
                     ! calculate the contribution to the corresponding
                     ! adiabatic population
                     staj=traj(i)%ista(j)
                     stak=traj(i)%ista(k)
                     if (staj.eq.stak) then
                        sjk=overlap(n,i,j,k)
                        adpop(staj,n)=adpop(staj,n)&
                             +sjk*conjg(cj)*ck/nintraj
                     endif

                  endif
               enddo
            enddo
            
         enddo

      enddo

      return

    end subroutine calcadpop

!#######################################################################
! function overlap: calculates the overlap of two multi-dimensional
!                   Gaussians for the jth and kth trajectories of the
!                   itraj'th IFG
!#######################################################################

    function overlap(istep,itraj,jindx,kindx)

      use sysdef
      use trajdef

      implicit none

      integer*8  :: istep,itraj,jindx,kindx,m      
      real*8     :: phdiff
      complex*16 :: overlap

      phdiff=traj(itraj)%phase(kindx,istep)-traj(itraj)%phase(jindx,istep)
      overlap=exp((0.0d0,1.0d0)*phdiff)      
      do m=1,natm*3
         overlap=overlap*overlap1d(istep,itraj,jindx,kindx,m)         
      enddo
      
      return

    end function overlap

!#######################################################################
! function overlap1d: calculates the overlap of two one-dimensional
!                     Gaussians for the mth dof, jth and kth
!                     trajectories of the itraj'th IFG:
!                     
!                     S_jk,m,itraj = P*exp(eta + i*xi)
!#######################################################################

    function overlap1d(istep,itraj,jindx,kindx,m)

      use sysdef
      use trajdef

      implicit none

      integer*8  :: istep,itraj,jindx,kindx,m,i
      complex*16 :: overlap1d
      real*8     :: pre,eta,xi,a,rj,rk,pj,pk,rdiff,pdiff,rcent

!-----------------------------------------------------------------------
! Widths: frozen Gaussian basis, so equal for both gaussians
!-----------------------------------------------------------------------
      a=alpha(m)

!-----------------------------------------------------------------------      
! Positions and momenta      
!-----------------------------------------------------------------------      
      rj=traj(itraj)%r(jindx,istep,m)
      pj=traj(itraj)%p(jindx,istep,m)
      rk=traj(itraj)%r(kindx,istep,m)
      pk=traj(itraj)%p(kindx,istep,m)

!-----------------------------------------------------------------------      
! Position and momentum differences, and position centroid
!-----------------------------------------------------------------------      
      rdiff=rj-rk
      pdiff=pj-pk
      rcent=(a*rj+a*rk)/(2.0d0*a)

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
      xi = (pj*rj-pk*rk)-rcent*pdiff

!-----------------------------------------------------------------------
! One-dimensional overlap
!-----------------------------------------------------------------------
      overlap1d=pre*exp(-eta+(0.0d0,1.0d0)*xi)

      return

    end function overlap1d

!#######################################################################

  end module adpop

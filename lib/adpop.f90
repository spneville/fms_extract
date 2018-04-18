  module adpop
    
    implicit none

  contains

!#######################################################################
! calcadpop: calculation of adiabatic state populations either with or
!            without a bootstrapping error estimation
!#######################################################################

    subroutine calcadpop
      use expec

      implicit none
      
      if (bootstrap) then
         call get_populations_bootstrap
      else
         call get_populations
      endif
      
      return
      
    end subroutine calcadpop

!#######################################################################
! get_populations: calculation of adiabatic populations
!#######################################################################
    
    subroutine get_populations

      use sysdef
      use trajdef
      use expec

      implicit none
      
      integer         :: i,j,k,n,itmp,staj,stak
      complex*16      :: sjk,cj,ck
      logical(kind=4) :: lkill

!----------------------------------------------------------------------
! Calculate the adiabatic populations
!----------------------------------------------------------------------
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
                  if (cj.ne.(0.0d0,0.0d0).and.ck.ne.(0.0d0,0.0d0)) then

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

!----------------------------------------------------------------------
! Output the adiabatic populations
!----------------------------------------------------------------------
      call wradpop
      
      return
      
    end subroutine get_populations

!#######################################################################
! get_populations_bootstrap: calculation of adiabatic populations with
!                            with a bootstrapping error estimation
!#######################################################################

    subroutine get_populations_bootstrap

      use sysdef
      use trajdef
      use expec

      implicit none
      
      integer                            :: i,i1,j,k,m,n,itmp,staj,stak
      integer, parameter                 :: nsample=2000
      integer, dimension(nsta)           :: unit
      real*8                             :: r
      real*8, allocatable                :: adps(:,:)
      real*8, dimension(nsta)            :: stdev
      complex*16                         :: sjk,cj,ck
      character(len=2)                   :: as
      character(len=18), dimension(nsta) :: aout
      logical(kind=4)                    :: lkill

!----------------------------------------------------------------------
! Calculate the adiabatic populations: this will fill in the adpop
! array
!----------------------------------------------------------------------
      call get_populations

!----------------------------------------------------------------------
! Allocation and initialisation
!----------------------------------------------------------------------
      ! Allocate arrays
      allocate(adps(nsta,nsample))
      adps=0.0d0

      ! Open output files
      do i=1,nsta
         unit(i)=20+i
         write(as,'(i2)') i
         write(aout(i),'(a)') 'adpop_bootstrap.'//trim(adjustl(as))
         open(unit(i),file=aout(i),form='formatted',status='unknown')
      enddo

!----------------------------------------------------------------------
! Bootstrapping error estimation
!----------------------------------------------------------------------
      ! Loop over timesteps
      do n=1,nstep,dstep

         ! Initialise the adps array
         adps=0.0d0
         
         ! Loop over samples
         do m=1,nsample
         
            ! Loop over randomly selected IFGs 
            do i1=1,nintraj
               call random_number(r)
               i=int(r*nintraj)

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
                     if (cj.ne.(0.0d0,0.0d0).and.ck.ne.(0.0d0,0.0d0)) then
                        
                        ! If the trajectory state indices are equal, then
                        ! calculate the contribution to the corresponding
                        ! adiabatic population
                        staj=traj(i)%ista(j)
                        stak=traj(i)%ista(k)
                        if (staj.eq.stak) then
                           sjk=overlap(n,i,j,k)
                           adps(staj,m)=adps(staj,m)&
                                +sjk*conjg(cj)*ck/nintraj
                        endif
                     
                     endif
                  enddo
               enddo
            
            enddo

         enddo

         ! Calculate the standard deviations of the populations in the
         ! adps array
         do i=1,nsta
            stdev(i)=standard_deviation(adps(i,:),nsample)
         enddo

         ! Write the adiabatic populations and standard deviations
         ! at the current timestep to file
         do i=1,nsta
            write(unit(i),'(F10.2,12x,2(2x,F6.4))') &
                 dt*(n-1)/41.341375d0,adpop(i,n),stdev(i)
         enddo
            
      enddo

!----------------------------------------------------------------------
! Deallocate arrays and close files
!----------------------------------------------------------------------
      deallocate(adps)

      do i=1,nsta
         close(unit(i))
      enddo
      
      return
      
    end subroutine get_populations_bootstrap

!#######################################################################
! standard_deviation: calculates the standard deviation of the values
!                     held in the array vals
!#######################################################################

    function standard_deviation(vals,dim) result(sd)

      implicit none

      integer                :: dim,i
      real*8, dimension(dim) :: vals
      real*8                 :: sd,mean

      mean=sum(vals)/dim

      sd=0.0d0
      do i=1,dim
         sd=sd+(vals(i)-mean)**2
      enddo
      sd=sd/dim
      sd=sqrt(sd)
      
      return
      
    end function standard_deviation
    
!#######################################################################
! wradpop: writes the adiabatic populations to file
!#######################################################################
    
    subroutine wradpop

      use sysdef
      use expec

      implicit none
      
      integer          :: i,n,unit
      character(len=8) :: aout

      unit=20

      ! Loop over states
      do i=1,nsta

         ! Write the current filename
         aout=''
         if (i.lt.10) then
            write(aout,'(a6,i1)') 'adpop.',i
         else
            write(aout,'(a6,i2)') 'adpop.',i
         endif

         ! Open current file
         open(unit,file=aout,form='formatted',status='unknown')

         ! Write adiabatic populations for the current state
         do n=1,nstep
            write(unit,'(F10.2,14x,F6.4)') dt*(n-1)/41.341375d0,&
                 adpop(i,n)
         enddo

         ! Close current file
         close(unit)

      enddo

      return

    end subroutine wradpop

    
!#######################################################################
! function overlap: calculates the overlap of two multi-dimensional
!                   Gaussians for the jth and kth trajectories of the
!                   itraj'th IFG
!#######################################################################

    function overlap(istep,itraj,jindx,kindx)

      use sysdef
      use trajdef

      implicit none

      integer    :: istep,itraj,jindx,kindx,m      
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

      integer    :: istep,itraj,jindx,kindx,m,i
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

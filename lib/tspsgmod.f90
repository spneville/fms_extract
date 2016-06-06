  module tspsgmod
    
    implicit none
    
    save
    
    integer                                   :: maxbas,maxbas_trim,basdim,&
                                                 nsample,maxdrop
    integer, dimension(:), allocatable        :: nbas,nbas_trim,ndrop
    integer, dimension(:,:), allocatable      :: indx,colldrop
    complex*16, dimension(:,:,:), allocatable :: smat
    
  contains

!#######################################################################

    subroutine tspsg_prep

      use sysdef
      use expec
      
      implicit none

!-----------------------------------------------------------------------
! Determine the total no. basis functions per electronic state
!-----------------------------------------------------------------------
      call get_nbas_tot

!-----------------------------------------------------------------------
! Prune the basis according to the overlap criterion in order to avoid
! linear dependencies in the sample basis.
!
! For the moment, if Lowdin's canonical orthogonlisation of the basis
! is to be performed, then set ovrthrsh s.t. all basis functions are
! sampled.
!-----------------------------------------------------------------------
      if (bastype.eq.2) ovrthrsh=1.0d0

      call prune_basis

!-----------------------------------------------------------------------
! Lowdin's canonical orthogonalisation of the Gaussian basis
!-----------------------------------------------------------------------
      if (bastype.eq.2) call lowdin_ortho
      
!-----------------------------------------------------------------------
! Write the TS-PSG basis file
!-----------------------------------------------------------------------
      call wrbasfile

!-----------------------------------------------------------------------
! Write the TS-PSG Hamiltonian file
!-----------------------------------------------------------------------
      call wrhamfile

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      deallocate(nbas)

      return

    end subroutine tspsg_prep

!#######################################################################

    subroutine get_nbas_tot

      use trajdef
      use sysdef
      
      implicit none

      integer :: itraj,ntraj,n,ista,istep
      real*8  :: tspawn,tkill

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------      
      allocate(nbas(nsta))
      
!-----------------------------------------------------------------------
! Determine the total number of basis functions per electronic state
!-----------------------------------------------------------------------
      nbas=0
      ! Loop over IFGs
      do itraj=1,nintraj         
         ! Loop over the trajectories for the current IFG
         ntraj=traj(itraj)%ntraj
         do n=1,ntraj
            tspawn=traj(itraj)%tspawn(n)
            tkill=traj(itraj)%tkill(n)
            ista=traj(itraj)%ista(n)
            ! Loop over timesteps
            do istep=1,nstep
               if (istep.le.tspawn) cycle
               if (istep.ge.tkill) cycle
               nbas(ista)=nbas(ista)+1
            enddo
         enddo
      enddo

      ! Maximum total no. basis functions over electronic states
      maxbas=maxval(nbas)

      ! Total no. basis functions across all electronic states
      basdim=sum(nbas(1:nsta))
      
      return
      
    end subroutine get_nbas_tot

!#######################################################################

    subroutine prune_basis

      use trajdef
      use sysdef
      use expec, only: ovrthrsh
      use gausstools
      
      integer           :: itraj,n,ntraj,istep,ista,k,i
      integer           :: ifgnum,trajnum,stepnum
      real*8            :: tspawn,tkill
      complex*16        :: ovr
      logical           :: laccept
      
!-----------------------------------------------------------------------
! Allocate and initialise arrays
!-----------------------------------------------------------------------
      allocate(indx(basdim,3))
      indx=0
      
!-----------------------------------------------------------------------
! Determine the indices of the basis functions (IFG no. and
! trajectory no. within that IFG) that have an acceptably small overlap
! with all other basis functions
!-----------------------------------------------------------------------
      nsample=0
      
      ! Loop over IFGs
      do itraj=1,nintraj

         ! Loop over the trajectories for the current IFG
         ntraj=traj(itraj)%ntraj
         do n=1,ntraj

            tspawn=traj(itraj)%tspawn(n)
            tkill=traj(itraj)%tkill(n)
            ista=traj(itraj)%ista(n)
            
            ! Loop over timesteps
            do istep=1,nstep

               ! Skip if the current basis function is either yet
               ! to spawn or is dead
               if (istep.le.tspawn) cycle
               if (istep.ge.tkill) cycle

               ! Check the overlaps of the current basis function with
               ! the already sampled basis functions
               laccept=.true.
               do k=1,nsample
                  ifgnum=indx(k,1)
                  trajnum=indx(k,2)
                  stepnum=indx(k,3)
                  ovr=overlap_general(itraj,ifgnum,n,trajnum,istep,stepnum)
                  if (abs(ovr).gt.ovrthrsh) then
                     laccept=.false.
                     exit
                  endif
               enddo

               ! Accept the current basis function if it satisfies the
               ! overlap criterion
               if (laccept) then
                  nsample=nsample+1
                  indx(nsample,1)=itraj
                  indx(nsample,2)=n
                  indx(nsample,3)=istep   
               endif
                  
            enddo
               
         enddo
         
      enddo

!-----------------------------------------------------------------------
! Determine the no. of sampled basis functions per electronic state
!-----------------------------------------------------------------------
      allocate(nbas_trim(nsta))

      nbas_trim=0
      do i=1,nsample
         itraj=indx(i,1)
         n=indx(i,2)
         ista=traj(itraj)%ista(n)
         nbas_trim(ista)=nbas_trim(ista)+1
      enddo

      ! Maximum no. of sampled basis functions over electronic states
      maxbas_trim=maxval(nbas_trim)
      
      return
      
    end subroutine prune_basis
      
!#######################################################################

    subroutine wrbasfile

      use trajdef
      use sysdef
      use expec

      implicit none

      integer                               :: unit,itraj,ntraj,n,&
                                               istep,tspawn,tkill,&
                                               ista,i,j
      real*8, dimension(:,:,:), allocatable :: r,p
      real*8, dimension(:,:), allocatable   :: phase
      logical                               :: lowdin
      
!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(r(nsta,maxbas_trim,natm*3))
      allocate(p(nsta,maxbas_trim,natm*3))
      allocate(phase(nsta,maxbas_trim))

!-----------------------------------------------------------------------
! Set up the state-specific basis arrays
!-----------------------------------------------------------------------
      nbas_trim=0
      do i=1,nsample
         itraj=indx(i,1)
         n=indx(i,2)
         istep=indx(i,3)
         ista=traj(itraj)%ista(n)
         nbas_trim(ista)=nbas_trim(ista)+1
         r(ista,nbas_trim(ista),:)=traj(itraj)%r(n,istep,:)
         p(ista,nbas_trim(ista),:)=traj(itraj)%p(n,istep,:)
         phase(ista,nbas_trim(ista))=traj(itraj)%phase(n,istep)         
      enddo

!-----------------------------------------------------------------------
! Open the basis file
!-----------------------------------------------------------------------
      unit=20
      open(unit,file='basis.dat',form='unformatted',status='unknown')

!-----------------------------------------------------------------------
! Write the basis file
!-----------------------------------------------------------------------
      ! Lowdin canonical orthogonalisation flag
      if (bastype.eq.2) then
         lowdin=.true.
      else
         lowdin=.false.
      endif
      write(unit) lowdin
      
      ! No. nuclear coordinates
      write(unit) natm*3

      ! No. electronic states
      write(unit) nsta

      ! Masses
      do i=1,natm
         do j=1,3
            write(unit) atmass(i)
         enddo
      enddo

      ! Frozen Gaussian widths
      write(unit) alpha(1:natm*3)

      ! No. sampled basis functions per electronic state
      do i=1,nsta
         write(unit) nbas_trim(i)
      enddo

      ! Basis function parameters
      do i=1,nsta
      
         ! Positions
         write(unit) r(i,1:nbas_trim(i),1:natm*3)

         ! Momenta
         write(unit) p(i,1:nbas_trim(i),1:natm*3)

         ! Phases
         write(unit) phase(i,1:nbas_trim(i))

      enddo

      ! Lowdin canonical orthogonalisation coefficients
      if (lowdin) then

         ! Numbers of dropped basis functions
         write(unit) ndrop

         ! Expansion coefficients
         do i=1,nsta
            write(unit) smat(i,1:nbas_trim(i),ndrop(i)+1:nbas_trim(i))
         enddo

         ! Indices of the dropped collocation points
         do i=1,nsta
            write(unit) colldrop(i,1:ndrop(i))
         enddo
            
      endif
      
!-----------------------------------------------------------------------
! Close the basis file
!-----------------------------------------------------------------------
      close(unit)

!-----------------------------------------------------------------------
! Ouput the basis information
!
! N.B., This really should also be written to a log file for future
! reference
!-----------------------------------------------------------------------
      write(6,'(/,a)') 'Number of sampled Gaussian basis functions:'
      do i=1,nsta
         write(6,'(a5,x,i2,a1,x,i6)') 'State',i,':',nbas_trim(i)
      enddo

      if (lowdin) then
         write(6,'(/,a)') 'Number of (Lowdin) orthogonalised basis &
              functions:'
         do i=1,nsta
            write(6,'(a5,x,i2,a1,x,i6)') 'State',i,':',nbas_trim(i)-ndrop(i)
         enddo
      endif
         
!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      deallocate(r)
      deallocate(p)
      deallocate(phase)

      return

    end subroutine wrbasfile

!#######################################################################
    
    subroutine wrhamfile

      use trajdef
      use sysdef

      implicit none
      
      integer                                 :: unit,itraj,ntraj,n,&
                                                 tspawn,tkill,ista,&
                                                 istep,i
      real*8, dimension(:,:,:), allocatable   :: ener
      real*8, dimension(:,:,:,:), allocatable :: nact

!-----------------------------------------------------------------------
! Open the Hamiltonian file
!-----------------------------------------------------------------------
      unit=20
      open(unit,file='hamiltonian.dat',form='unformatted',&
           status='unknown')

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(ener(nsta,maxbas_trim,nsta))
      allocate(nact(nsta,maxbas_trim,nsta,natm*3))

!-----------------------------------------------------------------------
! Set up the state-specific basis arrays
!-----------------------------------------------------------------------
      nbas_trim=0
      do i=1,nsample
         itraj=indx(i,1)
         n=indx(i,2)
         istep=indx(i,3)
         ista=traj(itraj)%ista(n)
         nbas_trim(ista)=nbas_trim(ista)+1
         ener(ista,nbas_trim(ista),:)=traj(itraj)%ener(n,istep,:)
         nact(ista,nbas_trim(ista),:,:)=traj(itraj)%nact(n,istep,:,:)
      enddo
      
!-----------------------------------------------------------------------
! Write the Hamiltonian file
!-----------------------------------------------------------------------
      do i=1,nsta

         ! Energies
         write(unit) ener(i,1:nbas_trim(i),:)

         ! NACTs
         write(unit) nact(i,1:nbas_trim(i),:,:)

      enddo

!-----------------------------------------------------------------------
! Close the Hamiltonian file
!-----------------------------------------------------------------------
      close(unit)

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      deallocate(ener)
      deallocate(nact)

      return

    end subroutine wrhamfile

!#######################################################################

    subroutine lowdin_ortho

      use sysdef
      use trajdef
      use gausstools
      use errormod
      
      implicit none

      integer                               :: i,j,k,sta1,sta2,ifg1,&
                                               ifg2,traj1,traj2,step1,&
                                               step2,dim,lwork,info,&
                                               itmp
      integer, dimension(nsta)              :: count1,count2
      integer, dimension(:), allocatable    :: sortindx
      real*8, dimension(:), allocatable     :: rwork,contrib
      real*8, dimension(:,:), allocatable   :: lambda
      real*8, parameter                     :: eigthrsh=1e-5
      complex*16, dimension(:), allocatable :: work
      character(len=120)                    :: errmsg
      
!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(smat(nsta,maxbas_trim,maxbas_trim))
      smat=(0.0d0,0.0d0)

      allocate(lambda(nsta,maxbas_trim))
      lambda=0.0d0

      allocate(ndrop(nsta))
      ndrop=0
      
!-----------------------------------------------------------------------
! Calculate the overlap matrices for each electronic state
!-----------------------------------------------------------------------
      count1=0
      do i=1,nsample

         ifg1=indx(i,1)
         traj1=indx(i,2)
         step1=indx(i,3)
         sta1=traj(ifg1)%ista(traj1)
         count1(sta1)=count1(sta1)+1

         count2=0
         do j=1,nsample
            
            ifg2=indx(j,1)
            traj2=indx(j,2)
            step2=indx(j,3)
            sta2=traj(ifg2)%ista(traj2)
            count2(sta2)=count2(sta2)+1
            
            ! Skip if the two sampled basis functions do not correspond
            ! to the same electronic state
            if (sta1.ne.sta2) cycle

            ! Calculate the current overlap integral
            smat(sta1,count1(sta1),count2(sta1))=&
                 overlap_general(ifg1,ifg2,traj1,traj2,step1,step2)

         enddo
         
      enddo

!-----------------------------------------------------------------------
! Diagonalise the overlap matrix for each electronic state
!-----------------------------------------------------------------------
      do i=1,nsta

         ! Set dimensions and allocate arrays
         dim=nbas_trim(i)
         lwork=max(1,2*dim-1)
         itmp=max(1,3*dim-2)
         allocate(work(lwork))
         allocate(rwork(itmp))

         ! Diagonalise the current overlap matrix
         call zheev('V','U',dim,smat(i,1:dim,1:dim),dim,&
              lambda(i,1:dim),work,lwork,rwork,info)
         if (info.ne.0) then
            errmsg=''
            write(errmsg,'(a,x,i2)') 'Diagonlisation of the overlap &
                 matrix failed for state ',i
            call errcntrl(errmsg)
         endif

         ! Deallocate arrays
         deallocate(work,rwork)
         
      enddo

!-----------------------------------------------------------------------
! Determine the number of basis functions to drop for each electronic
! state
!-----------------------------------------------------------------------
      do i=1,nsta
         do j=1,nbas_trim(i)
            if (lambda(i,j).ge.eigthrsh) exit
            if (lambda(i,j).lt.eigthrsh) ndrop(i)=ndrop(i)+1            
         enddo
      enddo

!-----------------------------------------------------------------------
! Normalise the new set of basis functions
!-----------------------------------------------------------------------
      ! Loop over electronic states
      do i=1,nsta

         ! Loop over the new basis functions for the current
         ! electronic state (omitting the dropped basis functions)
         do j=ndrop(i)+1,nbas_trim(i)
            smat(i,1:nbas_trim(i),j)=smat(i,1:nbas_trim(i),j)/sqrt(lambda(i,j))
         enddo

      enddo

!-----------------------------------------------------------------------
! Determine which collocation points are to be dropped
!
! Currently, we drop the points corresponding to the Gaussians that
! contribute the least to the retained Lowdin orthogonalised basis
! functions
!-----------------------------------------------------------------------
      maxdrop=maxval(ndrop)
      allocate(colldrop(nsta,maxdrop))
      colldrop=0

      ! Loop over electronic states
      do i=1,nsta

         ! Allocate arrays
         allocate(contrib(nbas_trim(i)))
         allocate(sortindx(nbas_trim(i)))
                  
         ! Loop over the sampled Gaussian basis functions
         contrib=0.0d0
         do k=1,nbas_trim(i)
            ! Loop over the RETAINED orthogonalised basis functions
            do j=1,ndrop(i)+1,nbas_trim(i)
               contrib(k)=contrib(k)+real(conjg(smat(i,k,j))*smat(i,k,j))
            enddo
         enddo

         ! Sort the contribution vector
         call dsortindxa1('A',nbas_trim(i),contrib,sortindx)

         ! Save the indices of the collocation points to drop
         do j=1,ndrop(i)
            colldrop(i,j)=sortindx(j)
         enddo

         ! Deallocate arrays
         deallocate(contrib,sortindx)
         
      enddo
      
      return
      
    end subroutine lowdin_ortho

!#######################################################################

    subroutine dsortindxa1(order,ndim,arrin,sortindx)
      
      implicit none

      character(1), intent(in)                :: order
      integer, intent(in)                     :: ndim
      real*8, dimension(ndim), intent(in)     :: arrin
      integer, dimension(ndim), intent(inout) :: sortindx
    
      integer                                 :: i,l,ir,sortindxt,j
      real*8                                  :: q

!!$ The subroutine is taken from the NR p233, employs heapsort.

      do i= 1,ndim
         sortindx(i)=i
      enddo
      
      l=ndim/2+1
      ir=ndim
    
      if (order.eq.'D') then
       
10       continue
         if (l.gt.1) then
            l=l-1
            sortindxt=sortindx(l)
            q=arrin(sortindxt)
         else
            sortindxt=sortindx(ir)
            q=arrin(sortindxt)
            sortindx(ir)=sortindx(1)
            ir=ir-1
            if (ir.eq.1) then
               sortindx(1)=sortindxt
               return
            endif
         endif
       
         i=l
         j=l+l
       
20       if (j.le.ir) then
            if (j.lt.ir) then
               if (arrin(sortindx(j)).gt.arrin(sortindx(j+1))) j=j+1 !
            endif
            if (q.gt.arrin(sortindx(j))) then !
               sortindx(i)=sortindx(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
            goto 20
         endif
         sortindx(i)=sortindxt
         goto 10
       
      else if (order.eq.'A') then
       
100      continue
      if (l.gt.1) then
         l=l-1
         sortindxt=sortindx(l)
         q=arrin(sortindxt)
      else
         sortindxt=sortindx(ir)
         q=arrin(sortindxt)
         sortindx(ir)=sortindx(1)
         ir=ir-1
         if(ir .eq. 1) then
            sortindx(1)=sortindxt
            return
         end if
      end if
      
      i=l
      j=l+l
      
200   if(j .le. ir) then
         if(j .lt. ir) then
            if(arrin(sortindx(j)) .lt. arrin(sortindx(j+1))) j=j+1 !
         end if
         if(q .lt. arrin(sortindx(j))) then !
            sortindx(i)=sortindx(j)
            i=j
            j=j+j
         else
            j=ir+1
         end if
         go to 200
      end if
      sortindx(i)=sortindxt
      go to 100
      
   end if
   
   return
   
  end subroutine dsortindxa1
    
!#######################################################################

  end module tspsgmod


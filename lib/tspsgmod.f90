  module tspsgmod
    
    implicit none
    
    save
    
    integer                              :: maxbas,maxbas2,basdim,nsample
    integer, dimension(:), allocatable   :: nbas,nbas2
    integer, dimension(:,:), allocatable :: indx
    
  contains

!#######################################################################

    subroutine tspsg_prep

      use sysdef

      implicit none

!-----------------------------------------------------------------------
! Determine the total no. basis functions per electronic state
!-----------------------------------------------------------------------
      call get_nbas_tot

!-----------------------------------------------------------------------
! Prune the basis according to the overlap criterion in order to avoid
! linear dependencies in the sample basis
!-----------------------------------------------------------------------
      call prune_basis
      
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
               if (istep.lt.tspawn) cycle
               if (istep.gt.tkill) cycle
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
      use gausstools
      
      integer           :: itraj,n,ntraj,istep,ista,k
      integer           :: ifgnum,trajnum,stepnum
      real*8            :: tspawn,tkill
      real*8, parameter :: ovrthrsh=0.75d0
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
               if (istep.lt.tspawn) cycle
               if (istep.gt.tkill) cycle

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

      return
      
    end subroutine prune_basis
      
!#######################################################################

    subroutine wrbasfile

      use trajdef
      use sysdef

      implicit none

      integer                               :: unit,itraj,ntraj,n,&
                                               istep,tspawn,tkill,&
                                               ista,i,j
      real*8, dimension(:,:,:), allocatable :: r,p
      real*8, dimension(:,:), allocatable   :: phase

!-----------------------------------------------------------------------
! Determine the no. of sampled basis functions per electronic state
!-----------------------------------------------------------------------
      allocate(nbas2(nsta))
      nbas2=0
      do i=1,nsample
         itraj=indx(i,1)
         n=indx(i,2)
         ista=traj(itraj)%ista(n)
         nbas2(ista)=nbas2(ista)+1
      enddo

      ! Maximum no. of sampled basis functions over electronic states
      maxbas2=maxval(nbas2)

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(r(nsta,maxbas2,natm*3))
      allocate(p(nsta,maxbas2,natm*3))
      allocate(phase(nsta,maxbas2))

!-----------------------------------------------------------------------
! Set up the state-specific basis arrays
!-----------------------------------------------------------------------
      nbas2=0
      do i=1,nsample
         itraj=indx(i,1)
         n=indx(i,2)
         istep=indx(i,3)
         ista=traj(itraj)%ista(n)
         nbas2(ista)=nbas2(ista)+1
         r(ista,nbas2(ista),:)=traj(itraj)%r(n,istep,:)
         p(ista,nbas2(ista),:)=traj(itraj)%p(n,istep,:)
         phase(ista,nbas2(ista))=traj(itraj)%phase(n,istep)         
      enddo

!-----------------------------------------------------------------------
! Open the basis file
!-----------------------------------------------------------------------
      unit=20
      open(unit,file='basis.dat',form='unformatted',status='unknown')

!-----------------------------------------------------------------------
! Write the basis file
!-----------------------------------------------------------------------
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

      ! No. basis functions per electronic state
      do i=1,nsta
         write(unit) nbas2(i)
      enddo

      ! Basis function parameters
      do i=1,nsta
      
         ! Positions
         write(unit) r(i,1:nbas2(i),1:natm*3)

         ! Momenta
         write(unit) p(i,1:nbas2(i),1:natm*3)

         ! Phases
         write(unit) phase(i,1:nbas2(i))

      enddo

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
      write(6,'(/,a)') 'Number of sampled basis functions:'
      do i=1,nsta
         write(6,'(a5,x,i2,a1,x,i6)') 'State',i,':',nbas2(i)
      enddo

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
      allocate(ener(nsta,maxbas2,nsta))
      allocate(nact(nsta,maxbas2,nsta,natm*3))

!-----------------------------------------------------------------------
! Set up the state-specific basis arrays
!-----------------------------------------------------------------------
      nbas2=0
      do i=1,nsample
         itraj=indx(i,1)
         n=indx(i,2)
         istep=indx(i,3)
         ista=traj(itraj)%ista(n)
         nbas2(ista)=nbas2(ista)+1
         ener(ista,nbas2(ista),:)=traj(itraj)%ener(n,istep,:)
         nact(ista,nbas2(ista),:,:)=traj(itraj)%nact(n,istep,:,:)
      enddo
      
!-----------------------------------------------------------------------
! Write the Hamiltonian file
!-----------------------------------------------------------------------
      do i=1,nsta

         ! Energies
         write(unit) ener(i,1:nbas2(i),:)

         ! NACTs
         write(unit) nact(i,1:nbas2(i),:,:)

!         do n=1,nbas2(i)
!            print*,i,n,nact(i,n,:,:)
!         enddo
         
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

  end module tspsgmod


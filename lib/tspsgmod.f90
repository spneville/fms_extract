  module tspsgmod
    
    implicit none
    
    save
    
    integer                            :: maxbas
    integer, dimension(:), allocatable :: nbas

  contains

!#######################################################################

    subroutine tspsg_prep

      use sysdef

      implicit none

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(nbas(nsta))

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
! Determine the number of basis functions per electronic state
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

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      maxbas=maxval(nbas)
      allocate(r(nsta,maxbas,natm*3))
      allocate(p(nsta,maxbas,natm*3))
      allocate(phase(nsta,maxbas))

!-----------------------------------------------------------------------
! Set up the state-specific basis arrays
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
               r(ista,nbas(ista),:)=traj(itraj)%r(n,istep,:)
               p(ista,nbas(ista),:)=traj(itraj)%p(n,istep,:)
               phase(ista,nbas(ista))=traj(itraj)%phase(n,istep)
            enddo
         enddo
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
      
      ! Frozen Gaussian widths
      write(unit) alpha(1:natm*3)

      ! No. basis functions per electronic state
      do i=1,nsta
         write(unit) nbas(i)
      enddo

      ! Basis function parameters
      do i=1,nsta
      
         ! Positions
         write(unit) r(i,1:nbas(i),1:natm*3)

         ! Momenta
         write(unit) p(i,1:nbas(i),1:natm*3)

         ! Phases
         write(unit) phase(i,1:nbas(i))

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
      write(6,'(/,a)') 'Number of basis functions:'
      do i=1,nsta
         write(6,'(a5,x,i2,a1,x,i6)') 'State',i,':',nbas(i)
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
      allocate(ener(nsta,maxbas,nsta))
      allocate(nact(nsta,maxbas,nsta,natm*3))

!-----------------------------------------------------------------------
! Set up the state-specific basis arrays
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
               ener(ista,nbas(ista),:)=traj(itraj)%ener(n,istep,:)
               nact(ista,nbas(ista),:,:)=traj(itraj)%nact(n,istep,:,:)
            enddo
         enddo
      enddo

!-----------------------------------------------------------------------
! Write the Hamiltonian file
!-----------------------------------------------------------------------
      do i=1,nsta

         ! Energies
         write(unit) ener(i,1:nbas(i),:)

         ! NACTs
         write(unit) nact(i,1:nbas(i),:,:)

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


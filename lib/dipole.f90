  module dipole

    implicit none

    integer                                       :: nmaindir,nifg,&
                                                     ntrajdir,ncentdir
    integer, dimension(:), allocatable            :: ifgindx_traj,&
                                                     ifgindx_cent
    complex*16, dimension(:,:), allocatable       :: dipexpec
    character(len=120), dimension(:), allocatable :: amaindir,&
                                                     amaindir_traj,&
                                                     amaindir_cent

  contains

!#######################################################################

    subroutine calc_dipole

      implicit none

!-----------------------------------------------------------------------
! Read the names of the main directories
!-----------------------------------------------------------------------
      call rdmaindirfile

!-----------------------------------------------------------------------
! Determine which of the main directories correspond to trajectories
! and which to centroids
!-----------------------------------------------------------------------
      call getdirtype

!-----------------------------------------------------------------------
! Calculate the dipole expectation values
!-----------------------------------------------------------------------
      call dipole_expec

!-----------------------------------------------------------------------
! Write the dipole expectation values to file
!-----------------------------------------------------------------------
      call wrdipole

      return
      
    end subroutine calc_dipole

!#######################################################################
      
    subroutine rdmaindirfile

      use expec
      use parsemod

      implicit none

      integer                              :: unit,n
      character(len=120)                   :: string

      integer                              :: inkw
      integer, parameter                   :: maxkw=200
      integer, dimension(maxkw)            :: ilkw
      character(len=120), dimension(maxkw) :: keyword

!-----------------------------------------------------------------------
! Open file
!-----------------------------------------------------------------------
      unit=20
      open(unit,file=colppdir_file,form='formatted',status='old')

!-----------------------------------------------------------------------
! Determine the number of main directories and allocate the amaindir
! array
!-----------------------------------------------------------------------
      nmaindir=0
5     continue
      call rdinp(unit,keyword,inkw,ilkw)

      if (keyword(1).ne.'end-file') then
         nmaindir=nmaindir+1
         goto 5
      endif

      allocate(amaindir(nmaindir))

!-----------------------------------------------------------------------
! Read the main directory names
!-----------------------------------------------------------------------
      rewind(unit)
      n=0
10    continue
      call rdinp(unit,keyword,inkw,ilkw)
      if (keyword(1).ne.'end-file') then
         n=n+1
         amaindir(n)=keyword(1)
         goto 10
      endif

!-----------------------------------------------------------------------
! Close file
!-----------------------------------------------------------------------
      close(unit)

      return

    end subroutine rdmaindirfile

!#######################################################################

    subroutine getdirtype

      implicit none

      integer :: i,kt,kc

!-----------------------------------------------------------------------
! Get the numbers of trajectory and centroid directories
!-----------------------------------------------------------------------
      ntrajdir=0
      ncentdir=0

      do i=1,nmaindir
         if (index(amaindir(i),'_traj').ne.0) then
            ntrajdir=ntrajdir+1
         else if (index(amaindir(i),'_cent').ne.0) then
            ncentdir=ncentdir+1
         endif
      enddo

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(amaindir_traj(ntrajdir))
      allocate(amaindir_cent(ncentdir))

!-----------------------------------------------------------------------
! Fill in the amaindir_traj and amaindir_cent arrays
!-----------------------------------------------------------------------
      kt=0
      kc=0
      do i=1,nmaindir
         if (index(amaindir(i),'_traj').ne.0) then
            kt=kt+1
            amaindir_traj(kt)=amaindir(i)
         else if (index(amaindir(i),'_cent').ne.0) then
            kc=kc+1
            amaindir_cent(kc)=amaindir(i)
         endif
      enddo

      return

    end subroutine getdirtype

!#######################################################################

    subroutine dipole_expec

      use expec
      use sysdef

      implicit none

      integer                                       :: i,itmp,nsubdir,&
                                                       tindx1,tindx2
      integer, dimension(:), allocatable            :: step,icontrib
      real*8, dimension(:,:,:,:), allocatable       :: dipmat
      complex*16, dimension(:), allocatable         :: coeff
      character(len=120), dimension(:), allocatable :: asubdir

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      itmp=nstep/dstep+1
      allocate(dipexpec(3,itmp))
      dipexpec=(0.0d0,0.0d0)

!-----------------------------------------------------------------------
! Determine the number of IFGs being considered, which may be less than
! nintraj
!-----------------------------------------------------------------------
      call getnifg

!-----------------------------------------------------------------------
! Trajectory contributions
!-----------------------------------------------------------------------
      ! Loop over the main trajectory directories
      do i=1,ntrajdir

         ! Ouput our progress
         write(6,'(2a)') 'Processing directory: ',trim(amaindir_traj(i))

         ! Get the list of timesteps/subdirectories
         call getsubdirs(amaindir_traj(i),nsubdir,asubdir,step)

         ! Read the coefficients for each timestep/subdirectory
         call getcoeff(coeff,amaindir_traj(i),asubdir,nsubdir)

         ! Read the dipole matrix elements for each timestep/subdirectory
         call rddipole(amaindir_traj(i),asubdir,nsubdir,dipmat)

         ! Determine the trajectory number
         call get_trajnum(tindx1,amaindir_traj(i))

         ! Nuclear contribution to the dipole expectation value
         call nuc_contrib(ifgindx_traj(i),tindx1,tindx1,nsubdir,step)
         
         ! Electronic contribution to the dipole expectation value
         call elec_contrib(ifgindx_traj(i),tindx1,tindx1,nsubdir,&
              step,dipmat)
         
      enddo

      ! Averaging over the IFGs
      dipexpec=dipexpec/nifg

      return

    end subroutine dipole_expec

!#######################################################################

    subroutine getnifg

      use sysdef

      implicit none

      integer                     :: i,k,unit
      integer, dimension(nintraj) :: cnt
      character(len=130)          :: ain
             
!-----------------------------------------------------------------------
! (1) Trajectories
!-----------------------------------------------------------------------
      allocate(ifgindx_traj(ntrajdir))

      unit=20
      cnt=0
      do i=1,ntrajdir
         ain=trim(amaindir_traj(i))//'/ifg_number'
         open(unit,file=ain,form='formatted',status='old')
         read(unit,*) k
         cnt(k)=cnt(k)+1
         ifgindx_traj(i)=k
         close(unit)
      enddo

!-----------------------------------------------------------------------
! (2) Centroids
!-----------------------------------------------------------------------
      allocate(ifgindx_cent(ncentdir))

      do i=1,ncentdir
         ain=trim(amaindir_cent(i))//'/ifg_number'
         open(unit,file=ain,form='formatted',status='old')
         read(unit,*) k
         cnt(k)=cnt(k)+1
         ifgindx_cent(i)=k
         close(unit)
      enddo

!-----------------------------------------------------------------------
! Determine the no. contributing IFGs
!-----------------------------------------------------------------------
      nifg=0
      do i=1,nintraj
         if (cnt(i).gt.0) nifg=nifg+1
      enddo

      return

    end subroutine getnifg

!#######################################################################

    subroutine getsubdirs(amaindir,nsubdir,asubdir,step)

      implicit none

      integer                                       :: unit,nsubdir,i
      integer, dimension(:), allocatable            :: step
      character(len=120)                            :: amaindir
      character(len=130)                            :: alist,string
      character(len=120), dimension(:), allocatable :: asubdir
      logical(kind=4)                               :: lexists
      
!-----------------------------------------------------------------------
! Read the sublist file
!-----------------------------------------------------------------------
      alist=trim(amaindir)//'/sublist'

      inquire(file=alist,exist=lexists)
      if (.not.lexists) then
         write(6,'(2a)') 'sublist file not found in the directory ',&
              trim(amaindir)
         STOP
      endif

      unit=20
      open(unit,file=alist,form='formatted',status='old')

      nsubdir=0
10    continue
      read(unit,'(a)',end=100) string
      nsubdir=nsubdir+1
      goto 10

100   continue
      if (allocated(asubdir)) deallocate(asubdir)
      allocate(asubdir(nsubdir))
      
      rewind(unit)
      do i=1,nsubdir
         read(unit,'(a)') asubdir(i)
      enddo

!-----------------------------------------------------------------------
! Determine the timesteps
!-----------------------------------------------------------------------
      if (allocated(step)) deallocate(step)
      allocate(step(nsubdir))

      do i=1,nsubdir
         read(asubdir(i)(5:9),*) step(i)
      enddo

      close(unit)

      return

    end subroutine getsubdirs 

!#######################################################################

    subroutine getcoeff(coeff,amaindir,asubdir,nsubdir)

      implicit none

      integer                                :: nsubdir,i,unit
      real*8                                 :: cr,ci
      complex*16, dimension(:), allocatable  :: coeff
      character(len=120)                     :: amaindir
      character(len=120), dimension(nsubdir) :: asubdir
      character(len=270)                     :: acoeff

!-----------------------------------------------------------------------
! Allocate the coeff array
!-----------------------------------------------------------------------
      if (allocated(coeff)) deallocate(coeff)
      allocate(coeff(nsubdir))
      coeff=0.0d0

!-----------------------------------------------------------------------
! Loop over timesteps/subdirectories and read the coefficents for each
!-----------------------------------------------------------------------
      unit=20
      do i=1,nsubdir
         acoeff=trim(amaindir)//'/'//trim(asubdir(i))//'coeff'
         open(unit,file=acoeff,form='formatted',status='old')
         read(unit,*)
         read(unit,*) cr,ci
         coeff(i)=dcmplx(cr,ci)
         close(unit)
      enddo

      return

    end subroutine getcoeff

!#######################################################################

    subroutine rddipole(amaindir,asubdir,nsubdir,dipmat)

      use sysdef

      implicit none
      
      integer                                 :: unit,nsubdir,i,j,k,&
                                                 s1,s2
      real*8, dimension(:,:,:,:), allocatable :: dipmat
      character(len=120)                      :: amaindir,string,atmp
      character(len=120), dimension(nsubdir)  :: asubdir
      character(len=300)                      :: filename

!-----------------------------------------------------------------------
! Allocate the dipole matrix array
!-----------------------------------------------------------------------
      if (allocated(dipmat)) deallocate(dipmat)
      allocate(dipmat(nsubdir,nsta,nsta,3))
      dipmat=0.0d0

!-----------------------------------------------------------------------
! Read the dipole matrix elements for each timestep/subdirectory
!-----------------------------------------------------------------------
      unit=20
      
      ! Loop over timesteps/subdirectories
      do i=1,nsubdir

         ! Open the dipole matrix file
         filename=trim(amaindir)//'/'//trim(asubdir(i)) &
              //'columbus/dipole.dat'
         open(unit,file=filename,form='formatted',status='old')
         
         ! Read the dipole matrix file
10       read(unit,'(a)',end=999) string
         
         ! On-diagonal elements
         if (index(string,'State ').ne.0) then
            ! State number
            read(string(6:),*) s1
            ! x,y,z components
            do k=1,4
               read(unit,*)
            enddo
            read(unit,*) atmp,(dipmat(i,s1,s1,j), j=1,3)

         ! Off-diagonal elements
         else if (index(string,'States:').ne.0) then

            ! State numbers
            read(string(8:),*) s1,s2 

            ! x,y,z components
            do k=1,4
               read(unit,*)
            enddo
            read(unit,*) atmp,(dipmat(i,s1,s1,j), j=1,3)

         endif

         goto 10

999      continue

         ! Close the dipole matrix file
         close(unit)

      enddo

      return

    end subroutine rddipole

!#######################################################################

    subroutine get_trajnum(tindx,amaindir)

      implicit none

      integer            :: tindx,k
      character(len=120) :: amaindir

      k=index(amaindir,'_traj')
      read(amaindir(k+5:k+6),'(i2)') tindx

      return

    end subroutine get_trajnum

!#######################################################################

    subroutine nuc_contrib(ifg,tindx1,tindx2,nsubdir,step)

      use expec
      use sysdef
      use trajdef
      use gausstools

      implicit none

      integer                     :: ifg,tindx1,tindx2,nsubdir,i,k
      integer, dimension(nsubdir) :: step
      complex*16, dimension(3)    :: gmug
      complex*16                  :: c1,c2,c12
      
      ! Loop over timesteps/subdirectories
      do i=1,nsubdir

         ! Calculate the nuclear dipole matrix elements
         ! <g_j | mu_N,c | g_k >, c=x,y,z, for the current timestep
         gmug=nucdip_integral(ifg,tindx1,tindx2,step(i))
         
         ! Coefficients
         c1=traj(ifg)%coe(tindx1,step(i))
         c2=traj(ifg)%coe(tindx2,step(i))
         c12=conjg(c1)*c2

         ! Add the current contribution to the dipole expectation
         ! value
         !k=1+(step(i)-1)/dstep         
         k=step(i)/dstep+1

         dipexpec(:,k)=dipexpec(:,k)+c12*gmug

      enddo
      
      return

    end subroutine nuc_contrib

!#######################################################################

    subroutine elec_contrib(ifg,tindx1,tindx2,nsubdir,step,dipmat)

      use expec
      use sysdef
      use trajdef
      use gausstools

      implicit none

      integer                                :: ifg,tindx1,tindx2,&
                                                nsubdir,i,k,s1,s2
      integer, dimension(nsubdir)            :: step
      real*8, dimension(nsubdir,nsta,nsta,3) :: dipmat
      complex*16                             :: c1,c2,c12,ovrlp

!-----------------------------------------------------------------------
! First-order saddlepoint approximation to the contributions of the
! electronic dipole operator to the total dipole expectation value
!-----------------------------------------------------------------------
      ! Loop over timesteps/subdirectories
      do i=1,nsubdir

         ! Overlap <g1|g2>
         if (tindx1.eq.tindx2) then
            ovrlp=(1.0d0,0.0d0)
         else
            ovrlp=overlap_general_nuc_only(ifg,ifg,tindx1,tindx2,&
                 step(i),step(i))
         endif

         ! Coefficients
         c1=traj(ifg)%coe(tindx1,step(i))
         c2=traj(ifg)%coe(tindx2,step(i))
         c12=conjg(c1)*c2

         ! State indices
         s1=traj(ifg)%ista(tindx1)
         s2=traj(ifg)%ista(tindx2)

         ! Add the current contribution to the dipole expectation
         ! value
         !k=1+(step(i)-1)/dstep
         k=step(i)/dstep+1

         dipexpec(:,k)=dipexpec(:,k)+c12*ovrlp*dipmat(i,s1,s2,:)
         
      enddo

      return

    end subroutine elec_contrib

!#######################################################################

    subroutine wrdipole

      use sysdef
      use expec

      implicit none

      integer :: unit,i,j,indx

!-----------------------------------------------------------------------
! Open the output file
!-----------------------------------------------------------------------
      unit=20
      open(unit,file='dipole.dat',form='formatted',status='unknown')

!-----------------------------------------------------------------------
! Write the dipole expectation values to file
!-----------------------------------------------------------------------
      write(unit,'(64a)') ('#', i=1,64)
      write(unit,'(a)') '# t (fs)         <mu_x>            &
           <mu_y>           <mu_z>'
      write(unit,'(64a)') ('#', i=1,64)

      do i=1,nstep,dstep
         
         indx=i/dstep+1

         write(unit,'(F8.3,3(3x,E15.7))') (i-1)*dt/41.341375d0,&
              (real(dipexpec(j,indx)), j=1,3)

      enddo

!-----------------------------------------------------------------------
! Close the output file
!-----------------------------------------------------------------------
      close(unit)

      return

    end subroutine wrdipole

!#######################################################################

  end module dipole

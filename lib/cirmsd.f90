  module cirmsd

  contains

!#######################################################################

    subroutine ci_spawn_rmsd
    
      use sysdef
      use expec, only: cifile

      implicit none

      integer                              :: nspawn,nci
      integer, dimension(:,:), allocatable :: spawnstates
      integer, dimension(:), allocatable   :: spawninc
      real*8, dimension(:,:), allocatable  :: xspawn,xci
      logical                              :: lproj

!-----------------------------------------------------------------------
! Read the spawn geometries from file
!-----------------------------------------------------------------------
      call rdspawnfile(xspawn,spawnstates,nspawn)

!-----------------------------------------------------------------------
! Read the conical intersection geometries from file
!-----------------------------------------------------------------------
      call rdcifile(xci,nci)

!-----------------------------------------------------------------------
! Determine whether the RMSD within the branching space is to be
! calculated, and, if so, calculated the projector onto the branching
! space
!-----------------------------------------------------------------------
      call is_proj(lproj,nci)
      if (lproj) call getbranchproj

!-----------------------------------------------------------------------
! Determine which spawn events are between the states of interest
!-----------------------------------------------------------------------
      call getspawninc(spawninc,nspawn,spawnstates)

!-----------------------------------------------------------------------
! Calculate the RMSD between each pair of spawn and conical intersection
! geometries using the Kabsch algorithm
!-----------------------------------------------------------------------
      call calc_rmsd(xspawn,xci,nspawn,nci,spawnstates,lproj,spawninc)

      return

    end subroutine ci_spawn_rmsd

!#######################################################################

    subroutine is_proj(lproj,nci)

      use expec

      implicit none

      integer :: nci
      logical :: lproj

      lproj=.false.

      if (hfile.ne.'') then
         lproj=.true.
         if (gfile(1).eq.''.or.gfile(2).eq.'') then
            write(6,'(/,2x,a,/)') 'The RMSD within the branching space &
                 has been requested, but not all gradient files have &
                 been given'
            STOP
         endif
      endif

      if (gfile(1).ne.'') then
         lproj=.true.
         if (hfile.eq.'') then
            write(6,'(/,2x,a,/)') 'The RMSD within the branching space &
                 has been requested, but the NACT file has not been &
                 given'
            STOP
         endif
      endif

      if (lproj.and.nci.gt.1) then
         write(6,'(/,2x,a,/)') 'The calculation of the RMSD within &
              the branching space is only currently supported for a &
              single CI geometry'
         STOP
      endif

      return

    end subroutine is_proj

!#######################################################################

    subroutine rdspawnfile(xspawn,spawnstates,nspawn)

      use sysdef

      implicit none

      integer                              :: unit,nspawn,i,j,k
      integer, dimension(:,:), allocatable :: spawnstates
      real*8, dimension(:,:), allocatable  :: xspawn

!-----------------------------------------------------------------------
! Determine the number of spawn geometries and allocate the xspawn
! array
!-----------------------------------------------------------------------
      unit=20
      open(unit,file='spawngeom.xyz',form='formatted',status='old')
      
      nspawn=0
10    continue
      do i=1,natm+2
         read(unit,*,end=15)
      enddo
      nspawn=nspawn+1
      goto 10

15    continue
      
      allocate(xspawn(nspawn,natm*3))
      allocate(spawnstates(nspawn,2))

!-----------------------------------------------------------------------
! Read the spawn geometries from file
!-----------------------------------------------------------------------
      rewind(unit)

      do i=1,nspawn
         read(unit,*)
         read(unit,'(23x,i2,27x,i2)') spawnstates(i,2),spawnstates(i,1)
         do j=1,natm
            read(unit,'(2x,3(3x,F13.9))') (xspawn(i,k), k=j*3-2,j*3)
         enddo
      enddo

      close(unit)

      return

    end subroutine rdspawnfile

!#######################################################################

    subroutine rdcifile(xci,nci)

      use sysdef
      use expec, only: cifile

      implicit none

      integer                             :: unit,nci,i,j,k
      real*8, dimension(:,:), allocatable :: xci
      character(len=10)                   :: atmp

!-----------------------------------------------------------------------
! Determine the number of CIs geometries and allocate the xci array
!-----------------------------------------------------------------------
      unit=20
      open(unit,file=cifile,form='formatted',status='old')
      
      nci=0
10    continue
      do i=1,natm+2
         read(unit,*,end=15)
      enddo
      nci=nci+1
      goto 10

15    continue
      
      allocate(xci(nci,natm*3))

!-----------------------------------------------------------------------
! Read the CI geometries from file
!-----------------------------------------------------------------------
      rewind(unit)

      do i=1,nci
         read(unit,*)
         read(unit,*)
         do j=1,natm
            read(unit,*) atmp,(xci(i,k), k=j*3-2,j*3)
         enddo
      enddo

      close(unit)

      return

    end subroutine rdcifile

!#######################################################################

    subroutine calc_rmsd(xspawn,xci,nspawn,nci,spawnstates,lproj,&
         spawninc)

      use expec
      use sysdef
      use kabschmod
      use permutemod
      use mathlib

      implicit none

      integer                          :: nspawn,nci,i,j,k,imin,nother,unit
      integer, dimension(nspawn,2)     :: spawnstates
      integer, dimension(nspawn)       :: spawninc
      integer, dimension(nci)          :: nmatch
      real*8, dimension(nspawn,natm*3) :: xspawn
      real*8, dimension(nci,natm*3)    :: xci
      real*8, dimension(natm*3)        :: currgeom,spawngeom,rotx,xcurr,&
                                          xmin,dvec
      real*8, dimension(nspawn,nci)    :: rmsd
      real*8                           :: fmin,pc,sum,currrmsd,minrmsd
      real*8, parameter                :: tol=2.5d0
      character(len=120)               :: atmp
      logical                          :: lproj
      
      integer, dimension(npermute)         :: indx
      integer, dimension(natm)             :: is_perm
      integer, dimension(:,:), allocatable :: P
      integer                              :: n,ncomb,count,spawnindx,&
                                              ciindx

      nmatch=0
      nother=0

      if (npermute.eq.0) then

         ! Loop over spawn geometries
         do i=1,nspawn
            if (spawninc(spawnindx).eq.0) cycle
            fmin=9999999.9d0
            ! Loop over ci geometries
            do j=1,nci
               ! Put the current spawn geometry into maximum coincidence
               ! with the current ci geometry
               currgeom=xci(j,:)
               spawngeom=xspawn(i,:)
               call kabsch(currgeom,spawngeom,rmsd(i,j),rotx)
               ! If the current ci geometry is the closest yet to the
               ! curret spawn geometry, then set imin=j
               if (rmsd(i,j).lt.fmin.and.rmsd(i,j).lt.tol) then
                  imin=j
                  fmin=rmsd(i,j)
               endif
            enddo
            if (fmin.lt.tol) then
               nmatch(imin)=nmatch(imin)+1
            else
               nother=nother+1
            endif
         enddo

      else
         
         ! Generate all permutations of the indices of the nuclei being
         ! permuted
         ncomb=factorial(npermute)
         allocate(P(ncomb,npermute))
         call permutate(pindx,P)
         ! Determine the indices of the nuclei that are to be permuted
         is_perm=0
         do i=1,natm
            do j=1,npermute
               if (pindx(j).eq.i) is_perm(i)=1
            enddo
         enddo
         ! Loop over spawn geometries
         do spawnindx=1,nspawn            
            if (spawninc(spawnindx).eq.0) cycle
            fmin=9999999.9d0
            ! Loop over ci geometries
            do ciindx=1,nci
               ! Set the current CI and spawn geometries
               currgeom=xci(ciindx,:)
               spawngeom=xspawn(spawnindx,:)               
               ! Loop over all permutations of the nuclei being permuted,
               ! calculate the RMSD from the CI geometry for each and save
               ! the geometry with the lowest RMSD
               minrmsd=1000000.0d0
               do i=1,ncomb
                  indx=P(i,:)
                  xcurr=currgeom
                  count=0
                  do j=1,natm
                     if (is_perm(j).eq.1) then
                        count=count+1
                        xcurr(j*3-2:j*3)=spawngeom(indx(count)*3-2:indx(count)*3)
                     endif
                  enddo                  
                  call kabsch(currgeom,xcurr,currrmsd,rotx)
                  if (currrmsd.lt.minrmsd) then
                     minrmsd=currrmsd
                     xmin=rotx
                  endif
               enddo

               ! If the RMSD within the branching space has been
               ! requested, then reset minrmsd to this
               if (lproj) then
                  dvec=currgeom-xmin
                  dvec=matmul(branchproj,dvec)
                  minrmsd=sqrt(dot_product(dvec,dvec)/real(natm))
               endif

               ! Set the minimum RMSD
               rmsd(spawnindx,ciindx)=minrmsd
               if (rmsd(spawnindx,ciindx).lt.fmin&
                    .and.rmsd(spawnindx,ciindx).lt.tol) then
                  imin=ciindx
                  fmin=rmsd(spawnindx,ciindx)
               endif
            enddo
            if (fmin.lt.tol) then
               nmatch(imin)=nmatch(imin)+1
            else
               nother=nother+1
            endif
         enddo

      endif

      ! Output results      
      write(6,*)
      sum=0.0d0
      do i=1,nci
         pc=100.0d0*real(nmatch(i))/real(nspawn)
         sum=sum+pc
         write(6,'(2x,a,x,i2,a1,2x,i3,x,a1,i4,a2)') &
              'Number of spawn geometries closest to CI geometry',&
              i,':',nmatch(i),'(',nint(pc),'%)'
      enddo
      write(6,'(2x,a,38x,i3,x,a1,i3,a2)') &
              'Other geometries:',nother,'(',nint(100.0d0-sum),'%)'
      write(6,*)

      ! Write the gnuplot file
      unit=20
      open(unit,file='cirmsd.gnu',form='formatted',status='unknown')
      atmp='# Units: Angstrom'
      write(unit,'(a)') adjustl(trim(atmp))
      atmp='# Spawn geometry | RMSD_1 | ... | RMSD_NCI'
      write(unit,'(a,/)') adjustl(trim(atmp))
      write(unit,'(a)') 'set xlabel "RMSD (Angstrom)"'
      write(unit,'(a)') 'set ylabel "Number of spawn geometries"'
      write(unit,'(a)')'set style fill solid border -1'
      write(unit,'(a)')'n=40 #number of intervals'
      write(unit,'(a)')'max=1.0 #max value'
      write(unit,'(a)')'min=0.0 #min value'
      write(unit,'(a)')'width=(max-min)/n #interval width'
      write(unit,'(a)')'#function used to map a value to the intervals'
      write(unit,'(a)')'hist(x,width)=width*floor(x/width)+width/2.0'
      write(unit,'(a)')'set boxwidth width*0.9'
      write(unit,'(a)')'plot "-" u (hist($2,width)):(1.0) smooth freq w boxes lc rgb "blue" notitle'
      do i=1,nspawn
         if (spawninc(i).eq.0) cycle
         write(unit,*) i,(rmsd(i,j),j=1,nci)
      enddo
      write(unit,'(a)') 'e'
      write(unit,'(a)') 'pause -1'
      close(unit)

      return

    end subroutine calc_rmsd

!#######################################################################

    subroutine getbranchproj

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
      allocate(branchvec(2,3*natm))

!-----------------------------------------------------------------------
! NACT vector
!-----------------------------------------------------------------------
      unit=358
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

!-----------------------------------------------------------------------
! Construction of the projector onto the branching space
!-----------------------------------------------------------------------
      allocate(branchproj(3*natm,3*natm))
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
      
    end subroutine getbranchproj

!#######################################################################

    subroutine getspawninc(spawninc,nspawn,spawnstates)

      use expec

      implicit none

      integer, dimension(:), allocatable :: spawninc
      integer, dimension(nspawn,2)       :: spawnstates
      integer                            :: nspawn,i

      allocate(spawninc(nspawn))

      spawninc=0

      if (cirmsdsta(1).eq.0) then
         spawninc=1
      else
         do i=1,nspawn
            if (spawnstates(i,1).eq.cirmsdsta(1)&
                 .and.spawnstates(i,2).eq.cirmsdsta(2)) spawninc(i)=1
            if (spawnstates(i,2).eq.cirmsdsta(1)&
                 .and.spawnstates(i,1).eq.cirmsdsta(2)) spawninc(i)=1
         enddo
      endif

      return

    end subroutine getspawninc

!#######################################################################

  end module cirmsd

!#######################################################################
! Subroutines and functions for the manipulation of geometries.
! Mainly to do with putting geometries into maximum coincidence wrt
! potential energy preserving operations.
!#######################################################################

  module geomtools

    implicit none

  contains
    
!#######################################################################
! kabsch: puts the Cartesian coordinates x2 into maximum coincidence
!         with the Cartesian coordinates x1 using the Kabsch algorithm
!#######################################################################
    
    subroutine kabsch(x1,x2,rmsd,rotx)

      use sysdef

      implicit none

      integer                          :: i,j,k,l
      real*8, dimension(natm*3)        :: x1,x2,rotx
      real*8, dimension(natm,3)        :: x1mat,x2mat,rmat
      real*8, dimension(3,3)           :: vut,rotmat,tmpmat
      real*8, dimension(3)             :: com1,com2,pm
      real*8                           :: totmass,det,d,rmsd,min

      integer                          :: lwork,info,dim
      double precision, dimension(3)   :: sigma
      double precision, dimension(3,3) :: covar,u,vt,orig
      double precision, dimension(15)  :: work

!-----------------------------------------------------------------------
! Translate the origins of both geometries to the centre of mass
!-----------------------------------------------------------------------
      totmass=0.0d0
      com1=0.0d0
      com2=0.0d0
      do i=1,natm
         totmass=totmass+atmass(i)
         do j=1,3
            com1(j)=com1(j)+x1(i*3-3+j)*atmass(i)
            com2(j)=com2(j)+x2(i*3-3+j)*atmass(i)
         enddo
      enddo

      com1=com1/totmass
      com2=com2/totmass

      do i=1,natm
         do j=1,3
            x1(i*3-3+j)=x1(i*3-3+j)-com1(j)
            x2(i*3-3+j)=x2(i*3-3+j)-com2(j)
         enddo
      enddo

      do i=1,natm
         do j=1,3
            x1mat(i,j)=x1(i*3-3+j)
            x2mat(i,j)=x2(i*3-3+j)
         enddo
      enddo

!-----------------------------------------------------------------------
! Calculate the 3X3 covariance matrix
!-----------------------------------------------------------------------
      covar=0.0d0
      do i=1,3
         do j=1,3
            do k=1,natm
               covar(i,j)=covar(i,j)+x1(k*3-3+i)*x2(k*3-3+j)
            enddo
         enddo
      enddo

!-----------------------------------------------------------------------
! Calculate the SVD of the covariance matrix
! N.B., This DOES work
!-----------------------------------------------------------------------
      orig=covar
      lwork=15
      dim=3

      call dgesvd('A','A',dim,dim,covar,dim,sigma,u,dim,vt,dim,work,lwork,info)

      covar=orig

      if (info.ne.0) then
         write(6,'(/,2x,a,/)') 'SVD of the covariance matrix in &
              subroutine kabsch failed'
         STOP
      endif

!-----------------------------------------------------------------------
! Calculate the rotation matrix
!-----------------------------------------------------------------------
      vut=0.0d0
      do i=1,3
         do j=1,3
            do k=1,3
               vut(i,j)=vut(i,j)+vt(k,i)*u(j,k)
            enddo
         enddo
      enddo

      det=finddet(covar,3)

      if (det.lt.0.0d0) then
         d=-1.0d0
      else
         d=1.0d0
      endif

      pm(1)=1.0d0
      pm(2)=1.0d0
      pm(3)=d

      rotmat=0.0d0
      do i=1,3
         do j=1,3
            do k=1,3
               rotmat(i,j)=rotmat(i,j)+vt(k,i)*pm(k)*u(j,k)
            enddo
         enddo
      enddo

!-----------------------------------------------------------------------
! For debugging purposes, calculate the minimum possible residual
!-----------------------------------------------------------------------
      min=0.0d0
      do k=1,natm
         do j=1,3
            min=min+x1mat(k,j)**2+x2mat(k,j)**2
         enddo
      enddo
      min=(min/natm)-(2.0d0/natm)*(sigma(1)+sigma(2)+sigma(3)*d)

!-----------------------------------------------------------------------
! Rotate the geometry x2
!-----------------------------------------------------------------------
      do i=1,natm
         do j=1,3
            rotx(i*3-3+j)=0.0
            do k=1,3
               rotx(i*3-3+j)=rotx(i*3-3+j)+&
                    rotmat(k,j)*x2(i*3-3+k)
            enddo
         enddo
      enddo

      rmsd=0.0d0
      do i=1,natm*3
         rmsd=rmsd+(rotx(i)-x1(i))**2
      enddo
      rmsd=rmsd/natm
      rmsd=sqrt(rmsd)

      return

    end subroutine kabsch

!#######################################################################
! finddet: returns the determinant of a matrix
!#######################################################################

    function finddet(matrix,n)

      implicit none

      integer*4              :: n,i,j,k,l
      real*8, dimension(n,n) :: matrix
      real*8                 :: m,temp,finddet
      logical(kind=4)        :: detexists = .TRUE.      

      l=1

      ! Convert to upper triangular form
      do k=1,n-1
         if (matrix(k,k).eq.0) then
            detexists = .false.
            do i=k+1,n
               if (matrix(i,k).ne.0) then
                  do j=1,n
                     temp=matrix(i,j)
                     matrix(i,j)=matrix(k,j)
                     matrix(k,j)=temp
                  enddo
                  detexists=.true.
                  l=-l
                  exit
               endif
            enddo
            if (.not.detexists) then
               finddet=0
               return
            endif
         endif
         do j=k+1,n
            m=matrix(j,k)/matrix(k,k)
            do i=k+1,n
               matrix(j,i)=matrix(j,i)-m*matrix(k,i)
            enddo
         enddo
      enddo
      
      ! Calculate determinant by finding product of diagonal elements
      finddet=l
      do i=1,n
         finddet=finddet*matrix(i,i)
      enddo
   
    end function finddet

!#######################################################################
! permutate: for the array of indices E, returns the array P of all
!            possible permutations of these indices
!#######################################################################
    
      recursive subroutine permutate(E,P) 

        implicit none
        
        integer, intent(in)  :: E(:)       ! array of objects 
        integer, intent(out) :: P(:,:)     ! permutations of E 
        integer              :: N, Nfac, i, k, &
                                S(size(P,1)/size(E), size(E)-1) 

        N = size(E); Nfac = size(P,1); 
        
        do i=1,N                           ! cases with E(i) in front 
           if( N>1 ) call permutate((/E(:i-1), E(i+1:)/), S) 
           forall(k=1:Nfac/N) P((i-1)*Nfac/N+k,:) = (/E(i), S(k,:)/) 
        end do
        
      end subroutine permutate

!#######################################################################
! maxcoinc: puts the Cartesian coordinates x2 into maximum coincidence
!           with the Cartesian coordinates x1 via the Kabsch
!           algorithm and optionally the permutation of a set of
!           identical nuclei
!#######################################################################
      
      function maxcoinc(x1,x2) result(xmax)

        use sysdef
        use expec
        use mathlib

        integer                              :: i,j
        integer, dimension(npermute)         :: indx
        integer, dimension(natm)             :: is_perm
        integer, dimension(:,:), allocatable :: P
        integer                              :: n,ncomb,count
        real*8, dimension(natm*3)            :: x1,x2,xmax,rotx,xcurr,&
                                                xtmp
        real*8                               :: rmsd,minrmsd
        
!-----------------------------------------------------------------------
! Rotations and centre of mass translations only
!-----------------------------------------------------------------------        
        if (npermute.eq.0) then
           call kabsch(x1,x2,rmsd,xmax)
        endif

!-----------------------------------------------------------------------        
! Rotations and centre of mass translations + permutations of a set
! of identical nuclei
!-----------------------------------------------------------------------           
        if (npermute.gt.0) then

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

           ! Loop over all permutations of the nuclei being permuted,
           ! calculate the RMSD from the reference geometry x1 for each and save
           ! the geometry with the lowest RMSD
           minrmsd=1e+6
           do i=1,ncomb

              indx=P(i,:)
              xcurr=x2

              count=0
              do j=1,natm
                 if (is_perm(j).eq.1) then
                    count=count+1
                    xcurr(j*3-2:j*3)=x2(indx(count)*3-2:indx(count)*3)
                 endif
              enddo

              call kabsch(x1,xcurr,rmsd,xtmp)              

              if (rmsd.lt.minrmsd) then
                 minrmsd=rmsd
                 xmax=xtmp
              endif

           enddo
              
        endif

        return
        
      end function maxcoinc
        
!#######################################################################
    
  end module geomtools

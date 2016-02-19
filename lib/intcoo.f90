  module intcoo

    implicit none

  contains

!#######################################################################
! x2int: for a given set of Cartesian coordinates xcoo returns the
!        value of a single internal coordinate intcoo
!
! ityp=1 <-> bond length
!      2 <-> angle
!      3 <-> dihedral
!      4 <-> twisting angle
!     -1 <-> Cartesian vector
!#######################################################################

    function x2int(xcoo) result(intcoo)
      
      use expec
      use sysdef

      implicit none

      real*8, dimension(natm*3) :: xcoo
      real*8                    :: intcoo

      select case(ityp)

      case(1) ! Bond length
         call calc_length(xcoo,intcoo)

      case(2) ! Bond angle
         call calc_angle(xcoo,intcoo)

      case(3) ! Dihedral angle
         call calc_dihedral(xcoo,intcoo)

      case(4) ! Twisting angle
         call calc_twist(xcoo,intcoo)

      end select

      return

    end function x2int

!#######################################################################

    function x2int_mom(pcoo,xcoo) result(intmom)

      use expec
      use sysdef

      implicit none
      
      real*8, dimension(natm*3) :: pcoo,xcoo
      real*8                    :: intmom

      select case(ityp)

      case(-1) ! Cartesian vector
         call calc_cart(pcoo,xcoo,intmom)
         return
      end select
         
      write(6,'(/,2x,a,/)') 'The selected nternal coordinate type is &
           not currently supported in the momentum representation'
      STOP

    end function x2int_mom
    
!#######################################################################

    subroutine calc_length(xcoo,intcoo)

      use expec
      use sysdef

      implicit none

      integer*8                 :: i,k1,k2
      real*8, dimension(natm*3) :: xcoo
      real*8                    :: intcoo
      
      intcoo=0.0d0

      k1=iatm(1)
      k2=iatm(2)
      
      do i=1,3
         intcoo=intcoo+(xcoo(k1*3-3+i)-xcoo(k2*3-3+i))**2
      enddo

      intcoo=sqrt(intcoo)
      
      return

    end subroutine calc_length

!#######################################################################

    subroutine calc_angle(xcoo,intcoo)

      use expec
      use sysdef

      implicit none

      integer*8                 :: i,k1,k2,k3
      real*8, dimension(natm*3) :: xcoo
      real*8                    :: intcoo,dp,len1,len2,pi
      real*8, dimension(3)      :: vec1,vec2

      pi=4.0d0*datan(1.0d0)

      k1=iatm(1)
      k2=iatm(2)
      k3=iatm(3)

      len1=0.0d0
      len2=0.0d0

      do i=1,3
         vec1(i)=xcoo(k2*3-3+i)-xcoo(k1*3-3+i)         
         vec2(i)=xcoo(k2*3-3+i)-xcoo(k3*3-3+i)
         len1=len1+vec1(i)**2
         len2=len2+vec2(i)**2
      enddo

      len1=sqrt(len1)
      len2=sqrt(len2)

      dp=dot_product(vec1,vec2)

      intcoo=dp/(len1*len2)
      intcoo=dacos(intcoo)
      intcoo=intcoo*180.0d0/pi

      return

    end subroutine calc_angle

!#######################################################################

    subroutine calc_dihedral(xcoo,intcoo)

      use expec
      use sysdef

      implicit none

      integer*8                 :: i,k1,k2,k3,k4
      real*8, dimension(natm*3) :: xcoo
      real*8                    :: intcoo,pi,ftmp1,ftmp2
      real*8, dimension(3)      :: vec1,vec2,vec3,tmp1,tmp2,tmp3,tmp4

      pi=4.0d0*datan(1.0d0)

      k1=iatm(1)
      k2=iatm(2)
      k3=iatm(3)
      k4=iatm(4)

      do i=1,3
         vec1(i)=xcoo(k2*3-3+i)-xcoo(k1*3-3+i)
         vec2(i)=xcoo(k3*3-3+i)-xcoo(k2*3-3+i)
         vec3(i)=xcoo(k4*3-3+i)-xcoo(k3*3-3+i)
      enddo

      tmp1=cross_product(vec1,vec2)
      tmp2=cross_product(vec2,vec3)
      
      ftmp1=dot_product(vec2,vec2)
      ftmp1=sqrt(ftmp1)
      tmp3=vec2/ftmp1
      
      tmp4=cross_product(tmp1,tmp2)
      ftmp1=dot_product(tmp4,tmp3)

      ftmp2=dot_product(tmp1,tmp2)

      intcoo=atan2(ftmp1,ftmp2)

      intcoo=intcoo*180.0d0/pi

      return

    end subroutine calc_dihedral

!#######################################################################

    subroutine calc_twist(xcoo,intcoo)

      use expec
      use sysdef
      
      integer*8                 :: i,k1,k2,k3,k4,k5,k6,k7,k8
      real*8, dimension(natm*3) :: xcoo
      real*8                    :: intcoo,pi,dot,len1,len2,len3,len4
      real*8, dimension(3)      :: vec1,vec2,vec3,vec4,cross1,cross2
      
      pi=4.0d0*datan(1.0d0)

      k1=iatm(1)
      k2=iatm(2)
      k3=iatm(3)
      k4=iatm(4)
      k5=iatm(5)
      k6=iatm(6)
      k7=iatm(7)
      k8=iatm(8)
      
      len1=0.0d0
      len2=0.0d0
      len3=0.0d0
      len4=0.0d0
      
      do i=1,3
         vec1(i)=xcoo(k2*3-3+i)-xcoo(k1*3-3+i)
         vec2(i)=xcoo(k4*3-3+i)-xcoo(k3*3-3+i)
         vec3(i)=xcoo(k6*3-3+i)-xcoo(k5*3-3+i)
         vec4(i)=xcoo(k8*3-3+i)-xcoo(k7*3-3+i)
         len1=len1+vec1(i)**2
         len2=len2+vec2(i)**2
         len3=len3+vec3(i)**2
         len4=len4+vec4(i)**2
      enddo
      len1=sqrt(len1)
      len2=sqrt(len2)
      len3=sqrt(len3)
      len4=sqrt(len4)
      vec1=vec1/len1
      vec2=vec2/len2
      vec3=vec3/len3
      vec4=vec4/len4

      cross1=cross_product(vec1,vec2)
      cross2=cross_product(vec3,vec4)

      dot=dot_product(cross1,cross2)
      intcoo=acos(dot)
      intcoo=intcoo*180.0d0/pi

      return

    end subroutine calc_twist

!#######################################################################

    subroutine calc_cart(pcoo,xcoo,intmom)

      use expec
      use sysdef

      implicit none

      integer*8                 :: i,j,k
      real*8, dimension(natm*3) :: pcoo,xcoo,tvec
      real*8                    :: intmom
      real*8, dimension(3,3)    :: rotmat
      
!-----------------------------------------------------------------------
! Exit if we are in the position representation
!-----------------------------------------------------------------------
      if (.not.lmomrep) then
         write(6,'(/,2x,a,/)') 'Reduced densities for Cartesian &
              displacements are currently only available in the &
              momentum representation.'
         STOP
      endif

!-----------------------------------------------------------------------
! Rotate the vector to the same frame as the current geometry
!-----------------------------------------------------------------------
      ! (1) Calculate the transformation matrix rotmat
      call minrmsd(xcoo,refgeom,rotmat)

      ! (2) Rotate the Cartesian vector cartvec
      do i=1,natm
         do j=1,3
            tvec(i*3-3+j)=0.0d0
            do k=1,3
               tvec(i*3-3+j)=tvec(i*3-3+j)+&
                    rotmat(k,j)*cartvec(i*3-3+k)
            enddo
         enddo
      enddo

!-----------------------------------------------------------------------
! Calculate the projection of the momentum onto the vector
!-----------------------------------------------------------------------
      intmom=dot_product(tvec,pcoo)

!      intmom=abs(intmom)

      return

    end subroutine calc_cart

!#######################################################################

    function cross_product(vec1,vec2) result(cp)

      implicit none

      real*8, dimension(3) :: vec1,vec2
      real*8, dimension(3) :: cp

      cp(1)=vec1(2)*vec2(3)-vec1(3)*vec2(2)
      cp(2)=vec1(3)*vec2(1)-vec1(1)*vec2(3)
      cp(3)=vec1(1)*vec2(2)-vec1(2)*vec2(1)

      return

    end function cross_product

!#######################################################################
! minrmsd: puts the Cartesian coordinates x2 into maximum coincidence
!          with the Cartesian coordinates x1 using the Kabsch algorithm
!#######################################################################

    subroutine minrmsd(x1,x2,rotmat)

      use sysdef

      implicit none

      integer*8                        :: i,j,k,l
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

      det=finddet2(covar,3)

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

      return

    end subroutine minrmsd

!#######################################################################

    function finddet2(matrix,n)

      implicit none

      integer*4              :: n,i,j,k,l
      real*8, dimension(n,n) :: matrix
      real*8                 :: m,temp,finddet2
      logical(kind=4)        :: detexists = .TRUE.      

      l=1

      !Convert to upper triangular form
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
               finddet2=0
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
      
      !Calculate determinant by finding product of diagonal elements
      finddet2=l
      do i=1,n
         finddet2=finddet2*matrix(i,i)
      enddo
   
    end function finddet2

!#######################################################################

  end module intcoo

  module permutemod

    contains

!#######################################################################
      
      recursive subroutine permutate(E, P) 

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

!You can test it with something like: 
!  read *, n 
!  allocate( P(product((/(i, i=1,n)/)),  n) ) 
!  call permutate( (/(i, i=1,n)/),  P ) 
!  do j=1,size(P,1) 
!    print *, P(j,:) 
!  end do

!#######################################################################

    end module permutemod
    

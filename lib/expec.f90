  module expec

    implicit none

    save

!-----------------------------------------------------------------------
! ijob:       integer specifying the job type:
!
!             ijob=1 <-> adiabatic populations
!                  2 <-> readind/writing of spawned geometries
!                  3 <-> reduced densities in terms of internal 
!                        coordinates
!                  4 <-> preparation of dyson orbital calculation input
!                  5 <-> calculation of RMSDs between spawn geometries
!                        and given conical intersection geometries
!
! adpop(s,n): adiabatic population of the sth state at the nth timestep
!
! ioutgeom:   specifies the format to be used in writing the spawn 
!             geometries:
!
!             0 <-> xyz, all geometries written to a single file
!             1 <-> Columbus geom format, single geometry per file
!
! iatm:       array containing integers that index the atoms involved
!             in the definition of an internal coordinate of interest
!
! ityp:       internal coordinate type
!
! vecfile:    name of the xyz file containing a reference geometry and
!             Cartesian vector. Used in the definition of the internal
!             coordinate if the 'cartvec' option is used
!
! cartvec:    Cartesian vector
!
! refgeom:    reference geometry that veccoo is defined wrt
!
! nmc:        no. samples to be used in the Monte Carlo procedure for
!             the calculation of reduced densities
!
! dstep:      time step for the calculation of reduced densities is
!             given by dstep*dt, where dt is the FMS timestep
!
! dgrid:      definition of the grid used in the discretistion of the
!             internal coordinate range used in the Monte Carlo
!             procedure for the calculation of reduced densities:
!             
!             dgrid(1): lower bound of the internal coordinate interval
!             dgrid(2): lower bound of the internal coordinate interval
!             dgrid(3): length of the partitions of the internal 
!                       coordinate interval
!             dgrid(4): no. partitions of the internal coordinate
!                       interval
!
! pfunc:      holds the values of the integrals of the reduced density
!             over each partion of the internal coordinate interval at
!             each timestep
!
! lrenorm:    logical flag to determine whether or not the coefficient
!             vector for each IFG is artificially renormalised at each
!             timestep
!
! nionize:    the no. of ionization channels to be considered in a
!             TRPES calculation
!
! iionize:    array holding the indices of the pairs neutral and
!             cationic states to be considered in a TRPES calculation
!
! ipshift:    constants by which to shift ionization potentials
!
! eprobe:     energy of the probe pulse for use in a TRPES calculation
!
! fwhm_e,fwhm_t: FHWMs of the Gaussian function used in the broadening
!                of the TRPES
!
! ngauss:     number of two-dimensional Gaussians to be admitted into
!             the construction of the TRPES
!
! gausspar:   parameters entering into the two-dimensional Gaussians
!             used in the construction of the TRPES:
!
!             gausspar(n,1): pre-multiplicative coefficient of the nth
!                            Gaussian
!             gausspar(n,2): width parameter of the nth Gaussian in
!                            the energy domain
!             gausspar(n,3): width parameter of the nth Gaussian in
!                            the time domain
!             gausspar(n,4): centre wrt E of the nth Gaussian
!             gausspar(n,5): centre wrt t of the nth Gaussian
!
! egrid:      upper and lower bounds for the photoelectron energy in
!             the calculation of the TRPES, as well as the no. points
!
! tgrid:      upper and lower bounds for the pump-probe delay in the
!             calculation of the TRPES, as well as the no. points
!-----------------------------------------------------------------------
    integer*8                              :: ijob,ioutgeom,nionize,&
                                              ngauss
    integer*8, dimension(10)               :: iatm
    integer*8                              :: ityp,nmc,dstep,dstate
    integer*8, dimension(:,:), allocatable :: iionize
    real*8, dimension(:), allocatable      :: cartvec,refgeom
    real*8, dimension(4)                   :: dgrid
    real*8, dimension(:,:), allocatable    :: adpop,gausspar
    real*8, dimension(:), allocatable      :: pfunc,ipshift
    real*8                                 :: eprobe,fwhm_e,fwhm_t
    real*8, dimension(3)                   :: egrid,tgrid
    logical(kind=4)                        :: lrenorm
    character(len=80)                      :: cifile,dnormfile,&
                                              vecfile

    ! Seam distance projection
    character(len=80)                      :: hfile
    character(len=80), dimension(2)        :: gfile
    real*8, dimension(:), allocatable      :: cigeom
    real*8, dimension(:,:), allocatable    :: branchvec,branchproj
    
    
    ! adcprep
    integer*8, dimension(3)                :: nlines
    real*8                                 :: thrsh_alive,tfinal
    character(len=80), dimension(3)        :: adcfile
    character(len=120), dimension(3,200)   :: adcinp
    logical(kind=4)                        :: ldummy

    ! adc_trxas
    real*8                                 :: gamma
    integer*8                              :: siord
    character(len=120)                     :: adcdir_file
   
  end module expec

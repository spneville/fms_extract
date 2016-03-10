  program fms_extract
    
    use expec, only: ityp,vecfile
    use sysdef, only: lmomrep

    implicit none
    
    character(len=80)                            :: adirfile
    character(len=80), dimension(:), allocatable :: adir

!-----------------------------------------------------------------------
! Read input
!-----------------------------------------------------------------------
    call rdinpfile(adirfile)
    
!-----------------------------------------------------------------------
! Read list of directories
!-----------------------------------------------------------------------
    call rddir(adirfile,adir)

!-----------------------------------------------------------------------
! Read atom labels and masses
!-----------------------------------------------------------------------
    call rdgeom(adir)

!-----------------------------------------------------------------------
! Assign Gaussian widths
!-----------------------------------------------------------------------
    call getwidths

!-----------------------------------------------------------------------
! Read the wavefunction and spawn geometries for each bundle of
! trajectories
!-----------------------------------------------------------------------
    call rdtraj(adir)

!-----------------------------------------------------------------------
! Determine whether or not a dummy atom is present (we assume that
! this is placed at the centre of mass)
!
! N.B., this only works if the INPUT directory is present, so for now
!       we will leave this as a user-specified argument
!-----------------------------------------------------------------------
!    call getdummy(adir(1))

!-----------------------------------------------------------------------
! Calculate the requested expectation values
!-----------------------------------------------------------------------
    call calc_expec

!-----------------------------------------------------------------------
! Output data
!-----------------------------------------------------------------------
    call wrout

    STOP
    
  contains

!#######################################################################

    subroutine rdinpfile(adirfile)

      use iomod
      use parsemod
      use errormod
      use expec
      use dysonmod

      implicit none

      integer                              :: unit,i,j,k,ndef,n,s,s1
      real*8, dimension(100)               :: tmpshift
      real*8                               :: ftmp,crosscorr,dtprobe
      character(len=80)                    :: adirfile,ainp
      character(len=60)                    :: msg

      integer                              :: inkw
      integer,   parameter                 :: maxkw=60
      integer,   dimension(maxkw)          :: ilkw
      character(len=120), dimension(maxkw) :: keyword

!-----------------------------------------------------------------------
! Set defaults
!-----------------------------------------------------------------------
      ijob=0
      iatm=0
      ityp=0
      ioutgeom=0
      nmc=20000
      dstep=20
      dstate=0
      vecfile=''
      ainp=''
      adirfile=''
      lrenorm=.false.
      lmomrep=.false.
      
      acol_n=''
      acol_c=''
      aprep_n=''
      aprep_c=''

      cifile=''

      dnormfile=''
      nionize=0
      tmpshift=0.0d0
      eprobe=0.0d0
      fwhm_e=0.0d0
      fwhm_t=0.0d0
      crosscorr=0.0d0
      dtprobe=0.0d0
      egrid(1)=-999
      egrid(2)=-999
      egrid(3)=100
      tgrid(1)=-999
      tgrid(2)=-999
      tgrid(3)=100

      adcfile=''
      
      ldummy=.false.

      thrsh_alive=0.0d0
      
      tfinal=0.0d0

      adcdir_file=''

      gamma=0.0d0
      
      siord=0

      hfile=''
      gfile=''
      
      npermute=0

      lcifilter=.false.
      cifstate=0

!-----------------------------------------------------------------------
! Read input file name
!-----------------------------------------------------------------------
      call getarg(1,ainp)
      if (ainp.eq.'') then
         write(6,'(/,2x,a,/)') 'Input file not given'
         STOP
      endif

!-----------------------------------------------------------------------
! Open input file
!-----------------------------------------------------------------------
      unit=20
      open(unit,file=ainp,form='formatted',status='old')

!-----------------------------------------------------------------------
! Read input file
!-----------------------------------------------------------------------
5     continue
      call rdinp(unit,keyword,inkw,ilkw)

      i=0
      if (keyword(1).ne.'end-input') then
10       continue
         i=i+1

         if (keyword(i).eq.'dir_file') then            
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),'(a)') adirfile
            else
               goto 100
            endif
            
         else if (keyword(i).eq.'job') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               if (keyword(i).eq.'adpop') then
                  ijob=1
               else if (keyword(i).eq.'spawngeom') then
                  ijob=2
               else if (keyword(i).eq.'density') then
                  ijob=3
               else if (keyword(i).eq.'dyson_prep') then
                  ijob=4
               else if (keyword(i).eq.'ci_rmsd') then
                  ijob=5
               else if (keyword(i).eq.'dyson_trpes') then
                  ijob=6
               else if (keyword(i).eq.'adc_prep') then
                  ijob=7
               else if (keyword(i).eq.'adc_trxas') then
                  if (keyword(i+1).eq.',') then
                     i=i+2
                     if (keyword(i).eq.'bound') then
                        ijob=8
                     else if (keyword(i).eq.'cont') then
                        ijob=9
                     else
                        msg='Unknown keyword: '//trim(keyword(i))
                        call errcntrl(msg)
                     endif
                  else
                     msg='The portion of the TR-XAS spectrum (bound or &
                          cont has not been specified)'
                     call errcntrl(msg)
                  endif
               else if (keyword(i).eq.'density_2d') then
                  ijob=10
               else
                  msg='Unknown job type: '//trim(keyword(i))
                  call errcntrl(msg)
               endif
            else
               goto 100
            endif

         else if (keyword(i).eq.'internal') then
            if (keyword(i+1).eq.'=') then
               i=i+2               
               if (keyword(i).eq.'length') then
                  ityp=1
                  ndef=2
               else if (keyword(i).eq.'angle') then
                  ityp=2
                  ndef=3
               else if (keyword(i).eq.'dihedral') then
                  ityp=3
                  ndef=4
               else if (keyword(i).eq.'twist') then
                  ityp=4
                  ndef=8
               else if (keyword(i).eq.'pyr') then
                  ityp=5
                  ndef=4
               else if (keyword(i).eq.'pyr2') then
                  ityp=6
                  ndef=6
               else if (keyword(i).eq.'cartvec') then
                  ityp=-1
               else if (keyword(i).eq.'seam') then
                  ityp=-2
               else
                  msg='Unknown internal coordinate type: '//trim(keyword(i))
                  call errcntrl(msg)
               endif
            else
               goto 100
            endif
            ! For an internal curvilinear coordinate, read indices of
            ! the atoms entering into the definition of the given
            ! internal coordinate
            if (ityp.gt.0) then
               if (keyword(i+1).ne.',') then
                  msg='Atom numbers for internal coordinate type '&
                       &//trim(keyword(i))//' not given'
                  call errcntrl(msg)
               endif               
               do j=1,ndef
                  i=i+2
                  read(keyword(i),*) iatm(j)
               enddo
            endif
            
            else if (keyword(i).eq.'internal2') then
            if (keyword(i+1).eq.'=') then
               i=i+2               
               if (keyword(i).eq.'length') then
                  ityp2=1
                  ndef=2
               else if (keyword(i).eq.'angle') then
                  ityp2=2
                  ndef=3
               else if (keyword(i).eq.'dihedral') then
                  ityp2=3
                  ndef=4
               else if (keyword(i).eq.'twist') then
                  ityp2=4
                  ndef=8
               else if (keyword(i).eq.'pyr') then
                  ityp2=5
                  ndef=4
               else if (keyword(i).eq.'pyr2') then
                  ityp2=6
                  ndef=6
               else if (keyword(i).eq.'cartvec') then
                  ityp2=-1
               else if (keyword(i).eq.'seam') then
                  ityp2=-2
               else
                  msg='Unknown internal coordinate type: '//trim(keyword(i))
                  call errcntrl(msg)
               endif
            else
               goto 100
            endif
            ! For an internal curvilinear coordinate, read indices of
            ! the atoms entering into the definition of the given
            ! internal coordinate
            if (ityp2.gt.0) then
               if (keyword(i+1).ne.',') then
                  msg='Atom numbers for internal coordinate type '&
                       &//trim(keyword(i))//' not given'
                  call errcntrl(msg)
               endif               
               do j=1,ndef
                  i=i+2
                  read(keyword(i),*) iatm2(j)
               enddo
            endif

            ! For a rectilinear Cartesian vector, read the name of the
            ! xyz file containing the reference geometry and the
            ! vector
            if (ityp.lt.0) then
               i=i+2
               read(keyword(i),'(a)') vecfile
            endif

         else if (keyword(i).eq.'nmc') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) nmc
            else
               goto 100
            endif

         else if (keyword(i).eq.'dgrid') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) dgrid(1)
               if (keyword(i+1).eq.',') then
                  i=i+2
                  read(keyword(i),*) dgrid(2)
               else
                  msg='Internal coordinate upper bound has not been given'
                  call errcntrl(msg)
               endif
               if (keyword(i+1).eq.',') then
                  i=i+2
                  read(keyword(i),*) dgrid(4)
               else
                  msg='Number of internal coordinate partitions has not '&
                       //'been given'
                  call errcntrl(msg)
               endif               
               dgrid(3)=(dgrid(2)-dgrid(1))/dgrid(4)
            else
               goto 100
            endif

         else if (keyword(i).eq.'dgrid2') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) dgrid2(1)
               if (keyword(i+1).eq.',') then
                  i=i+2
                  read(keyword(i),*) dgrid2(2)
               else
                  msg='Internal coordinate upper bound has not been given'
                  call errcntrl(msg)
               endif
               if (keyword(i+1).eq.',') then
                  i=i+2
                  read(keyword(i),*) dgrid2(4)
               else
                  msg='Number of internal coordinate partitions has not '&
                       //'been given'
                  call errcntrl(msg)
               endif               
               dgrid2(3)=(dgrid2(2)-dgrid2(1))/dgrid2(4)
            else
               goto 100
            endif

         else if (keyword(i).eq.'dstep') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) dstep
            else
               goto 100
            endif

         else if (keyword(i).eq.'dstate') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) dstate
            else
               goto 100
            endif

         else if (keyword(i).eq.'renorm') then
            lrenorm=.true.
            
         else if (keyword(i).eq.'geom_format') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               if (keyword(i).eq.'columbus') then
                  ioutgeom=1
               else
                  msg='Geometry format: '//trim(keyword(i))//' not recognised'
                  call errcntrl(msg)
               endif
            else
               goto 100
            endif

         else if (keyword(i).eq.'col_dir_cation') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               acol_c=keyword(i)
            else
               goto 100
            endif

         else if (keyword(i).eq.'col_dir_neutral'.or.keyword(i).eq.'col_dir') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               acol_n=keyword(i)
            else
               goto 100
            endif

         else if (keyword(i).eq.'prep_dir_cation') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               aprep_c=keyword(i)
            else
               goto 100
            endif

         else if (keyword(i).eq.'prep_dir_neutral') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               aprep_n=keyword(i)
            else
               goto 100
            endif

         else if (keyword(i).eq.'ci_file') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               cifile=keyword(i)
            else
               goto 100
            endif

         else if (keyword(i).eq.'dyson_norm_file') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               dnormfile=keyword(i)
            else
               goto 100
            endif

         else if (keyword(i).eq.'trpes_states') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               ! First determine the no. ionization channels
               ! to be considered and allocate the iionize array
               nionize=0
               k=i-1
20             continue
               k=k+1
               if (keyword(k).eq.'(') then
                  nionize=nionize+1
                  k=k+3
               endif
               if (k.lt.inkw) goto 20
               allocate (iionize(nionize,2))
               ! Second, read the state indices
               k=i-1
               n=0
25             continue
               k=k+1
               if (keyword(k).eq.'(') then
                  n=n+1
                  k=k+1
                  read(keyword(k),*) iionize(n,1)
                  k=k+2
                  read(keyword(k),*) iionize(n,2)
               endif
               if (k.lt.inkw) goto 25               
               ! Set i=inkw so that we can move on to the next line
               i=inkw
            else
               goto 100
            endif            

         else if (keyword(i).eq.'ipshift_section') then
            tmpshift=0.0d0
30          continue
            call rdinp(unit,keyword,inkw,ilkw)
            if (keyword(1).ne.'end-ipshift_section') then
               read(keyword(1),*) s
               read(keyword(2),*) s1
               read(keyword(3),*) ftmp
               do k=1,nionize
                  if (iionize(k,1).eq.s.and.iionize(k,2).eq.s1) then
                     tmpshift(k)=ftmp
                  endif
               enddo
               goto 30
            endif
            
         else if (keyword(i).eq.'eprobe') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) eprobe
            else
               goto 100
            endif

         else if (keyword(i).eq.'cross_corr') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) crosscorr
            else
               goto 100
            endif

         else if (keyword(i).eq.'dtprobe') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) dtprobe
            else
               goto 100
            endif

         else if (keyword(i).eq.'e_bounds') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) egrid(1)
               if (keyword(i+1).eq.',') then
                  i=i+2
                  read(keyword(i),*) egrid(2)
               else
                  msg='Upper energy bound has not been given'
                  call errcntrl(msg)
               endif
               ! Optional: no. grid points
               if (keyword(i+1).eq.',') then
                  i=i+2
                  read(keyword(i),*) egrid(3)
               endif
            else
               goto 100
            endif

         else if (keyword(i).eq.'t_bounds') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) tgrid(1)
               if (keyword(i+1).eq.',') then
                  i=i+2
                  read(keyword(i),*) tgrid(2)
               else
                  msg='Upper time delay bound has not been given'
                  call errcntrl(msg)
               endif
               ! Optional: no. grid points
               if (keyword(i+1).eq.',') then
                  i=i+2
                  read(keyword(i),*) tgrid(3)
               endif
            else
               goto 100
            endif

         else if (keyword(i).eq.'momentum_rep') then
            lmomrep=.true.

         else if (keyword(i).eq.'adcfile_ip') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),'(a)') adcfile(1)
            else
               goto 100
            endif

         else if (keyword(i).eq.'adcfile_dv') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),'(a)') adcfile(2)
            else
               goto 100
            endif

         else if (keyword(i).eq.'adcfile_si') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),'(a)') adcfile(3)
            else
               goto 100
            endif
            
         else if (keyword(i).eq.'col_dummy') then
            ldummy=.true.

         else if (keyword(i).eq.'thrsh_alive') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) thrsh_alive
               thrsh_alive=thrsh_alive*41.341375d0
            else
               goto 100
            endif
            
         else if (keyword(i).eq.'tfinal') then
             if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) tfinal
               tfinal=tfinal*41.341375d0
            else
               goto 100
            endif

         else if (keyword(i).eq.'adc_dir_file') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),'(a)') adcdir_file
            else
               goto 100
            endif

         else if (keyword(i).eq.'gamma') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) gamma
               gamma=gamma/27.2113845d0
            else
               goto 100
            endif

         else if (keyword(i).eq.'si_order') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) siord
            else
               goto 100
            endif

         else if (keyword(i).eq.'hfile') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),'(a)') hfile
            else
               goto 100
            endif

         else if (keyword(i).eq.'gfile') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),'(a)') gfile(1)
               if (keyword(i+1).eq.',') then
                  i=i+2
                  read(keyword(i),'(a)') gfile(2)
               else
                  msg='Only one gradient file has been given with &
                       the gfile keyword'
                  call errcntrl(msg)
               endif               
            else
               goto 100
            endif

         else if (keyword(i).eq.'permute') then
             if (keyword(i+1).eq.'=') then
                i=i+2
                k=i
35              continue
                npermute=npermute+1
                if (keyword(k+1).eq.',') then
                   k=k+2
                   goto 35
                endif
                allocate(pindx(npermute))
                do k=1,npermute
                   read(keyword(i),*) pindx(k)
                   i=i+2
                enddo                
             else
                goto 100
             endif

          else if (keyword(i).eq.'ci_filter') then
             lcifilter=.true.
             ityp=-2
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) cifdthrsh
               ! Optional state index
               if (keyword(i+1).eq.',') then
                  i=i+2
                  read(keyword(i),*) cifstate
               endif
            else
               goto 100
            endif
                
         else
            ! Exit if the keyword is not recognised
            msg='Unknown keyword: '//trim(keyword(i))
            call errcntrl(msg)            
         endif

         ! If there are more keywords to be read on the current line,
         ! then read them, else read the next line
         if (i.lt.inkw) then
            goto 10
         else
            goto 5
         endif

         ! Exit if a required argument has not been given with a keyword
100      continue
         msg='No argument given with the keyword '//trim(keyword(i))
         call errcntrl(msg)
      endif

!-----------------------------------------------------------------------
! Close input file
!-----------------------------------------------------------------------
      close(unit)

!-----------------------------------------------------------------------
! Check that all required information has been specified
!-----------------------------------------------------------------------
      if (ijob.eq.0) then
         msg='Job type not specified'
         call errcntrl(msg)
      endif
      
      if (adirfile.eq.'') then
         msg='Directory file not specified'
         call errcntrl(msg)
      endif

      if (ijob.eq.3) then
         if (ityp.eq.0) then
            msg='Internal coordinate definition not given'
            call errcntrl(msg)
         endif
         if (ityp.eq.-2) then
            if (hfile.eq.'') then
               msg='The NACT vector file has not been given'
               call errcntrl(msg)
            endif
            if (gfile(1).eq.''.or.gfile(2).eq.'') then
               msg='The gradient vector files have not been given'
               call errcntrl(msg)
            endif
            if (cifile.eq.'') then
               msg='The CI file has not been given'
               call errcntrl(msg)
            endif
         endif
      endif
      
      msg=''
      if (ijob.eq.4) then
         if (acol_n.eq.'') msg='Neutral Columbus directory name not given'
         if (acol_c.eq.'') msg='Cation Columbus directory name not given'
         if (aprep_n.eq.'') msg='Neutral superdyson prep directory name not given'
         if (aprep_c.eq.'') msg='Cation superdyson prep directory name not given'
         if (msg.ne.'') call errcntrl(msg)
      endif

      if (ijob.eq.5.and.cifile.eq.'') then
         msg='Conical intersection file not specified'
         call errcntrl(msg)
      endif

      if (ijob.eq.6) then
         if (nionize.eq.0) then
            msg='Ionization channels have not been specified'
            call errcntrl(msg)
         endif
         if (dnormfile.eq.'') then
            msg='Dyson norm file name not given'
            call errcntrl(msg)
         endif
         if (eprobe.eq.0) then
            msg='Probe energy not given'
            call errcntrl(msg)
         endif
         if (dtprobe.eq.0) then
            msg='Probe FWHM wrt t not given'
            call errcntrl(msg)
         endif
         if (crosscorr.eq.0) then
            msg='Pump/probe cross correlation not given'
            call errcntrl(msg)
         endif
         if (egrid(1).eq.-999) then
            msg='Bounds on the photoelectron energy have not been given'
            call errcntrl(msg)
         endif
         if (tgrid(1).eq.-999) then
            msg='Bounds on the pump-probe delay have not been given'
            call errcntrl(msg)
         endif
      endif

      if (ijob.eq.7) then
         if (adcfile(1).eq.''.and.adcfile(2).eq.''.and.adcfile(3).eq.'') then
            msg='No ADC template files have been given'
            call errcntrl(msg)
         endif
      endif
      
      if (ijob.eq.8.or.ijob.eq.9) then
         if (adcdir_file.eq.'') then
            msg='The name of the ADC directory file has not been given'
            call errcntrl(msg)
         endif
         if (egrid(1).eq.-999) then
            msg='Bounds on the photoelectron energy have not been given'
            call errcntrl(msg)
         endif
         if (tgrid(1).eq.-999) then
            msg='Bounds on the pump-probe delay have not been given'
            call errcntrl(msg)
         endif
         if (crosscorr.eq.0.0d0) then
            msg='Pump/probe cross correlation not given'
            call errcntrl(msg)
         endif
         if (gamma.eq.0.0d0) then
            msg='The Gamma value has not been given'
            call errcntrl(msg)
         endif
      endif

      if (ijob.eq.9) then
         if (siord.eq.0) then
            msg='The Stieltjes imaging order has not been given'
            call errcntrl(msg)
         endif
      endif

!-----------------------------------------------------------------------
! If the job type is the calculation of a TRPES, then:
!
! (1) allocate and fill in the ipshift array
! (2) calculate the Gaussian width parameters for use in energetic and
!     temporal broadening of the spectrum
!-----------------------------------------------------------------------
      if (ijob.eq.6) then
         ! (1) ipshift array
         allocate(ipshift(nionize))
         ipshift=0.0d0
         do i=1,nionize
            ipshift(i)=tmpshift(i)
         enddo
         ! (2) Gaussian broadening parameters
         fwhm_t=crosscorr
         fwhm_e=en_fwhm(dtprobe)
      endif

      ! TR-TXAS
      if (ijob.eq.8.or.ijob.eq.9) then
         fwhm_t=crosscorr
      endif

      return

    end subroutine rdinpfile

!#######################################################################

    function en_fwhm(deltat)

      implicit none

      real*8 :: en_fwhm,deltat

      ! Calculate the FWHM in the frequency domain using units of eV
      ! (from units of fs)
      en_fwhm=1.823832d0/deltat

      return
      
    end function en_fwhm

!#######################################################################

    subroutine rddir(adirfile,adir)
      
      use sysdef
      use parsemod
      use errormod

      implicit none

      integer                                      :: unit,n
      character(len=80)                            :: adirfile,string
      character(len=80), dimension(:), allocatable :: adir

      integer                                      :: i,k
      integer                                      :: inkw
      integer, parameter                           :: maxkw=60
      integer, dimension(maxkw)                    :: ilkw
      character(len=120), dimension(maxkw)         :: keyword
      character(len=60)                            :: msg

!-----------------------------------------------------------------------
! Open directories file
!-----------------------------------------------------------------------
      unit=20
      open(unit,file=adirfile,form='formatted',status='old')

!-----------------------------------------------------------------------
! First pass: determine the no. directories and allocate the adir array
!-----------------------------------------------------------------------      
      nintraj=0
5     continue
      call rdinp(unit,keyword,inkw,ilkw)
      if (keyword(1).ne.'end-file') then
         nintraj=nintraj+1
         goto 5
      endif

      allocate(adir(nintraj))
      adir=''

!-----------------------------------------------------------------------
! Second pass: read the directory names
!-----------------------------------------------------------------------
      rewind(unit)

      n=0
10    continue
      call rdinp(unit,keyword,inkw,ilkw)
      if (keyword(1).ne.'end-file') then
         n=n+1
         read(keyword(1),'(a)') adir(n)
         goto 10
      endif

!-----------------------------------------------------------------------
! Open directories file
!-----------------------------------------------------------------------
      close(unit)

      return

    end subroutine rddir

!#######################################################################

    subroutine rdgeom(adir)

      use sysdef

      implicit none
      
      integer                               :: unit,i,n
      character(len=80), dimension(nintraj) :: adir
      character(len=100)                    :: ageom

!-----------------------------------------------------------------------
! Open the Geometry.dat file
!-----------------------------------------------------------------------
      unit=20
      ageom=trim(adir(1))//'/Geometry.dat'
      open(unit,file=ageom,form='formatted',status='old')

!-----------------------------------------------------------------------
! Read the number of atoms and allocate arrays
!-----------------------------------------------------------------------
      read(unit,*)
      read(unit,*) n
      natm=n
      allocate(atlbl(n),atnum(n),atmass(n))

!-----------------------------------------------------------------------
! Read the atom labels and determine the atomic masses/numbers
!-----------------------------------------------------------------------
      do i=1,n
         read(unit,'(a2)') atlbl(i)
         if (atlbl(i).eq.'C') then
            atnum(i)=6.0d0
            atmass(i)=12.0d0
         else if (atlbl(i).eq.'H') then
            atnum(i)=1.0d0
            atmass(i)=1.00782504d0
         else if (atlbl(i).eq.'N') then
            atnum(i)=7.0d0
            atmass(i)=14.00307401d0
         else
            write(6,'(/,2(2x,a),/)') 'Unknown atom type:',&
                 trim(atlbl(i))
            STOP
         endif         
      enddo

!-----------------------------------------------------------------------
! Close the Geometry.dat file
!-----------------------------------------------------------------------
      close(unit)

      return

    end subroutine rdgeom

!#######################################################################

    subroutine rdtraj(adir)

      use trajdef
      use sysdef
      use expec, only: lrenorm

      implicit none

      integer                               :: i,k,nspawn
      character(len=80), dimension(nintraj) :: adir
      character(len=150)                    :: atmp

      write(6,'(/,70a)') ('-', i=1,70)
      write(6,'(21x,a)') 'Reading the FMS output files'
      write(6,'(70a)') ('-', i=1,70)

!-----------------------------------------------------------------------   
! Allocate the traj derived data type
!-----------------------------------------------------------------------   
      allocate(traj(nintraj))

!-----------------------------------------------------------------------   
! Loop over initial trajectories (i.e., directories), create temporary
! directories/files and read the trajectory information
!-----------------------------------------------------------------------
      do i=1,nintraj

         ! Output where we are up to
         k=len_trim(adir(i))
         write(6,'(a21,i3,1x,a3,i3,a2,1x,a)') 'Processing directory',&
              i,'of',nintraj,':',adir(i)(1:k)

         ! Make temporary directory to work from
         call mktmpdir(adir(i),atmp)

         ! Determine propagation length, dt, nspawn, etc
         nspawn=traj(i)%ntraj-1
         call propinfo(i,nspawn,atmp,adir(i))         

         ! Allocate and initialise arrays
         call alloc(i)

         ! Read the TrajDump files
         call rdtrajdump(i,atmp)

         ! Read and output the spawn geometries.
         ! We also read the parent/spawn trajectory information here.
         call getspawngeom(atmp,nspawn,i)

         ! Clean up
         call cleanup(adir(i))
         
      enddo

!-----------------------------------------------------------------------
! If requested, renormalise the coefficient vector for each IFG at
! each timestep
!-----------------------------------------------------------------------
      if (lrenorm) call renorm_cvec

      return

    end subroutine rdtraj

!#######################################################################

    subroutine getdummy(dir1)

      use sysdef
      use expec, only: ldummy
      use iomod
      
      implicit none

      integer            :: unit,i
      character(len=80)  :: dir1
      character(len=120) :: filename,string

      filename=trim(dir1)//'INPUT/geom'

      call getfreeunit(unit)

      open(unit,file=filename,form='formatted',status='old')

      ldummy=.false.
      do i=1,natm
         read(unit,'(a)') string
         if (index(string,'X').ne.0) ldummy=.true.
      enddo
         
      close(unit)

      return
      
    end subroutine getdummy
      
!#######################################################################

    subroutine mktmpdir(currdir,atmp)

      implicit none

      integer            :: k
      character(len=80)  :: currdir
      character(len=150) :: stem,atmp
      character(len=250) :: acmnd
      logical            :: ldir,found

!-----------------------------------------------------------------------
! Make temporary directory
!-----------------------------------------------------------------------
      atmp=''
      stem=''

      ! Set directory stem
      ! If the OUTPUT directory is not present, then we assume that
      ! all the output files are in the main directory
      inquire(file=trim(currdir)//'/OUTPUT/.',exist=found)
      if (found) then
         stem=trim(currdir)//'/OUTPUT'
      else 
         stem=trim(currdir)
      endif
      
      ! If TMP directory already exists, remove
      inquire(file=trim(stem)//'/TMP/.',exist=ldir)
      if (ldir) then
         acmnd='rm -rf '//trim(stem)//'/TMP/'
         call system(acmnd)
      endif
      
      ! Create TMP directory
      atmp=trim(stem)//'/TMP'
      acmnd='mkdir '//trim(stem)//'/TMP'
      call system(acmnd)

!-----------------------------------------------------------------------
! Copy all gzipped files to TMP and unzip them
!-----------------------------------------------------------------------
!      acmnd='cp '//trim(stem)//'/*.gz '//trim(stem)//'/TMP'
      acmnd='cp '//trim(stem)//'/*.* '//trim(stem)//'/TMP'      
      call system(acmnd)

      acmnd='gunzip -f '//trim(stem)//'/TMP/*.gz'
      call system(acmnd)

      return
      
    end subroutine mktmpdir

!#######################################################################

    subroutine propinfo(itraj,nspawn,atmp,currdir)

      use trajdef
      use sysdef

      implicit none

      integer                              :: itraj,nspawn,k,unit
      real*8                               :: ftmp
      character(len=80)                    :: currdir,string
      character(len=150)                   :: atmp,afile

!-----------------------------------------------------------------------
! Read the FMS log file and determine the no. spawned trajectories, tf
! and dt
!-----------------------------------------------------------------------
      ! (1) Spawn.log: no. spawned trajectories
      unit=20
      afile=trim(atmp)//'/Spawn.log'
      open(unit,file=afile,form='formatted',status='old')
      
      nspawn=0
10    continue
      read(unit,'(a)',end=11) string
      if (index(string,'#').eq.0) nspawn=nspawn+1
      goto 10
11    continue

      close(unit)

      ! (2) Control.dat: system and propagation info
      afile=''
      k=len_trim(currdir)
      afile=currdir(1:k)//'/Control.dat'
      open(unit,file=afile,form='formatted',status='old')

20    continue
      read(unit,'(a)',end=25) string
      call getstart(string,k)
      if (string(k:k+13).eq.'SimulationTime') then
         call getval(string,tf)
      else if (string(k:k+7).eq.'TimeStep'.and.&
           index(string,'Rejection').eq.0) then
         call getval(string,dt)
      else if (string(k:k+8).eq.'NumStates') then
         call getval(string,ftmp)
         nsta=int(ftmp)
      else if (string(k:k+11).eq.'NumParticles') then
         call getval(string,ftmp)
         natm=int(ftmp)
      endif
      goto 20

25    continue
      close(unit)

      traj(itraj)%ntraj=nspawn+1
      ftmp=tf/dt
      nstep=int(ftmp)+1

      return

    end subroutine propinfo

!#######################################################################

    subroutine alloc(itraj)

      use sysdef
      use trajdef
      use expec

      implicit none

      integer :: itraj,ntraj

      ntraj=traj(itraj)%ntraj

!-----------------------------------------------------------------------
! Trajectory parameters
!-----------------------------------------------------------------------
      ! Positions
      if (allocated(traj(itraj)%r)) deallocate(traj(itraj)%r)
      allocate(traj(itraj)%r(ntraj,nstep,natm*3))
      traj(itraj)%r=0.0d0

      ! Momenta
      if (allocated(traj(itraj)%p)) deallocate(traj(itraj)%p)
      allocate(traj(itraj)%p(ntraj,nstep,natm*3))
      traj(itraj)%p=0.0d0

      ! Phases
      if (allocated(traj(itraj)%phase)) deallocate(traj(itraj)%phase)
      allocate (traj(itraj)%phase(ntraj,nstep))
      traj(itraj)%phase=0.0d0

      ! Coefficients
      if (allocated(traj(itraj)%coe)) deallocate(traj(itraj)%coe)
      allocate(traj(itraj)%coe(ntraj,nstep))
      traj(itraj)%coe=0.0d0

      ! Timesteps at which trajectories were killed
      if (allocated(traj(itraj)%tkill)) deallocate(traj(itraj)%tkill)
      allocate(traj(itraj)%tkill(ntraj))
      traj(itraj)%tkill=0.0d0

      ! State indices for each trajectory
      if (allocated(traj(itraj)%ista)) deallocate(traj(itraj)%ista)
      allocate(traj(itraj)%ista(ntraj))
      traj(itraj)%ista=0

      ! ID's of the parent trajectories for each spawned trajectory
      if (allocated(traj(itraj)%ispawn)) deallocate(traj(itraj)%ispawn)
      allocate(traj(itraj)%ispawn(ntraj))

      ! Spawn times
      if (allocated(traj(itraj)%tspawn)) deallocate(traj(itraj)%tspawn)
      allocate(traj(itraj)%tspawn(ntraj))

!-----------------------------------------------------------------------
! Expectation values
!-----------------------------------------------------------------------
      if (allocated(adpop)) deallocate(adpop)
      allocate(adpop(nsta,nstep))
      adpop=0

      return

    end subroutine alloc

!#######################################################################

    subroutine getwidths

      use sysdef

      implicit none
      
      integer :: i
      real*8  :: width

!-----------------------------------------------------------------------
! Default Gaussian widths
!-----------------------------------------------------------------------
      allocate(alpha(natm*3))

      do i=1,natm

         select case(atlbl(i))

         case('H','h')
            width=4.5d0

         case('C','c')
            width=22.5d0

         case('N','n')
            width=19.5d0

         case('O','o')
            width=13.0d0

         case('S','s')
            width=17.5d0

         case('F','f')
            width=8.5d0

         case default
            write(6,'(/,2(2x,a),/)') 'Unknown atom type:',&
                 trim(atlbl(i))
            STOP

         end select

         alpha(i*3-2:i*3)=width

      enddo

      return
      
    end subroutine getwidths

!#######################################################################

    subroutine rdtrajdump(itraj,atmp)
      
      use sysdef
      use trajdef

      implicit none

      integer                     :: itraj,ntraj,unit,i,j,k,&
                                     n,ncoo
      real*8                      :: t,gtmp,crtmp,citmp,dum,stmp
      real*8, dimension(3*natm)   :: rtmp,ptmp
      character(len=150)          :: atmp,a2
      character(len=11)           :: a1
      character(len=50)           :: fmat

      ntraj=traj(itraj)%ntraj

!-----------------------------------------------------------------------
! Read the amplitudes, positions, momenta and phases for each
! trajectory
!-----------------------------------------------------------------------
      unit=20
      ncoo=3*natm
      call wrfmat(fmat,ncoo)

      ! Loop over trajectories
      do i=1,ntraj

         ! Open current TrajDump file
         a1=''
         if (i.lt.10) then
            write(a1,'(a9,i1)') 'TrajDump.',i
         else
            write(a1,'(a9,i2)') 'TrajDump.',i
         endif
         a2=trim(atmp)//'/'//trim(a1)
         open(unit,file=a2,form='formatted',status='old')

         ! Read current TrajDump file
         ! First read the spawn time
         read(unit,*)
         read(unit,'(F10.2)') t
         traj(itraj)%tspawn(i)=int(t/dt)+1
         backspace(unit)

10       continue
         read(unit,fmat,end=15) t,(rtmp(j),j=1,ncoo),(ptmp(j),j=1,ncoo),&
              gtmp,crtmp,citmp,dum,stmp
         n=1+(t/dt)
         ! Positons
         traj(itraj)%r(i,n,:)=rtmp(:)
         ! Momenta
         traj(itraj)%p(i,n,:)=ptmp(:)
         ! Phases
         traj(itraj)%phase(i,n)=gtmp
         ! Coefficients
         traj(itraj)%coe(i,n)=cmplx(crtmp,citmp)         
         goto 10

15       continue

         ! Close current TrajDump file
         close(unit)

         ! State index for the current trajectory
         traj(itraj)%ista(i)=int(stmp)
         
         ! Time step of death for the current trajectory
         traj(itraj)%tkill(i)=n

      enddo

      return

    end subroutine rdtrajdump

!#######################################################################

    subroutine wrfmat(fmat,ncoo)

      implicit none

      integer           :: ncoo
      character(len=50) :: fmat

      fmat=''
      if (ncoo.lt.10) then
         write(fmat,'(a7,i1,a6,i1,a13)') &
              '(F10.2,',ncoo,'F10.4,',ncoo,'F10.4,5F10.4)'
      else if (ncoo.lt.100) then
         write(fmat,'(a7,i2,a6,i2,a13)') &
              '(F10.2,',ncoo,'F10.4,',ncoo,'F10.4,5F10.4)'
      else
         write(fmat,'(a7,i3,a6,i3,a13)') &
              '(F10.4,',ncoo,'F10.4,',ncoo,'F10.4,5F10.2)'
      endif

      return

    end subroutine wrfmat

!#######################################################################

    subroutine getspawngeom(currdir,nspawn,ipass)

      use expec
      use sysdef
      use trajdef

      implicit none

      integer                           :: nspawn,unit,i,ipass,iout,k,m
      real*8, dimension(natm*3)         :: currxcoo
      character(len=2), dimension(natm) :: atmlbl
      character(len=8)                  :: atmp
      character(len=150)                :: currdir,aspawn,comment

!-----------------------------------------------------------------------   
! Open the output xyz file if this is the 1st pass
!-----------------------------------------------------------------------   
      iout=15
      if (ipass.eq.1) then
         open(iout,file='spawngeom.xyz',form='formatted',status='unknown')
      endif

      ! Loop over spawned trajectories
      unit=20
      do i=1,nspawn

         ! Open Spawn file
         atmp=''
         write(atmp(1:6),'(a)') 'Spawn.'
         if (i+1.lt.10) then
            write(atmp(7:7),'(i1)') i+1
         else
            write(atmp(7:8),'(i2)') i+1
         endif
         aspawn=''

         aspawn=trim(currdir)//'/'//trim(atmp)

         open(unit,file=aspawn,form='formatted',status='old')

         ! Read Spawn file
         read(unit,*)
         read(unit,'(a)') comment
         read(comment,'(37x,i3)') traj(ipass)%ispawn(i+1)

         do k=1,natm
            read(unit,*) atmlbl(k),(currxcoo(m), m=k*3-2,k*3)
         enddo

         ! Write to the spawngeom.xyz file
         write(iout,'(i3)') natm
         write(iout,'(a)') trim(comment)
         do k=1,natm
            write(iout,'(a2,3(3x,F13.9))') atmlbl(k),(currxcoo(m), m=k*3-2,k*3)
         enddo

         ! If the user has requested individual spawn geometry files,
         ! then write these
         if (ioutgeom.ne.0) call wr1spawngeom(ioutgeom,currxcoo,natm,&
              atmlbl,ipass,i)

         ! Close Spawn file
         close(unit)

      enddo

      ! If this is the last pass, then close the output xyz file
      if (ipass.eq.nintraj) close(iout)

      return

    end subroutine getspawngeom

!#######################################################################
    
    subroutine wr1spawngeom(ioutgeom,currxcoo,natm,atmlbl,ipass,igeom)

      implicit none

      integer                           :: ioutgeom,natm,ipass,igeom
      real*8, dimension(natm*3)         :: currxcoo
      character(len=2), dimension(natm) :: atmlbl
      character(len=80)                 :: acmnd
      logical(kind=4)                   :: ldir

!-----------------------------------------------------------------------
! If this is the 1st pass, then create the directory where the spawn
! geometries will be written
!-----------------------------------------------------------------------
      if (ipass.eq.1.and.igeom.eq.1) then
         inquire(file='spawn_geoms/.',exist=ldir)
         if (ldir) then
            acmnd='rm -rf spawn_geoms'
            call system(acmnd)
         endif
         acmnd='mkdir spawn_geoms'
         call system(acmnd)
      endif

!-----------------------------------------------------------------------
! Write the current spawn geometry to file
!-----------------------------------------------------------------------
      select case(ioutgeom)

      case(1) ! Columbus
         call wrspawngeom_columbus(currxcoo,natm,atmlbl,ipass,igeom)

      end select

      return

    end subroutine wr1spawngeom

!#######################################################################

    subroutine wrspawngeom_columbus(currxcoo,natm,atmlbl,ipass,igeom)

      implicit none

      integer                           :: natm,ipass,igeom,unit,i,j
      real*8, dimension(natm*3)         :: currxcoo
      real*8                            :: fmass,fatnum
      character(len=2), dimension(natm) :: atmlbl
      character(len=150)                :: aout,atmp1,atmp2

!-----------------------------------------------------------------------
! Write filename
!-----------------------------------------------------------------------      
      write(atmp1,'(i3)') ipass
      write(atmp2,'(i3)') igeom      
      aout='spawn_geoms/geom_'//trim(adjustl(atmp1))//'_'//trim(adjustl(atmp2))

!-----------------------------------------------------------------------
! Open file
!-----------------------------------------------------------------------
      unit=30
      open(unit,file=aout,form='formatted',status='new')

!-----------------------------------------------------------------------
! Write to file
!-----------------------------------------------------------------------
      do i=1,natm
         fmass=mass(atmlbl(i))
         fatnum=atnum(atmlbl(i))
         write(unit,'(1x,a1,4x,F4.1,3(2x,F12.8),3x,F11.8)') &
              atmlbl(i),fatnum,(currxcoo(j)/0.529177249d0,j=i*3-2,i*3),&
              fmass
      enddo

!-----------------------------------------------------------------------
! Close file
!-----------------------------------------------------------------------
      close(unit)

      return
      
    end subroutine wrspawngeom_columbus

!#######################################################################

    function mass(label)

      implicit none
      
      real*8           :: mass
      character(len=*) :: label

      if (label.eq.'H') then
         mass=1.00782504d0
      else if (label.eq.'C') then
         mass=12.00000000d0
      else
         write(6,'(2(2x,a))') 'Unknown atom type:',trim(label)
         STOP
      endif

      return

    end function mass

!#######################################################################

    function atnum(label)

      implicit none

      real*8           :: atnum
      character(len=*) :: label

      if (label.eq.'H') then
         atnum=1.0d0
      else if (label.eq.'C') then
         atnum=6.0d0
      else
         write(6,'(2(2x,a))') 'Unknown atom type:',trim(label)
         STOP
      endif

      return

    end function atnum

!#######################################################################

    subroutine getstart(string,k)

      implicit none

      integer           :: k,i
      character(len=80) :: string

      k=0
      do i=1,len_trim(string)
         if (string(i:i).ne.'') then
            if (k.eq.0) k=i
         endif
      enddo

      return

    end subroutine getstart

!#######################################################################
    
    subroutine getval(string,ftmp)

      implicit none
      
      integer           :: i,itmp,ilbl,jlbl
      real*8            :: ftmp
      character(len=80) :: string

      do i=1,len_trim(string)
         if (string(i:i).eq.'=') itmp=i
      enddo

      ilbl=0
      do i=itmp+1,len_trim(string)
         if (string(i:i).ne.' '.and.string(i:i).ne.achar(9)) then
            if (ilbl.eq.0) ilbl=i
         endif
      enddo

      jlbl=0
      do i=ilbl+1,len_trim(string)
         if (string(i:i).eq.' '.or.string(i:i).eq.achar(9)) then
            if (jlbl.eq.0) jlbl=i-1
         endif
      enddo

      read(string(ilbl:jlbl),*) ftmp

      return

    end subroutine getval

!#######################################################################

    subroutine cleanup(currdir)
      
      implicit none

      character(len=80)  :: currdir
      character(len=150) :: atmp
      logical            :: found

!-----------------------------------------------------------------------
! Remove the current temporary directory
!-----------------------------------------------------------------------
      inquire(file=trim(currdir)//'/OUTPUT/.',exist=found)
      
      if (found) then
         atmp='rm -r '//trim(currdir)//'/OUTPUT/TMP'
      else
         atmp='rm -r '//trim(currdir)//'/TMP'
      endif

      call system(atmp)

      return

    end subroutine cleanup

!#######################################################################

    subroutine renorm_cvec

      use trajdef
      use sysdef
      use gausstools, only: ispop

      implicit none

      integer         :: i,j,n
      real*8          :: sumsq
      complex*16      :: cj
      logical(kind=4) :: lalive

!-----------------------------------------------------------------------
! For each IFG at each timestep, renormalise the coefficient vector
!-----------------------------------------------------------------------
      do n=1,nstep ! Loop over timesteps
         do i=1,nintraj ! Loop over IFGs

            ! Skip dead IFGs
            lalive=ispop(i,n)
            if (.not.lalive) cycle

            sumsq=0.0d0
            do j=1,traj(i)%ntraj ! Loop over trajectories of the current IFG
               cj=traj(i)%coe(j,n)
               sumsq=sumsq+real(conjg(cj)*cj)               
            enddo
            traj(i)%coe(:,n)=traj(i)%coe(:,n)/sqrt(sumsq)
         enddo
      enddo

      return

    end subroutine renorm_cvec

!#######################################################################

    subroutine calc_expec

      use expec, only: ijob
      use density
      use density_2d
      use density_mom
      use adpop
      use dysonprep
      use cirmsd
      use trpes
      use adcprep
      use adctas

      implicit none

!-----------------------------------------------------------------------   
! Calculate the requested expectation values
!
! ijob=1 <-> Adiabatic populations
!      2 <-> Spawned geometries (do nothing here as this has already
!            been dealt with in rdtraj)
!      3 <-> One-dimensional reduced densities
!      4 <-> Preparation of input files for the calculation of Dyson
!            orbital norms
!      5 <-> Calculation of RMSDs between the spawn geometries and
!            given conical intersection geometries
!      6 <-> Calculation of TRPES using Dyson orbital norms
!      7 <-> Preparation of input files for ADC absorption or
!            photoionization cross-section calculations
!      8 <-> Calculation of the bound part of the TR-TXAS using 
!            ADC cross-sections
!      9 <-> Calculation of the continuum part of the TR-TXAS using 
!            ADC cross-sections
!     10 <-> Two-dimensional reduced densities
!-----------------------------------------------------------------------   
      if (ijob.eq.1) then
         call calcadpop
      else if (ijob.eq.3) then
         if (lmomrep) then
            call reddens_mom
         else
            call reddens
         endif
      else if (ijob.eq.4) then
         call mkdysoninp
      else if (ijob.eq.5) then
         call ci_spawn_rmsd
      else if (ijob.eq.6) then
         call trpes_dnorm
      else if (ijob.eq.7) then
         call mkadcinp
      else if (ijob.eq.8.or.ijob.eq.9) then
         call adc_trtxas
      else if (ijob.eq.10) then
         call reddens_2d
      endif

      return

    end subroutine calc_expec

!#######################################################################

    subroutine wrout
      
      use sysdef
      use expec, only: ijob

      implicit none

      if (ijob.eq.1) then
         call wradpop
      else if (ijob.eq.3) then
         call wrreddens
      endif

      return

    end subroutine wrout

!#######################################################################

    subroutine wradpop

      use sysdef
      use expec

      implicit none
      
      integer                       :: i,n,unit
      character(len=8)              :: aout

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

    subroutine wrreddens

      implicit none

      integer :: iout

!-----------------------------------------------------------------------
! Write gnuplot file
!-----------------------------------------------------------------------
      iout=30
      open(iout,file='dens.gnu',form='formatted',status='unknown')

      write(iout,'(a)') '# ~/.gnuplot'
      write(iout,'(a,/)') 'set palette @MATLAB'
      write(iout,'(a)') 'set pm3d map interpolate 0,0'
      write(iout,'(a)') 'set xlabel ''Time (fs)'''
      write(iout,'(a,/)') 'set ylabel ''x (a.u.)'''
      write(iout,'(a,/)') 'splot ''dens.dat'''
      write(iout,'(a)') 'pause -1'

      close(iout)

      return

    end subroutine wrreddens

!#######################################################################

  end program fms_extract

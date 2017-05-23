  program fms_extract
    
    use expec, only: ityp,vecfile,ijob,tfinal
    use sysdef, only: lmomrep,tf

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
! Read reference geometry Cartesian coordinates, atom labels and masses
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
! If it has not been set by the user, set tfinal to tf
!
! N.B., tf is the propagation time, whereas tfinal is the final time
!       to be used in the expectation value or spectrum calculations
!-----------------------------------------------------------------------
    if (tfinal.eq.0.0d0) tfinal=tf

!-----------------------------------------------------------------------
! Set up the centroids
!
! N.B., we should only be doing this if the job type requires
! centroid information
!-----------------------------------------------------------------------
    if (ijob.eq.15.or.ijob.eq.16) then
       call getcentroids
    endif

!-----------------------------------------------------------------------
! Calculate the requested expectation values
!-----------------------------------------------------------------------
    call calc_expec

!-----------------------------------------------------------------------
! Output data
!-----------------------------------------------------------------------
    call wrout

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
      integer,   parameter                 :: maxkw=200
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

      tdmmask=1

      gamma=0.0d0
      
      siord=0

      hfile=''
      gfile=''
      
      npermute=0

      lcifilter=.false.
      cifstate=0

      cirmsdsta=0

      ldiff=.false.
      osc0file=''

      ovrthrsh=0.75d0
      bastype=1
      
      gamessfile=''
      gmsdir_file=''

      ! Columbus post. prep. directory file
      colppdir_file=''
      
!-----------------------------------------------------------------------
! Read input file name
!-----------------------------------------------------------------------
      call getarg(1,ainp)
      if (ainp.eq.'') then
         write(6,'(/,2x,a,/)') 'Input file not given'
         STOP
      endif

      k=index(ainp,'.inp')
      if (k.eq.0) ainp=trim(ainp)//'.inp'

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
                     else if (keyword(i).eq.'trpes') then
                        ijob=11
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

               else if (keyword(i).eq.'ts-psg_prep' &
                    .or.keyword(i).eq.'tspsg_prep') then
                  ijob=12

               else if (keyword(i).eq.'gamess_prep') then
                  ijob=13
                  
               else if (keyword(i).eq.'gamess_trxas') then
                  ijob=14

               else if (keyword(i).eq.'dipole_prep') then
                  ijob=15
                  
               else if (keyword(i).eq.'dipole') then
                  ijob=16
                  
               else if (keyword(i).eq.'adc_orben') then
                  ijob=17

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
               else if (keyword(i).eq.'maxangle') then
                  ityp=7
                  ndef=6
               else if (keyword(i).eq.'minangle') then
                  ityp=8
                  ndef=6
               else if (keyword(i).eq.'anglediff') then
                  ityp=9
                  ndef=6
               else if (keyword(i).eq.'pyrdiff') then
                  ityp=10
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
               else if (keyword(i).eq.'maxangle') then
                  ityp2=7
                  ndef=6
               else if (keyword(i).eq.'minangle') then
                  ityp2=8
                  ndef=6
               else if (keyword(i).eq.'anglediff') then
                  ityp2=9
                  ndef=6
               else if (keyword(i).eq.'pyrdiff') then
                  ityp2=10
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
            
         else if (keyword(i).eq.'adcfile_dys') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),'(a)') adcfile(4)
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

         else if (keyword(i).eq.'tdm') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               tdmmask=0
31             continue
               if (keyword(i).eq.'x') then
                  tdmmask(1)=1
               else if (keyword(i).eq.'y') then
                  tdmmask(2)=1
               else if (keyword(i).eq.'z') then
                  tdmmask(3)=1
               else
                  msg='Unknown tdm component: '//trim(keyword(i))
                  call errcntrl(msg)
               endif
               if (keyword(i+1).eq.',') then
                  i=i+2
                  goto 31
               endif
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
            
         else if (keyword(i).eq.'cirmsd_states') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) cirmsdsta(1)
               if (keyword(i+1).eq.',') then
                  i=i+2
                  read(keyword(i),*) cirmsdsta(2)
               else
                  msg='Only one state index has been given with &
                       the cirmsd_states keyword'
                  call errcntrl(msg)
               endif
            else
               goto 100
            endif

         else if (keyword(i).eq.'tspsg_basis') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               if (keyword(i).eq.'gaussian') then
                  bastype=1
               else if (keyword(i).eq.'lowdin') then
                  bastype=2
               else
                  msg='Unknown TS-PSG basis type: '//trim(keyword(i))
                  call errcntrl(msg)
               endif
            else
               goto 100
            endif

         else if (keyword(i).eq.'ovrthrsh') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) ovrthrsh
            else
               goto 100
            endif

         else if (keyword(i).eq.'gamessfiles') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               ngamessinp=0
40             continue
               ngamessinp=ngamessinp+1
               read(keyword(i),'(a)') gamessfile(ngamessinp)
               if (keyword(i+1).eq.',') then
                  i=i+2
                  goto 40
               endif
            else
               goto 100
            endif

         else if (keyword(i).eq.'gamess_dir_file') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),'(a)') gmsdir_file
            else
               goto 100
            endif

         else if (keyword(i).eq.'colpp_dir_file') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),'(a)') colppdir_file
            else
               goto 100
            endif

         else if (keyword(i).eq.'diff') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               ldiff=.true.
               read(keyword(i),'(a)') osc0file
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
         if (adcfile(1).eq.''.and.adcfile(2).eq.''.and.adcfile(3).eq.''&
              .and.adcfile(4).eq.'') then
            msg='No ADC template files have been given'
            call errcntrl(msg)
         endif
      endif
      
      if (ijob.eq.8.or.ijob.eq.9.or.ijob.eq.11) then
         if (adcdir_file.eq.'') then
            msg='The name of the ADC directory file has not been given'
            call errcntrl(msg)
         endif
         if (egrid(1).eq.-999) then
            msg='The energy bounds have not been given'
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

      if (ijob.eq.14) then
         if (gmsdir_file.eq.'') then
            msg='The gamess directory file has not been given'
            call errcntrl(msg)
         endif
         if (egrid(1).eq.-999) then
            msg='The energy bounds have not been given'
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

      if (ijob.eq.16) then
         if (colppdir_file.eq.'') then
            msg='The name of the Columbus post-processing directory file &
                 has not been given'
            call errcntrl(msg)
         endif
      endif

      if (ijob.eq.17) then
         if (adcdir_file.eq.'') then
            msg='The name of the ADC directory file has not been given'
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

      ! ADC Dyson norm TRPES
      if (ijob.eq.11) then
         fwhm_t=crosscorr
         fwhm_e=en_fwhm(dtprobe)
      endif

!-----------------------------------------------------------------------
! If the job type is the calculation of a TRXAS, then:
!
! (1) Set the FWHM wrt time
!-----------------------------------------------------------------------
      if (ijob.eq.14.or.ijob.eq.8.or.ijob.eq.9) then
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
      integer, parameter                           :: maxkw=200
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
      
      integer                               :: unit,i,j,n
      character(len=80), dimension(nintraj) :: adir
      character(len=100)                    :: ageom
      character(len=2)                      :: atmp
      
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
      allocate(atlbl(n),atnum(n),atmass(n),r0(n*3))
      
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
! Read the reference Cartesian coordinates
!-----------------------------------------------------------------------
      rewind(unit)
      read(unit,*)
      read(unit,*)
      do i=1,n
         read(unit,*) atmp,(r0(j), j=i*3-2,i*3)
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
      use expec, only: lrenorm,ijob

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
         call alloc_traj(i)

         ! Determine the spawn times
         call getspawntime(i,atmp)

         ! Read the TrajDump files
         call rdtrajdump(i,atmp)         

         ! Read and output the spawn geometries.
         ! We also read the parent/spawn trajectory information here.
         call getspawngeom(atmp,nspawn,i)

         ! If required, read the logfileCLS file
         if (ijob.eq.12) call rdlogfilecls(atmp,nspawn+1,i)         
         
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
      ! If we have hit an empty line, go to the next line...
      if (string.eq.'') goto 20
      ! ... else read any pertinant information
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

    subroutine alloc_traj(itraj)

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
      traj(itraj)%ispawn=0

      ! Spawn times (actually the timestep at which the trajectory
      ! was spawned)
      if (allocated(traj(itraj)%tspawn)) deallocate(traj(itraj)%tspawn)
      allocate(traj(itraj)%tspawn(ntraj))
      traj(itraj)%tspawn=0

      ! Energies at the Gaussian centres
      if (allocated(traj(itraj)%ener)) deallocate(traj(itraj)%ener)
      allocate(traj(itraj)%ener(ntraj,nstep,nsta))
      traj(itraj)%ener=0.0d0
      
      ! NACTs at the Gaussian centres
      if (allocated(traj(itraj)%nact)) deallocate(traj(itraj)%nact)
      allocate(traj(itraj)%nact(ntraj,nstep,nsta,natm*3))
      traj(itraj)%nact=0.0d0

      ! Dipoles and transition dipoles at the Gaussian centres
      if (allocated(traj(itraj)%dipole)) deallocate(traj(itraj)%dipole)
      allocate(traj(itraj)%dipole(ntraj,nstep,3,nsta,nsta))
      traj(itraj)%dipole=0.0d0

!-----------------------------------------------------------------------
! Expectation values
!-----------------------------------------------------------------------
      if (allocated(adpop)) deallocate(adpop)
      allocate(adpop(nsta,nstep))
      adpop=0

      return

    end subroutine alloc_traj

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
                                     n,ncoo,spawnstep
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

         ! Time step at which the trajectory was spawned
         spawnstep=traj(itraj)%tspawn(i)

         ! Open current TrajDump file
         a1=''
         if (i.lt.10) then
            write(a1,'(a9,i1)') 'TrajDump.',i
         else
            write(a1,'(a9,i2)') 'TrajDump.',i
         endif
         a2=trim(atmp)//'/'//trim(a1)
         open(unit,file=a2,form='formatted',status='old')

         ! Skip past the comment line
         read(unit,*)

         ! Read the times corresponding to integer multiples of dt
10       continue
         rtmp=0.0d0
         read(unit,fmat,end=15) t,(rtmp(j),j=1,ncoo),(ptmp(j),j=1,ncoo),&
              gtmp,crtmp,citmp,dum,stmp

         ! Timestep
         n=1+int(t/dt)

         ! Skip if we are not at a full timestep
         ! N.B., this only works if we hit all the full timesteps,
         ! which isn't actually gauranteed by the propagation
         ! algotithm used in the FMS90 code...
         !if (mod(t,dt).ne.0.0d0) goto 10

         ! Skip if we are before the spawn time
         if (n.lt.spawnstep) goto 10

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

    subroutine getspawntime(itraj,atmp)
      
      use sysdef
      use trajdef

      implicit none

      integer            :: itraj,ntraj,unit,i,k
      real*8             :: t
      character(len=150) :: atmp,aspawn

      ntraj=traj(itraj)%ntraj

      unit=20

      ! First trajectory
      traj(itraj)%tspawn(1)=1
      
      ! Spawned trajectories
      do i=2,ntraj
         
         ! Open the Spawn file
         aspawn=''
         aspawn=trim(atmp)//'/'
         k=len_trim(aspawn)
         write(aspawn(k+1:k+6),'(a)') 'Spawn.'
         if (i.lt.10) then
            write(aspawn(k+7:k+7),'(i1)') i
         else
            write(aspawn(k+7:k+8),'(i2)') i
         endif
         open(unit,file=aspawn,form='formatted',status='old')

         ! Read the spawn time
         read(unit,*)
         read(unit,'(F11.2)') t
         traj(itraj)%tspawn(i)=int(t/dt)+1

         if (mod(t,dt).ne.0.0d0) then
            traj(itraj)%tspawn(i)=traj(itraj)%tspawn(i)+1
         endif

          ! Close the Spawn file
         close(unit)
         
      enddo

      return

    end subroutine getspawntime

!#######################################################################

    subroutine rdlogfilecls(currdir,ntraj,itraj)
      
      use expec
      use sysdef
      use trajdef

      implicit none

      integer                   :: ntraj,itraj,unit,trajnum,i,k,&
                                   nproc,nadtype,modus
      integer, dimension(ntraj) :: npnts
      real*8                    :: ftmp
      character(len=150)        :: currdir,string
      character(len=170)        :: filename

!-----------------------------------------------------------------------
! Check to see whether nadtype=1 has been set in col.inp
! If not, then we cannot get the NACTs from the logfileCLS files      
!-----------------------------------------------------------------------
      call checknadtype(currdir,nadtype)

      ! If the job type is the preparation of TS-PSG input, then
      ! we cannot proceed unless nadtype=1
      if (ijob.eq.12.and.nadtype.ne.1) then
         write(6,'(/,2x,a,/)') 'NACTs can only be read from &
              logfileCLS if nadtype=1'
         STOP
      endif
         
!-----------------------------------------------------------------------
! Determine the no. processors used
!-----------------------------------------------------------------------
      unit=20      
      filename=trim(currdir)//'/FMS.out'
      open(unit,file=filename,form='formatted',status='old')
      
      nproc=1      
1     read(unit,'(a)') string
      if (index(string,'Slaves:').eq.0) goto 1      
2     read(unit,'(a)') string
      if (index(string,'Generating initial conditions').eq.0) then
         nproc=nproc+1
         goto 2
      endif

      close(unit)

!-----------------------------------------------------------------------
! Read the energies and properites at the Gaussian centres from the
! logfileCLS files.
!
! modus = 1 <-> Read energies
!         2 <-> Read properties
!
! N.B., due to the way in which the transition dipoles are written to
!       file, the energies must be read first
!-----------------------------------------------------------------------
      do modus=1,2
         
         ! (1) logfileCLS
         filename=trim(currdir)//'/logfileCLS'
         open(unit,file=filename,form='formatted',status='old')
         call rdlogfilecls_1file(modus,0,unit,itraj,nadtype)
         close(unit)
         
         ! (2) logfileCLS.n, n=1,...,nproc
         !
         ! Loop over processor numbers
         do i=2,nproc
            
            ! Open the current logfileCLS file
            filename=trim(currdir)//'/logfileCLS'
            k=len_trim(filename)
            if (i-1.lt.10.and.i-1.ne.0) then
               write(filename(k+1:k+2),'(a1,i1)') '.',i-1
            else if (i-1.ge.10) then
               write(filename(k+1:k+3),'(a1,i2)') '.',i-1
            endif
            open(unit,file=filename,form='formatted',status='old')
            
            ! Read the current logfileCLS file
            call rdlogfilecls_1file(modus,i-1,unit,itraj,nadtype)
            
            ! Close the current logfileCLS file
            close(unit)

         enddo
         
      enddo

      ! ENERGY CHECK
      !do i=1,nstep
      !   if (traj(itraj)%ener(1,i,1).ne.0.0d0) print*,(i-1)*dt,traj(itraj)%ener(1,i,1)
      !enddo
         
      return

    end subroutine rdlogfilecls

!#######################################################################

    subroutine checknadtype(currdir,nadtype)

      use parsemod

      implicit none

      integer                              :: inkw,unit,nadtype
      integer,   parameter                 :: maxkw=200
      integer,   dimension(maxkw)          :: ilkw
      character(len=120), dimension(maxkw) :: keyword
      character(len=150)                   :: currdir
      character(len=170)                   :: filename

!-----------------------------------------------------------------------
! Open the col.inp file
!-----------------------------------------------------------------------
      unit=25
      filename=trim(currdir)//'/col.inp'
      open(unit,file=filename,form='formatted',status='old')

!-----------------------------------------------------------------------
! Determine what nadtype has been set to
!-----------------------------------------------------------------------
10    call rdinp(unit,keyword,inkw,ilkw)
      if (keyword(1).ne.'nadtype') goto 10

      read(keyword(3),*) nadtype

!-----------------------------------------------------------------------
! Close the col.inp file
!-----------------------------------------------------------------------
      close(unit)
      
      return

    end subroutine checknadtype

!#######################################################################
! rdlogfilecls_1file: reads a single logfilecls file.
!                     modus = 1 <-> Read MRCI energies
!                             2 <-> Read properties
!
! Note that, in the current incarnation, the MRCI energies have to be
! read before the transition dipole moments can be read.
!#######################################################################
    
    subroutine rdlogfilecls_1file(modus,iproc,unit,itraj,nadtype)

      use expec
      use sysdef
      use trajdef

      implicit none

      integer                           :: modus,iproc,unit,itraj,&
                                           trajnum,istep,s1,s2,&
                                           trajsta,currsta,i,j,k,&
                                           nadtype,npairs,found
      real*8                            :: time,ftmp,e1,e2
      real*8, dimension(:), allocatable :: tlast_ener,tlast_dm,&
                                           tlast_tdm,tlast_nact
      real*8, dimension(3*natm)         :: tmpvec
      real*8, parameter                 :: dethrsh=1e-6
      character(len=180)                :: string

!-----------------------------------------------------------------------
! Initialisation
!-----------------------------------------------------------------------
      ! Times at which various sections were last read: it is necessary
      ! to keep track of this in order to identify spawning events,
      ! which appear with the same timestamp as the previous 'good'
      ! entry
      allocate(tlast_ener(traj(itraj)%ntraj))
      allocate(tlast_dm(traj(itraj)%ntraj))
      allocate(tlast_tdm(traj(itraj)%ntraj))
      allocate(tlast_nact(traj(itraj)%ntraj))
      tlast_ener=-1.0d0
      tlast_dm=-1.0d0
      tlast_tdm=-1.0d0
      tlast_nact=-1.0d0
      
!-----------------------------------------------------------------------
! If we are reading logfileCLS (not logfileCLS.n), then skip past the
! first 'initialisation' calculation
!-----------------------------------------------------------------------
      if (iproc.eq.0) then

         found=0

5        read(unit,'(a)') string

         if (index(string,'Dalton Summary').ne.0) found=found+1
            
         if (found.lt.2) goto 5
         
      endif
      
!-----------------------------------------------------------------------
! Read the electronic structure information
!-----------------------------------------------------------------------
10    read(unit,'(a)',end=999) string

      !-----------------------------------------------------------------
      ! MRCI energies
      !-----------------------------------------------------------------
      if (modus.eq.1.and.index(string,'CIUDG  Summary').ne.0) then

         ! Determine the trajectory index and only proceed if we
         ! are not at a centroid calculation
         read(string,'(48x,i4)') trajnum
         if (trajnum.gt.0) then

            ! Read the time and determine the timestep
            read(string,'(59x,F10.2)') time


            istep=int(time/dt)+1

            ! Only proceed if we are at a full timestep and we are not
            ! at a spawning event
            if (mod(time,dt).eq.0.0d0.and.time.ne.tlast_ener(trajnum)) then

               ! Read the MRCI energies
               found=0
15             read(unit,'(a)',end=999) string
               if (index(string,'MRCI convergence criteria').eq.0) then

                  if (index(string,'ROOT').ne.0) then
                     found=found+1
                     read(unit,*)
                     read(unit,'(a)') string
                     read(string,'(30x,F16.10)') ftmp
                     traj(itraj)%ener(trajnum,istep,found)=ftmp
                  endif

                  ! Continue reading if we haven't found all the roots yet
                  if (found.lt.nsta) goto 15
               endif

               ! Save the current time in order to identify and skip
               ! spawning events
               tlast_ener(trajnum)=time
               
            endif
            
         endif

         ! Read the next line in the logfileCLS file
         goto 10

      endif

      !-----------------------------------------------------------------
      ! Dipole matrix elements: on-diagonal (electronic component only)
      !-----------------------------------------------------------------
      if (modus.eq.2.and.index(string,'EXPLVL Summary').ne.0) then

         ! Determine the trajectory index and only proceed if we
         ! are not at a centroid calculation
         read(string,'(48x,i4)') trajnum
         if (trajnum.gt.0) then

            ! Read the time and determine the timestep
            read(string,'(69x,F10.2)') time
            istep=int(time/dt)+1

            ! Only proceed if we are at a full timestep and we are not
            ! at a spawning event
            if (mod(time,dt).eq.0.0d0.and.time.ne.tlast_dm(trajnum)) then

               ! Read the state number
               read(string,'(60x,i2)') s1

               ! Read the dipole moment
               do i=1,3
                  read(unit,*)
               enddo
               read(unit,'(a)') string
               read(string,'(11x,3(5x,F11.8))') &
                    (traj(itraj)%dipole(trajnum,istep,i,s1,s1), i=1,3)

               ! Save the current time in order to identify and skip
               ! spawning events
               tlast_dm(trajnum)=time
               
            endif
            
         endif

         ! Read the next line in the logfileCLS file
         goto 10
      endif

      !-----------------------------------------------------------------
      ! Dipole matrix elements: off-diagonal
      !-----------------------------------------------------------------
      if (modus.eq.2.and.index(string,'TRNCI Summary').ne.0) then

         ! Determine the trajectory index and only proceed if we
         ! are not at a centroid calculation
         read(string,'(47x,i4)') trajnum
         if (trajnum.gt.0) then

            ! Read the time and determine the timestep
            read(string,'(58x,F10.2)') time
            istep=int(time/dt)+1

            ! Only proceed if we are at a full timestep and we are not
            ! at a spawning event
            if (mod(time,dt).eq.0.0d0.and.time.ne.tlast_tdm(trajnum)) then
            
               ! Read the state energies and match these to the
               ! corresponding state numbers
               ! N.B. We here assume that the corresponding MRCI energies
               ! have already been written...
               do i=1,5
                  read(unit,*)
               enddo
               read(unit,'(a)') string
               read(string,'(18x,F19.8,F23.8)') e1,e2
               s1=0
               s2=0
               do i=1,nsta
                  if (abs(e1-traj(itraj)%ener(trajnum,istep,i)).lt.dethrsh) s1=i
                  if (abs(e2-traj(itraj)%ener(trajnum,istep,i)).lt.dethrsh) s2=i
               enddo
               if (s1.eq.0.or.s2.eq.0) then

                  print*,trajnum,time
                  
                  write(6,'(/,a,/)') &
                       'Error reading the transition dipole moments...'
                  STOP
               endif
               
               ! Read the transition dipole moment
               do i=1,13
                  read(unit,*)
               enddo
               read(unit,'(a)') string
               read(string,'(13x,3(3x,F10.6))') &
                    (traj(itraj)%dipole(trajnum,istep,i,s1,s2), i=1,3)
               do i=1,3
                  traj(itraj)%dipole(trajnum,istep,i,s2,s1)=&
                       traj(itraj)%dipole(trajnum,istep,i,s1,s2)
               enddo

               ! Save the current time in order to identify and skip
               ! spawning events
               tlast_tdm(trajnum)=time
               
            endif
               
         endif

         ! Read the next line in the logfileCLS file
         goto 10

      endif
      
      !-----------------------------------------------------------------
      ! NACTs: only required if the job type is a TS-PSG preparation
      !-----------------------------------------------------------------
      if (modus.eq.2.and.index(string,'NAD Summary').ne.0.and.ijob.eq.12) then
         print*,"WRITE THE NACT PARSING CODE!"
         STOP
      endif
      
      
      ! Read the next line in the logfileCLS file
      goto 10

      
      ! End of file
999   continue
      
      return

    end subroutine rdlogfilecls_1file

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

      use utils

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

    subroutine getcentroids

      use trajdef
      use centdef
      use sysdef
      use gausstools
      
      implicit none

      integer           :: i,j,n,n1,n2,indx,istep,stepi,stepf
      
!-----------------------------------------------------------------------
! Output where we are at
!-----------------------------------------------------------------------
      write(6,'(/,2x,a,/)') 'Calculating the centroid information...'
      
!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      ! Numbers of centroids per IFG
      allocate(ncent(nintraj))

      ! Determine the maximum possible no. centroids across all
      ! IFGs
      maxcent=0
      do i=1,nintraj
         n=traj(i)%ntraj
         ncent(i)=n*(n-1)/2
         if (ncent(i).gt.maxcent) maxcent=ncent(i)
      enddo

      ! Allocate the cent array of centroid derived types
      allocate(cent(nintraj,maxcent))

!-----------------------------------------------------------------------
! Allocate the arrays of the cetroid derived types
!-----------------------------------------------------------------------
      do i=1,nintraj
         do j=1,ncent(i)
            call alloc_cent(i,j)
         enddo
      enddo

!-----------------------------------------------------------------------
! Fill in the centroid derived types
!-----------------------------------------------------------------------
      ! Loop over IFGs
      do i=1,nintraj
         
         ! Loop over the 1st trajectory index
         do n1=1,traj(i)%ntraj-1
            
            ! Loop over the 2nd trajectory index
            do n2=n1+1,traj(i)%ntraj
               
               ! Get the centroid index for the current pair
               ! of trajectories
               indx=centindx(n1,n2,i)
               
               ! Trajectory indices: lower-triangle, column-major order
               ! convention
               cent(i,indx)%indx1=n2
               cent(i,indx)%indx1=n1
               
               ! Determine the timestep interval for which both
               ! trajectories are alive
               call get_interval_alive_2traj(i,n1,n2,stepi,stepf)
               cent(i,indx)%tspawn=stepi
               cent(i,indx)%tkill=stepf

               ! If the two trajectories are never simultaneously alive,
               ! then flag the centroid as inactive this and cycle
               if (stepi.gt.stepf) then
                  cent(i,indx)%lcontrib=.false.
                  cycle
               else
                  cent(i,indx)%lcontrib=.true.
               endif
               
               ! Loop over the timesteps for which both trajectories
               ! are alive
               do istep=stepi,stepf
                  
                  ! Centroid position vector
                  cent(i,indx)%r(istep,:)=rcent(n1,n2,i,istep)
                  
                  ! Centroid momentum vector
                  cent(i,indx)%p(istep,:)=pcent(n1,n2,i,istep)
                  
                  ! Overlap of the trajectories
                  cent(i,indx)%ovrlp(istep)=&
                       overlap_general_nuc_only(i,i,n1,n2,istep,istep)
                  
               enddo

               ! If the overlap of the two trajectories if never above
               ! threshold, then flag the centroid as inactive
               if (maxval(abs(cent(i,indx)%ovrlp)).lt.othrsh) &
                    cent(i,indx)%lcontrib=.false.
            enddo
            
         enddo

      enddo
         
      return
      
    end subroutine getcentroids

!#######################################################################

    subroutine alloc_cent(ifg,cindx)

      use sysdef
      use centdef
      use expec
      
      implicit none

      integer :: ifg,cindx

!-----------------------------------------------------------------------
! Centroid information
!-----------------------------------------------------------------------
      ! Positions
      allocate(cent(ifg,cindx)%r(nstep,natm*3))
      cent(ifg,cindx)%r=0.0d0

      ! Momenta
      allocate(cent(ifg,cindx)%p(nstep,natm*3))
      cent(ifg,cindx)%p=0.0d0

      ! Widths
      allocate(cent(ifg,cindx)%alpha(natm*3))
      cent(ifg,cindx)%alpha=0.0d0

      ! Overlaps
      allocate(cent(ifg,cindx)%ovrlp(nstep))
      cent(ifg,cindx)%ovrlp=0.0d0
      
      return
      
    end subroutine alloc_cent

!#######################################################################
! get_interval_alive_2traj: determines the first and last timesteps
!                           (stepi and stepf) at which two trajectories
!                           are simultaneously alive
!#######################################################################
    
    subroutine get_interval_alive_2traj(ifg,n1,n2,stepi,stepf)

      use trajdef
      
      implicit none

      integer :: ifg,n1,n2,stepi,stepf,tspawn1,tspawn2,tkill1,tkill2

      tspawn1=traj(ifg)%tspawn(n1)
      tspawn2=traj(ifg)%tspawn(n2)      

      tkill1=traj(ifg)%tkill(n1)
      tkill2=traj(ifg)%tkill(n2)

      stepi=max(tspawn1,tspawn2)
      stepf=min(tkill1,tkill2)
      
      return
      
    end subroutine get_interval_alive_2traj
      
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
      use adcspec
      use tspsgmod
      use gamessprep
      use gamesstrxas
      use dipoleprep
      use dipole
      use adc_orben
      
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
!            calculated using columbus/superdyson
!      7 <-> Preparation of input files for ADC absorption or
!            photoionization cross-section calculations
!      8 <-> Calculation of the bound part of the TRXAS using 
!            ADC cross-sections
!      9 <-> Calculation of the continuum part of the TRXAS using 
!            ADC cross-sections
!     10 <-> Two-dimensional reduced densities
!     11 <-> Calculation of TRPES using Dyson orbital norms
!            calculated using the ADC program
!     12 <-> Preparation of input for a TS-PSG dynamics calculation
!     13 <-> Preparation of input for Columbus-MRCI/GAMESS-ORMAS X-ray
!            absorption cross-section calculations
!     14 <-> Calculation of the bound part of the TRXAS using
!            Columbus-MRCI/GAMESS-ORMAS cross-sections
!     15 <-> Preparation of input for the calculation of dipole
!            matrix elements (a pre-requesite for the calculation
!            dipole expectation values)
!     16 <-> Calculation of dipole expectation values
!     17 <-> Incoherently averaged ADC orbital energies
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
      else if (ijob.eq.8.or.ijob.eq.9.or.ijob.eq.11) then
         call adc_spec
      else if (ijob.eq.10) then
         call reddens_2d
      else if (ijob.eq.12) then
         call tspsg_prep
      else if (ijob.eq.13) then
         call mkgamessinp
      else if (ijob.eq.14) then
         call gamess_trxas
      else if (ijob.eq.15) then
         call prep_dipole_inp
      else if (ijob.eq.16) then
         call calc_dipole
      else if (ijob.eq.17) then
         call calc_adc_orben
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

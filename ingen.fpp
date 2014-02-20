program ingen

! 
! Reads in namelists, checks for consistency of inputs, writes a report
! to runname.report, checks for requested variations, 
! generates new input files, and exits.
! 
! Consistency checks/reports:
!
! wstar_units incompatible with nonlinear
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!                                  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!                                  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!        Declarations              !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!                                  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!                                  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use mp, only: init_mp, finish_mp
  use constants 
  use file_utils
  use text_options, only: text_option, get_option_value
  implicit none
  
  character (100) :: pythonin
  integer :: interactive_record, interactive_input

  integer :: in_file, i, ierr, unit, is, report_unit, iunit, ncut, npmax
  logical :: exist, scan, stdin
  logical :: coll_on = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!                                  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!                                  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!          Namelists               !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!                                  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!                                  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  namelist /ingen_knobs/ ncut, scan, stdin, npmax

!CMR
  logical, parameter :: debug=.false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!                                  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!                                  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!      Main code starts here       !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!                                  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!                                  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call init_mp
 
  if (debug) write(6,*) 'ingen: call get_namelists'
  
  call get_namelists
  if (debug) write(6,*) 'ingen: call report'
  call report
  if (debug) write(6,*) 'ingen: call write_namelists'
  call write_namelists

  if (debug) write(6,*) 'ingen: call interactive, scan=',scan
  if (scan) call interactive

  call finish_mp

contains

  subroutine interactive
    use species, only: spec, nspec, has_electron_species
    use geometry, only: beta_prime_input, bishop
    use run_parameters, only: beta, fapar, fbpar
    integer :: sel, nbeta, j, ise
    real :: beta_low, beta_high, dbeta, beta_save
    real :: fapar_save, fbpar_save, pri, pe, alpi, tpe_save, ptot, alp, dbdr
    real :: alt, aln, fac, beta_prime_save, bishop_save
    real, dimension (:), allocatable :: tp_save, fp_save
    character (500) :: tag1, tag2
    logical :: first = .true.

    if (first) then
       call get_unused_unit (interactive_record)
       open (unit=interactive_record, file='.'//trim(run_name)//".record")
       first = .false.

       if (.not. stdin) then
          call get_unused_unit (interactive_input)
          open (unit=interactive_input, file=trim(pythonin))
       else
          interactive_input = 5
       end if
    end if

    call tell ('Interactive specification of a parameter scan')

100 continue
    
    call text
    call text ('Choose a parameter that you would like to vary (1-6):')
    call text ('(1) beta            (4) temperature gradient')
    call text ('(2) beta_prime      (5) density gradient')
    call text ('(3) collisionality  (6) Z_effective')
    call text
    call get_choice (6, sel)
    
    select case (sel)
       
    case default
       call tell ('Try again.  Choose an integer between 1 and 6, inclusively.')
       goto 100

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  1.0  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    case (1) 
       call tell ('You have chosen to vary beta.')

101    continue

       call text
       call text ('Choose from the following:')
       call text ('(1) Vary beta self-consistently')
       call text ('(2) Vary beta with all other parameters held fixed (not self-consistently).')
       call text
       call get_choice (2, sel)

       select case (sel)

       case default
          call tell ('Try again.  Choose an integer between 1 and 2, inclusively.')
          goto 101

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  1.1  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       case (1) 
          call tell ('You have chosen to vary beta self-consistently.')

102       continue
          call text
          call text ('Choose from the following:')
          call text ('(1) Hold beta_prime fixed, vary electron temperature gradient scale length')
          call text ('(2) Hold beta_prime fixed, vary all temperature gradient scale lengths by same factor')
          call text ('(3) Hold beta_prime fixed, vary all density gradient scale lengths by same factor')
          call text

          call get_choice (3, sel)

          select case (sel)
             
          case default
             call tell ('Try again.  Choose an integer between 1 and 2, inclusively.')
             goto 102
             
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  1.1.1  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          case (1)  
             call tell ('You have chosen to vary beta and electron temperature gradient at fixed beta_prime.')

             call beta_range_low (beta_low, 'le', 0.)
             call beta_range_high (beta_high, 'le', beta_low)
             call num_runs (nbeta)

             call tell ('Preparing a self-consistent beta scan at fixed beta_prime.', &
                  'The electron temperature gradient scale length will be varied', &
                  'to maintain consistency.')

             call run_number (sel, nbeta)

             write (tag1, fmt='(i3," runs prepared with beta_min = ",e16.10,&
                  &" and beta_max = ",e16.10)') nbeta, beta_low, beta_high
             write (tag2, fmt='("Files are numbered from ",i3," to ",i3)') sel, sel+nbeta-1

             call tell (tag1, tag2)

             ptot = 0.
             alp = 0.
             pri = 0.
             pe = 0.
             alpi = 0.
             ise = 0
             do is=1,nspec
                if (spec(is)%type == 2) then
                   pe = spec(is)%dens * spec(is)%temp
                   ise = is
                else
                   pri = pri + spec(is)%dens * spec(is)%temp
                   alpi = alpi + spec(is)%dens * spec(is)%temp *(spec(is)%fprim + spec(is)%tprim)
                endif
                ptot = ptot + spec(is)%dens * spec(is)%temp
                alp = alp + spec(is)%dens * spec(is)%temp *(spec(is)%fprim + spec(is)%tprim)
             end do
             
             if (.not. has_electron_species(spec)) call tell ('You really should use electrons for electromagnetic runs.')

             alp = alp/ptot
             dbdr = - beta*ptot*alp

             dbeta = (beta_high-beta_low)/(nbeta-1)
             do j = sel, sel+nbeta-1
                
                beta_save = beta
                beta = beta_low + (j - sel)*dbeta
                
                tpe_save = spec(ise)%tprim
                spec(ise)%tprim = - (spec(ise)%fprim + alpi/pe + dbdr/beta/pe)

                write (tag1, fmt='("Varying beta and L_Te self-consistently with& 
                     & beta_prime fixed")') 

                write (tag2, fmt='("beta = ",e16.10," and electron tprim = ",e16.10)') beta, spec(ise)%tprim 

                call write_namelists (j, tag1, tag2)
                spec(ise)%tprim = tpe_save
             end do
             beta = beta_save

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  1.1.2  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          case (2)

             call tell ('You have chosen to vary beta and all temperature &
                  &gradient scale lengths together at fixed beta_prime.')

             call beta_range_low (beta_low, 'le', 0.)
             call beta_range_high (beta_high, 'le', beta_low)
             call num_runs (nbeta)

             call tell ('Preparing a self-consistent beta scan at fixed beta_prime.', &
                  'All temperature gradient scale lengths will be varied', &
                  'by the same factor to maintain consistency.')

             call run_number (sel, nbeta)

             write (tag1, fmt='(i3," runs prepared with beta_min = ",e16.10,&
                  &" and beta_max = ",e16.10)') nbeta, beta_low, beta_high
             write (tag2, fmt='("Files are numbered from ",i3," to ",i3)') sel, sel+nbeta-1

             call tell (tag1, tag2)

             allocate (tp_save (nspec))

             ptot = 0.
             alt = 0.
             aln = 0.
             do is=1,nspec
                ptot = ptot + spec(is)%dens * spec(is)%temp
                alt = alt + spec(is)%dens * spec(is)%temp *(spec(is)%tprim)
                aln = aln + spec(is)%dens * spec(is)%temp *(spec(is)%fprim)
             end do
             
             if (.not. has_electron_species(spec)) call tell ('You really should use electrons for electromagnetic runs.')

             alp = (alt+aln)/ptot
             dbdr = - beta*ptot*alp

             dbeta = (beta_high-beta_low)/(nbeta-1)
             do j = sel, sel+nbeta-1
                
                beta_save = beta
                beta = beta_low + (j - sel)*dbeta
                
                fac = -(dbdr/beta+aln)/alt
                tp_save = spec%tprim
                spec%tprim = fac*spec%tprim

                write (tag1, fmt='("Varying beta and all L_T values self-consistently")')
                write (tag2, fmt='("beta = ",e16.10," and tprim values scaled by ",e16.10)') beta, fac

                call write_namelists (j, tag1, tag2) 
                spec%tprim = tp_save
             end do
             beta = beta_save

             deallocate (tp_save)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  1.1.3  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          case (3)

             call tell ('You have chosen to vary beta and all density &
                  &gradient scale lengths together at fixed beta_prime.')

             call beta_range_low (beta_low, 'le', 0.)
             call beta_range_high (beta_high, 'le', beta_low)
             call num_runs (nbeta)

             call tell ('Preparing a self-consistent beta scan at fixed beta_prime.', &
                  'All density gradient scale lengths will be varied', &
                  'by the same factor to maintain consistency.')

             call run_number (sel, nbeta)

             write (tag1, fmt='(i3," runs prepared with beta_min = ",e16.10,&
                  &" and beta_max = ",e16.10)') nbeta, beta_low, beta_high
             write (tag2, fmt='("Files are numbered from ",i3," to ",i3)') sel, sel+nbeta-1

             call tell (tag1, tag2)

             allocate (fp_save (nspec))

             ptot = 0.
             alt = 0.
             aln = 0.
             do is=1,nspec
                ptot = ptot + spec(is)%dens * spec(is)%temp
                alt = alt + spec(is)%dens * spec(is)%temp *(spec(is)%tprim)
                aln = aln + spec(is)%dens * spec(is)%temp *(spec(is)%fprim)
             end do
             
             if (.not. has_electron_species(spec)) call tell ('You really should use electrons for electromagnetic runs.')

             alp = (alt+aln)/ptot
             dbdr = - beta*ptot*alp

             dbeta = (beta_high-beta_low)/(nbeta-1)
             do j = sel, sel+nbeta-1
                
                beta_save = beta
                beta = beta_low + (j - sel)*dbeta
                 
                fac = -(dbdr/beta+alt)/aln
                fp_save = spec%fprim
                spec%fprim = fac*spec%fprim

                write (tag1, fmt='("Varying beta and all L_n values self-consistently")')
                write (tag2, fmt='("beta = ",e16.10," and tprim values scaled by ",e16.10)') beta, fac

                call write_namelists (j, tag1, tag2) 
                spec%fprim = fp_save
             end do
             beta = beta_save

             deallocate (fp_save)
          end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  1.2  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       case (2)
          call tell ('You have selected to vary beta non-self-consistently.')
          
          call beta_range_low (beta_low, 'lt', 0.)
          call beta_range_high (beta_high, 'le', beta_low)
          call num_runs (nbeta)

          call tell ('Preparing a non-self-consistent beta scan.')

          call run_number (sel, nbeta)

          write (tag1, fmt='(i3," runs prepared with beta_min = ",e16.10,&
               &" and beta_max = ",e16.10)') nbeta, beta_low, beta_high
          write (tag2, fmt='("Files are numbered from ",i3," to ",i3)') sel, sel+nbeta-1
          
          call tell (tag1, tag2)

          dbeta = (beta_high-beta_low)/(nbeta-1)
          do j = sel, sel+nbeta-1
             
             beta_save = beta
             fapar_save = fapar 
             fbpar_save = fbpar

             beta = beta_low + (j - sel)*dbeta
             if (beta == 0.) then 
                fapar = 0.
                fbpar = 0.
             else
                if (fapar == 0. .and. fbpar == 0.) then
                   fapar = 1.0 ;  fbpar = 1.0
                end if
             end if
             
             write (tag1, fmt='("Varying beta, all else fixed")')
             write (tag2, fmt='("beta = ",e16.10)') beta

             call write_namelists (j, tag1, tag2)

             fapar = fapar_save 
             fbpar = fbpar_save
             beta = beta_save
          end do
          
       end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  2.0  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    case (2) ! beta_prime
       
       call tell ('You have chosen to vary beta_prime.')

115    continue
       call text
       call text ('Choose from the following:')
       call text ('(1) Vary beta_prime self-consistently')
       call text ('(2) Vary beta_prime with ALL other parameters held fixed (non-self-consistently).')
       call text
       call get_choice (2, sel)

       select case (sel)

       case default
          call tell ('Try again.  Choose an integer between 1 and 2, inclusively.')
          goto 115

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  2.1  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       case (1)
          call tell ('You have chosen to vary beta_prime self-consistently.')

116       continue
          call text
          call text ('Choose from the following:')
          call text ('(1) Hold gradient scale lengths fixed, vary beta')
          call text

          call get_choice (1, sel)

          select case (sel)

          case default
             call tell ('Try again.  Choose an integer between 1 and 1, inclusively.')
             goto 116

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  2.1.1  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          case (1)  
             call tell ('You have chosen to vary beta_prime while holding gradient scale lengths fixed.')

             call beta_prime_range_low (beta_low)
             call beta_prime_range_high (beta_high, beta_low)
             call num_runs (nbeta)

             call tell ('Preparing a self-consistent beta_prime scan.', &
                  'Beta will be varied to maintain consistency.')

             call run_number (sel, nbeta)

             write (tag1, fmt='(i3," runs prepared with beta_prime_min = ",e16.10,&
                  &" and beta_prime_max = ",e16.10)') nbeta, beta_low, beta_high
             write (tag2, fmt='("Files are numbered from ",i3," to ",i3)') sel, sel+nbeta-1
             
             call tell (tag1, tag2)

             ptot = 0.
             alp = 0.
             do is=1,nspec
                ptot = ptot + spec(is)%dens * spec(is)%temp
                alp = alp + spec(is)%dens * spec(is)%temp *(spec(is)%fprim + spec(is)%tprim)
             end do
            
             alp = alp/ptot
             dbdr = - beta*ptot*alp

             if (alp == 0.) then
                call tell ('Cannot proceed, because Lp = infinity', &
                     'No input files for scan written')
                return
             end if

             beta_save = beta
             beta_prime_save = dbdr

             fac = -1./(ptot*alp)

             dbeta = (beta_high-beta_low)/(nbeta-1)   ! actually, this is dbeta_prime
             do j = sel, sel+nbeta-1
                
                beta_prime_save = beta_prime_input
                beta_prime_input = beta_low + (j - sel)*dbeta
                
                beta_save = beta
                beta = beta_prime_input*fac

                fapar_save = fapar ; fbpar_save = fbpar
                if (beta == 0.) then
                   fapar = 0.      ; fbpar = 0.
                else
                   if (fapar == 0. .and. fbpar == 0.) then
                      fapar = 1.0 ;  fbpar = 1.0
                   end if
                end if

                select case (bishop)
                case default
                   bishop_save = bishop
                   bishop = 6
                case (4)
                   ! nothing, continue to use bishop = 4
                end select

             
                write (tag1, fmt='("Varying beta_prime and beta self-consistently")')
                write (tag2, fmt='("beta_prime = ",e16.10," and beta = ",e16.10)') beta_prime_input, beta

                call write_namelists (j, tag1, tag2)

                fapar = fapar_save 
                fbpar = fbpar_save
                beta = beta_save
                beta_prime_input = beta_prime_save

                select case (bishop)
                case default
                   bishop = bishop_save
                case (4)
                   ! nothing, continue to use bishop = 4
                end select

             end do
          end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  2.2   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       case (2)

          call tell ('You have selected to vary beta_prime non-self-consistently.')          
          
          call beta_prime_range_low (beta_low)
          call beta_prime_range_high (beta_high, beta_low)
          call num_runs (nbeta)

          call tell ('Preparing a non-self-consistent beta_prime scan.')

          call run_number (sel, nbeta)
          
          write (tag1, fmt='(i3," runs prepared with beta_prime_min = ",e16.10,& 
               &" and beta_prime_max = ",e16.10)') nbeta, beta_low, beta_high 
          write (tag2, fmt='("Files are numbered from ",i3," to ",i3)') sel, sel+nbeta-1 
          
          call tell (tag1, tag2)

          dbeta = (beta_high-beta_low)/(nbeta-1) 
          do j = sel, sel+nbeta-1
             
             beta_prime_save = beta_prime_input
             beta_prime_input = beta_low + (j - sel)*dbeta
             
             select case (bishop)
             case default
                bishop_save = bishop
                bishop = 6
             case (4)
                ! nothing, continue to use bishop = 4
             end select
             
             write (tag1, fmt='("Varying beta_prime only (non-self-consistently)")')
             write (tag2, fmt='("beta_prime = ",e16.10)') beta_prime_input

             call write_namelists (j, tag1, tag2)
             
             beta_prime_input = beta_prime_save
             
             select case (bishop)
             case default
                bishop = bishop_save
             case (4)
                ! nothing, continue to use bishop = 4
             end select
             
          end do
          
       end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   1.3  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    case (3) ! collisionality
       call tell ('You have chosen to vary collisionality.')
       call text ('Not yet implemented (sorry).')
       call text

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   1.4  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    case (4) ! temperature gradient
       call tell ('You have chosen to vary temperature gradient.')
       call text ('Not yet implemented (sorry).')
       call text

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   1.5  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    case (5) ! density gradient
       call tell ('You have chosen to vary density gradient.')
       call text ('Not yet implemented (sorry).')
       call text

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   1.6  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    case (6) ! Z_effective
       call tell ('You have chosen to vary Z_effective.')
       call text ('Not yet implemented (sorry).')
       call text
       
    end select

  end subroutine interactive

  subroutine get_namelists
!CMR, 17/11/2009: use gs2_diagnostics module to pick up public variables
!                 and routines to reduces maintenance of ingen.
!                 Should replicate this for other modules in future.
!CMR, 2/2/2011:   Have extended this to include use of theta_grid module:
!                 Strategy is simply to add two types of routines to modules:
!                      wnml_xxxxx   to write the modules namelists
!                      check_xxxxx  to perform the ingen checking inside the module
!                 More object oriented strategy, easier maintenance of ingen.f90 
!                 which is gradually shrinking.!                  
!
    use antenna, only: init_antenna
    use collisions, only: coll_read_parameters=>read_parameters
    use collisions, only: collision_model_switch, collision_model_none
    use collisions, only: collision_model_full, collision_model_ediffuse
    use collisions, only: collision_model_lorentz, collision_model_lorentz_test
    use fields, only : fields_read_parameters=>read_parameters
    use hyper, only : hyper_read_parameters=>read_parameters
    use species,only : spec, nspec
    use gs2_layouts, only: init_gs2_layouts
    use gs2_reinit, only: init_reinit
    use dist_fn, only: dist_fn_read_parameters=>read_parameters
    use nonlinear_terms, only: nonlinear_terms_read_parameters=>read_parameters
    use run_parameters, only: init_run_parameters
    use init_g, only: init_init_g
    use species, only: init_species, nspec
    use gs2_diagnostics,only: gs2diag_read_parameters=>read_parameters
    use theta_grid, only: init_theta_grid
    use kt_grids, only: init_kt_grids
    use le_grids, only: init_le_grids
    use theta_grid, only: init_theta_grid, nbset, ntgrid
    use constants, only: pi
    implicit none
    logical :: list, accelx=.false., accelv=.false.
    logical, parameter :: debug=.false.
    integer :: is
    call init_file_utils (list, name="template")
if (debug) write(6,*) 'get_namelists: called init_file_utils'
    ncut= 100000
    npmax=100000
    scan = .false.
    stdin = .true.
    pythonin = "."//trim(run_name)//".pythonin"
    in_file=input_unit_exist("ingen_knobs", exist)
    if (exist) read(unit=input_unit("ingen_knobs"), nml=ingen_knobs)

if (debug) write(6,*) 'get_namelists: if (scan), scan=',scan
    if (scan) then
       if (.not. stdin) then
          if (pythonin == "") then
             write (*,*) 'Need to specify pythonin in ingen_knobs.'
          end if 
       end if 
    end if 

if (debug) write(6,*) 'get_namelists: layouts'
    call init_gs2_layouts
    call init_antenna

    ! run_parameters:
    call init_run_parameters

if (debug) write(6,*) 'get_namelists: collisions'
    call coll_read_parameters

if (debug) write(6,*) 'get_namelists: init_g'
    call init_init_g

if (debug) write(6,*) 'get_namelists: fields'
    call fields_read_parameters
    
if (debug) write(6,*) 'get_namelists: gs2_reinit'
    call init_reinit

if (debug) write(6,*) 'get_namelists: hyper'
    call hyper_read_parameters

if (debug) write(6,*) 'get_namelists: kt_grids'
    call init_kt_grids

if (debug) write(6,*) 'get_namelists: le_grids'
    call init_le_grids(accelx, accelv)

if (debug) write(6,*) 'get_namelists: nonlinear terms'
    call nonlinear_terms_read_parameters

if (debug) write(6,*) 'get_namelists: init_species'
    call init_species

    do is=1, nspec
        coll_on = spec(is)%vnewk > epsilon(0.0) .or. coll_on
    end do
    coll_on = coll_on .and. (collision_model_switch /= collision_model_none)

!CMR, 2/2/2011:  reduce much duplication by calling init_theta_grid
if (debug) write(6,*) 'get_namelists: call init_theta_grid, ntgrid=',ntgrid
    call init_theta_grid
if (debug) write(6,*) 'get_namelists: done init_theta_grid, ntgrid=',ntgrid

    call dist_fn_read_parameters
    ! gs2_diagnostics

! CMR 18/11/2009:  reduce duplication by calling gs2diag_read_parameters 
!
   call gs2diag_read_parameters(.false.)

if (debug) write(6,*) 'get_namelists: returning'

  end subroutine get_namelists

  subroutine write_namelists (jr, tag1, tag2)
    use antenna, only: wnml_antenna
    use collisions, only: wnml_collisions
    use dist_fn, only: wnml_dist_fn, wnml_dist_fn_species 
    use fields, only: wnml_fields
    use gs2_reinit, only: wnml_gs2_reinit
    use gs2_layouts, only: wnml_gs2_layouts
    use gs2_diagnostics, only: wnml_gs2_diagnostics
    use hyper, only: wnml_hyper
    use init_g, only : wnml_init_g
    use kt_grids, only: wnml_kt
    use le_grids, only: wnml_le_grids
    use nonlinear_terms, only: nonlin, wnml_nonlinear_terms
    use run_parameters, only: wnml_run_parameters
    use species, only: wnml_species, nspec, spec, has_electron_species
    use theta_grid, only: wnml_theta_grid
    use theta_grid_params, only: wnml_theta_grid_params

    integer, intent (in), optional :: jr
    character (*), intent (in), optional :: tag1, tag2
    logical :: electrons

    integer :: h, t, u
    integer :: i
    character (4) :: suffix
    character(20) :: datestamp, timestamp, zone
    
    call get_unused_unit (unit)

    if (present(jr)) then

       h = jr / 100
       t = (jr - h * 100) / 10
       u = (jr - h * 100 - t * 10)
       suffix = '_'//achar(48+h)//achar(48+t)//achar(48+u)
       open (unit=unit, file=trim(run_name)//suffix//".in")
    else
       open (unit=unit, file=trim(run_name)//".inp")
    endif

    write (unit, *)

    write (unit, fmt="('gs2')")
    datestamp(:) = ' '
    timestamp(:) = ' '
    zone(:) = ' '
    call date_and_time (datestamp, timestamp, zone)
    write (unit=unit, fmt="('Date: ',a,' Time: ',a,1x,a)") &
         trim(datestamp), trim(timestamp), trim(zone)

    if (present(tag1)) then
       write (unit, *) '*****************************************************'
       write (unit, *) trim(tag1)
       if (present(tag2)) write (unit, *) trim(tag2)
       write (unit, *) '*****************************************************'
    end if

    call wnml_gs2_layouts(unit)
    call wnml_collisions(unit)
    call wnml_init_g(unit)
    call wnml_dist_fn(unit)
    call wnml_fields(unit)
    call wnml_gs2_diagnostics(unit)
    if (nonlin) call wnml_gs2_reinit(unit)
    call wnml_hyper(unit)
    call wnml_kt(unit)
    call wnml_le_grids(unit)
    call wnml_nonlinear_terms(unit)
    electrons=has_electron_species(spec)
    call wnml_run_parameters(unit,electrons,coll_on)
    call wnml_species(unit)
    call wnml_dist_fn_species(unit)

    call wnml_theta_grid_params(unit)
    call wnml_theta_grid(unit)

    call wnml_antenna(unit)

    write(unit, fmt=*)
    close (unit)

   end subroutine write_namelists

   subroutine factors (n, j, div)
!CMR, 3/10/13:
!    Find all the factors of n and return in div(j)
     integer, intent (in) :: n
     integer, intent (out) :: j
     integer, dimension (:), intent (out) :: div
     integer :: i, imax

! find: i = lowest factor of n
! and therefore imax=n/1 is the HIGHEST factor of n
     do i=2,n
        if (mod(n,i)==0) exit
     end do
     imax = n/i
     j=1
! loop over all possible factors of n, and return in div(j)
     do i=1,imax
        if (mod(n,i)==0) then
           div(j) = i
           j=j+1
        end if
     end do
     div(j) = n
   end subroutine factors

   subroutine pfactors (n, div)
     integer, intent (in) :: n
     integer, dimension (:), intent (out) :: div
     integer, dimension (50), parameter :: primes = (/ &
          2, 3, 5, 7, 11, &
          13, 17, 19, 23, 29, &
          31, 37, 41, 43, 47, &
          53, 59, 61, 67, 71, &
          73, 79, 83, 89, 97, &
          101, 103, 107, 109, 113, &
          127, 131, 137, 139, 149, &
          151, 157, 163, 167, 173, &
          179, 181, 191, 193, 197, &
          199, 211, 223, 227, 229 /)

     integer :: i, ntmp

     ntmp = n 
     i=1
     do while (ntmp > 1 .and. i < 51)
        do while (mod(ntmp, primes(i)) == 0)
           if (i < 4) div(i) = div(i) + 1
           if (i > 3) div(4) = primes(i)
           ntmp = ntmp / primes(i)
        end do
        i=i+1
     end do
   end subroutine pfactors

   subroutine report
     use antenna, only: check_antenna
     use collisions, only: check_collisions
     use dist_fn, only: check_dist_fn
     use fields, only: check_fields
     use gs2_diagnostics, only: check_gs2_diagnostics
     use gs2_diagnostics, only: dump_fields_periodically, save_for_restart
     use gs2_diagnostics, only: nsave, make_movie, nmovie, exit_when_converged, nwrite, omegatol
     use gs2_reinit, only: delt_adj
     use gs2_time, only: code_dt, user_dt, code_dt_min
     use hyper, only: check_hyper
     use init_g, only: check_init_g
     use kt_grids, only: check_kt_grids, grid_option, gridopt_switch
     use kt_grids, only: gridopt_box, naky, ntheta0, nx, ny
!     use le_grids, only: leok_le_grids, check_le_grids
     use le_grids, only: negrid, nlambda
     use nonlinear_terms, only: nonlin, cfl, check_nonlinear_terms
     use run_parameters, only: check_run_parameters
     use run_parameters, only: beta, tite, margin, code_delt_max
     use run_parameters, only: nstep, wstar_units
     use species, only: check_species, spec, nspec, has_electron_species
     use theta_grid, only: check_theta_grid
     use theta_grid, only: gb_to_cv, nbset, ntgrid
     use theta_grid_params, only: nperiod, ntheta, eps, epsl, rmaj, r_geo
     use theta_grid_params, only: pk, qinp, rhoc, shift, shat
     use theta_grid_params, only: akappa, akappri, tri, tripri
     use collisions, only: collision_model_switch, collision_model_none
     use collisions, only: collision_model_full, collision_model_ediffuse
     use collisions, only: collision_model_lorentz, collision_model_lorentz_test
     use collisions, only: use_le_layout

     implicit none
     real :: alne, dbetadrho_spec
     real :: kxfac, drhodpsi
     character (20) :: datestamp, timestamp, zone
     character (200) :: line
     logical :: le_ok = .true.
     integer :: j, nmesh
     integer, dimension(4) :: pfacs

     call get_unused_unit (report_unit)
     call open_output_file (report_unit, ".report")

     write (report_unit, *) 

     write (report_unit, fmt="('GS2')")
     datestamp(:) = ' '
     timestamp(:) = ' '
     zone(:) = ' '
     call date_and_time (datestamp, timestamp, zone)
     write (report_unit, fmt="('Date: ',a,'/',a,'/',a,&
          &  ' Time: ',a,':',a,1x,a)") datestamp(5:6), &
          &  datestamp(7:8), datestamp(1:4), timestamp(1:2), timestamp(3:4), trim(zone)

     write (report_unit, fmt="(/'------------------------------------------------------------'/)")

!     le_ok=leok_le_grids(report_unit)
     le_ok = .true.  ! temporary fix until leok_le_grids is restored to the distibution.  BD (sorry about that)
     
     if (le_ok) then

        if (nonlin) then
           if (gridopt_switch /= gridopt_box) then
              write (report_unit, *) 
              write (report_unit, fmt="('################# WARNING #######################')")
              write (report_unit, fmt="('Nonlinear runs must be carried out in a box.')") 
              write (report_unit, fmt="('Set grid_option to box in the kt_grids_knobs namelist')") 
              write (report_unit, fmt="('or set nonlinear_mode to off in the nonlinear_knobs namelist.')") 
              write (report_unit, fmt="('THIS IS AN ERROR.')") 
              write (report_unit, fmt="('################# WARNING #######################')")
              write (report_unit, *) 
           end if

           write (report_unit, fmt="('nmesh=(2*ntgrid+1)*2*nlambda*negrid*nx*ny*nspec')")
           nmesh = (2*ntgrid+1)*2*nlambda*negrid*nx*ny*nspec
           write (report_unit, fmt="('Number of meshpoints:    ',i16)") nmesh

 !
 ! check that nx, ny have no large prime factors
 !
         pfacs = 0
           call pfactors (ny, pfacs)
           if (pfacs(4) /= 0) then
              write (report_unit, *) 
              write (report_unit, fmt="('################# WARNING #######################')")
              write (report_unit, fmt="('ny is a multiple of ',i4)") pfacs(4)
              write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')")
              write (report_unit, fmt="('ny should have only 2, 3, and/or 5 as prime factors.')")
              write (report_unit, fmt="('################# WARNING #######################')")
              write (report_unit, *) 
           end if
           i = 1
           if (pfacs(1) > 0) i=2**pfacs(1)
           if (pfacs(2) > 0) i=3**pfacs(2)*i
           if (pfacs(3) > 0) i=5**pfacs(3)*i
           if (i /= ny) then
              write (report_unit, *) 
              write (report_unit, fmt="('################# WARNING #######################')")
              write (report_unit, fmt="('ny = ',i3)") ny
              write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')")
              write (report_unit, fmt="('ny should have only 2, 3, and/or 5 as prime factors.')")
              write (report_unit, fmt="('################# WARNING #######################')")
              write (report_unit, *) 
           end if

           pfacs = 0
           call pfactors (nx, pfacs)
           if (pfacs(4) /= 0) then
              write (report_unit, *) 
              write (report_unit, fmt="('################# WARNING #######################')")
              write (report_unit, fmt="('nx is a multiple of ',i3)") pfacs(4)
              write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')")
              write (report_unit, fmt="('nx should have only 2, 3, and/or 5 as prime factors.')")
              write (report_unit, fmt="('################# WARNING #######################')")
              write (report_unit, *) 
           end if
           i = 1
           if (pfacs(1) > 0) i=2**pfacs(1)
           if (pfacs(2) > 0) i=3**pfacs(2)*i
           if (pfacs(3) > 0) i=5**pfacs(3)*i
           if (i /= nx) then
              write (report_unit, *) 
              write (report_unit, fmt="('################# WARNING #######################')")
              write (report_unit, fmt="('nx = ',i3)") nx
              write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')")
              write (report_unit, fmt="('nx should have only 2, 3, and/or 5 as prime factors.')")
              write (report_unit, fmt="('################# WARNING #######################')")
              write (report_unit, *) 
           end if
        else
           write (report_unit, fmt="('nmesh=(2*ntgrid+1)*2*nlambda*negrid*ntheta0*naky*nspec')")
           nmesh = (2*ntgrid+1)*2*nlambda*negrid*ntheta0*naky*nspec
           write (report_unit, fmt="('Number of meshpoints:    ',i14)") nmesh
        end if

        write (report_unit, fmt="(T12,' ntgrid=',i12)") ntgrid
        write (report_unit, fmt="(T9,'2*ntgrid+1=',i12)") 2*ntgrid+1
        write (report_unit, fmt="(T12,'nlambda=',i12)") nlambda
        write (report_unit, fmt="(T12,' negrid=',i12)") negrid
        write (report_unit, fmt="(T12,'ntheta0=',i12)") ntheta0
        write (report_unit, fmt="(T12,'     nx=',i12)") nx
        write (report_unit, fmt="(T12,'   naky=',i12)") naky
        write (report_unit, fmt="(T12,'     ny=',i12)") ny
        write (report_unit, fmt="(T12,'  nspec=',i12,/)") nspec

        call nprocs (nmesh)
        if (nonlin) then
           write (report_unit, fmt="(/'Nonlinear run => consider #proc sweetspots for xxf+yxf objects!')") 
           call nprocs_xxf(nmesh)
           call nprocs_yxf(nmesh)
        endif

        if(use_le_layout) then
          write (report_unit, fmt="(/'Collisions using le_lo')") 
          call nprocs_le(nmesh)
        else
          select case (collision_model_switch)
          case (collision_model_full)
             write (report_unit, fmt="(/'Collisions using lz_lo')") 
             call nprocs_lz(nmesh)
             write (report_unit, fmt="(/'Collisions using e_lo')") 
             call nprocs_e(nmesh)
          case(collision_model_lorentz,collision_model_lorentz_test)
             write (report_unit, fmt="(/'Collisions using lz_lo')") 
             call nprocs_lz(nmesh)
          case (collision_model_ediffuse)
             write (report_unit, fmt="(/'Collisions using e_lo')") 
             call nprocs_e(nmesh)
          end select
       end if
     endif

     write (report_unit, *) 
     write (report_unit, fmt="('------------------------------------------------------------')")
     write (report_unit, *) 

     call check_species(report_unit,beta,tite,alne,dbetadrho_spec)

     write (report_unit, *) 
     write (report_unit, fmt="('------------------------------------------------------------')")
     write (report_unit, *) 
     call check_theta_grid(report_unit,alne,dbetadrho_spec)

     write (report_unit, fmt="(/'------------------------------------------------------------'/)")

     if (coll_on) then 
        call check_collisions(report_unit) 
     else
       write (report_unit, fmt="('All collisionality parameters (vnewk) are zero.')")
       write (report_unit, fmt="('No collision operator will be used.')")
     end if
       
    write (report_unit, *) 
    write (report_unit, fmt="('------------------------------------------------------------')")
    write (report_unit, *) 
    call check_fields(report_unit) 

    
    write (report_unit, *) 
    write (report_unit, fmt="('------------------------------------------------------------')")
    write (report_unit, *) 

    call check_init_g(report_unit)
       
    write (report_unit, *) 
    write (report_unit, fmt="('------------------------------------------------------------')")
    write (report_unit, *) 

    call check_nonlinear_terms(report_unit,delt_adj)
    if (nonlin) then
    else
       write (report_unit, *) 
       write (report_unit, fmt="('This is a linear calculation.')")
       write (report_unit, *)        
       write (report_unit, fmt="('The time step (code_dt) = ',e10.4)") code_dt
       write (report_unit, fmt="('The time step (user_dt) = ',e10.4)") user_dt
       write (report_unit, fmt="('The maximum number of time steps is ',i7)") nstep
       if (exit_when_converged) then
          write (report_unit, fmt="('When the frequencies for each k have converged, the run will stop.')")
          write (report_unit, fmt="('The convergence has to be better than one part in ',e10.4)") 1./omegatol
       end if

       if (wstar_units) then
          write (report_unit, *) 
          write (report_unit, fmt="('The timestep for each ky is scaled by a factor of 1/ky.')")
       end if
    end if

    write (report_unit, *) 
    write (report_unit, fmt="('Data will be written to ',a,' every ',i4,' timesteps.')") trim(run_name)//'.out.nc', nwrite
    write (report_unit, *) 

    if (dump_fields_periodically) then
       write (report_unit, *) 
       write (report_unit, fmt="('Data will be written to dump.fields.t=* every ',i4,' timesteps.')") 10*nmovie
       write (report_unit, *) 
    end if

    if (make_movie) then
       write (report_unit, *) 
       write (report_unit, fmt="('Movie data will be written to runname.movie.nc every ',i4,' timesteps.')") nmovie
       write (report_unit, *) 
    end if

    if (save_for_restart .and. nsave > 0) then
       write (report_unit, *) 
       write (report_unit, fmt="('Restart data will be written every ',i4,' timesteps.')") nsave
       write (report_unit, *) 
    end if
    
    write (report_unit, *) 
    write (report_unit, fmt="('------------------------------------------------------------')")
    write (report_unit, *) 

    call check_kt_grids(report_unit)

    write (report_unit, *) 
    write (report_unit, fmt="('------------------------------------------------------------')")
    write (report_unit, *) 

    call check_dist_fn(report_unit)

    write (report_unit, *) 
    write (report_unit, fmt="('------------------------------------------------------------')")
    write (report_unit, *) 

    call check_antenna(report_unit)

    call check_hyper(report_unit)

! BD broke this accidentally and will fix it.  July 22, 2011
!
!    write (report_unit, fmt="(/'------------------------------------------------------------'/)")
!
!    call check_le_grids(report_unit,le_ok)   

    write (report_unit, fmt="(/'------------------------------------------------------------'/)")
    call check_run_parameters(report_unit)

! diagnostic controls:
    call check_gs2_diagnostics(report_unit)

    call close_output_file (report_unit)

  end subroutine report
 
  subroutine wsweetspots(sym,sdim,nfac,facs,npmax,LUN)
! writes out sweetspot core counts for a layout
!   sym: character string label for each dimension 
!   sdim:  size of each dimension 
!   nfac:  #factors of each sdim
!   facs(i,j):  ith factor of jth dimension
 
   character(3), dimension(:), intent(in):: sym
   integer, dimension(:), intent(in):: sdim, nfac
   integer, dimension(:,:), intent(in):: facs
   integer, intent(in):: npmax
   integer, optional:: LUN
   integer :: lout=6, npe, i, j, lcores, nfacs
   if (present(LUN)) lout=LUN
   nfacs=size(nfac)
   lcores=1
   write (lout, fmt="('  npe = ',i8,'  (',a,')')") 1,trim(sym(1))
   do i=1,nfacs
      do j=2,nfac(i)
         npe=facs(j,i)*lcores
         if (npe .gt. npmax) exit
         write (lout, fmt="('  npe = ',i8,'  (',a,')')") npe,trim(sym(i))
      enddo
      lcores=lcores*sdim(i)
   enddo
  end subroutine wsweetspots

  subroutine nprocs_xxf(nmesh)
    use nonlinear_terms, only : nonlin
    use species, only : nspec
    use kt_grids, only: gridopt_switch, gridopt_single, gridopt_range, gridopt_specified, gridopt_box
    use kt_grids, only: naky, ntheta0
    use le_grids, only: negrid, nlambda
    use theta_grid, only: ntgrid
    use gs2_layouts, only: layout
    implicit none
    real :: fac
    integer, intent (in) :: nmesh
    integer :: nefacs, nlfacs, nkxfacs, nkyfacs, nsgfacs, nspfacs, ntgfacs
    integer, dimension(:,:), allocatable :: facs
    integer :: npe
    real :: time
    integer :: maxfacs
    integer, allocatable, dimension(:):: spfacs, efacs, lfacs, sgfacs, tgfacs, kyfacs
    integer, dimension(6) :: nfac, sdim
    character(3), dimension(6):: sym

    if (.not.nonlin) return
    write (report_unit, fmt="('xxf sweetspot #proc up to:',i8)") npmax
    maxfacs=max(nspec,negrid,nlambda,2,2*ntgrid+1,naky)/2+1
    allocate (spfacs(maxfacs),efacs(maxfacs),lfacs(maxfacs),sgfacs(maxfacs),tgfacs(maxfacs),kyfacs(maxfacs),facs(maxfacs,6))
    call factors (nspec, nspfacs, spfacs)
    call factors (negrid, nefacs, efacs)
    call factors (nlambda, nlfacs, lfacs)
    call factors (2, nsgfacs, sgfacs)
    call factors (2*ntgrid+1, ntgfacs, tgfacs)
    call factors (naky, nkyfacs, kyfacs)

    select case (layout)
       case ('lexys','lxyes','lyxes','yxles','xyles')
          sym = (/ "s  ","e  ","l  ","sgn","tg ","y  " /)
          sdim = (/ nspec, negrid, nlambda, 2, 2*ntgrid+1, naky /)
          nfac= (/ nspfacs, nefacs, nlfacs, nsgfacs, ntgfacs, nkyfacs /)
          facs(:,1)=spfacs; facs(:,2)=efacs; facs(:,3)=lfacs
          facs(:,4)=sgfacs; facs(:,5)=tgfacs; facs(:,6)=nkyfacs
          call wsweetspots(sym,sdim,nfac,facs,npmax,LUN=report_unit)
       case ('yxels')
          sym = (/ "s  ","l  ","e  ","sgn","tg ","y  " /)
          sdim = (/ nspec, nlambda, negrid, 2, 2*ntgrid+1, naky /)
          nfac= (/ nspfacs, nlfacs, nefacs, nsgfacs, ntgfacs, nkyfacs /)
          facs(:,1)=spfacs; facs(:,2)=lfacs; facs(:,3)=efacs
          facs(:,4)=sgfacs; facs(:,5)=tgfacs; facs(:,6)=nkyfacs
          call wsweetspots(sym,sdim,nfac,facs,npmax,LUN=report_unit)
       end select
       deallocate (facs,spfacs,efacs,lfacs,sgfacs,tgfacs,kyfacs)
  end subroutine nprocs_xxf

  subroutine nprocs_yxf(nmesh)
    use nonlinear_terms, only : nonlin
    use species, only : nspec
    use kt_grids, only: gridopt_switch, gridopt_single, gridopt_range, gridopt_specified, gridopt_box
    use kt_grids, only: naky, ntheta0, nx
    use le_grids, only: negrid, nlambda
    use theta_grid, only: ntgrid

    use gs2_layouts, only: layout
    implicit none
    real :: fac
    integer, intent (in) :: nmesh
    integer :: nefacs, nlfacs, nkxfacs, nkyfacs, nsgfacs, nspfacs, ntgfacs
    integer, dimension(:,:), allocatable :: facs
    integer :: npe
    real :: time
    integer :: maxfacs
    integer, allocatable, dimension(:):: spfacs, efacs, lfacs, sgfacs, tgfacs, kxfacs
    integer, dimension(6) :: nfac, sdim
    character(3), dimension(6):: sym

    if (.not.nonlin) return
    write (report_unit, fmt="('yxf sweetspot #proc up to:',i8)") npmax 
    maxfacs=max(nspec,negrid,nlambda,2,2*ntgrid+1,naky)/2+1
    allocate (spfacs(maxfacs),efacs(maxfacs),lfacs(maxfacs),sgfacs(maxfacs),tgfacs(maxfacs),kxfacs(maxfacs),facs(maxfacs,6))
    call factors (nspec, nspfacs, spfacs)
    call factors (negrid, nefacs, efacs)
    call factors (nlambda, nlfacs, lfacs)
    call factors (2, nsgfacs, sgfacs)
    call factors (2*ntgrid+1, ntgfacs, tgfacs)
    call factors (nx, nkxfacs, kxfacs)

    select case (layout)
       case ('lexys','lxyes','lyxes','yxles','xyles')
          sym = (/ "s  ","e  ","l  ","sgn","tg ","x  " /)
          sdim = (/ nspec, negrid, nlambda, 2, 2*ntgrid+1, nx /)
          nfac= (/ nspfacs, nefacs, nlfacs, nsgfacs, ntgfacs, nkxfacs /)
          facs(:,1)=spfacs; facs(:,2)=efacs; facs(:,3)=lfacs
          facs(:,4)=sgfacs; facs(:,5)=tgfacs; facs(:,6)=kxfacs
          call wsweetspots(sym,sdim,nfac,facs,npmax,LUN=report_unit)
       case ('yxels')
          sym = (/ "s  ","l  ","e  ","sgn","tg ","x  " /)
          sdim = (/ nspec, nlambda, negrid, 2, 2*ntgrid+1, nx /)
          nfac= (/ nspfacs, nlfacs, nefacs, nsgfacs, ntgfacs, nkxfacs /)
          facs(:,1)=spfacs; facs(:,2)=lfacs; facs(:,3)=efacs
          facs(:,4)=sgfacs; facs(:,5)=tgfacs; facs(:,6)=nkxfacs
          call wsweetspots(sym,sdim,nfac,facs,npmax,LUN=report_unit)
     end select
     deallocate (facs,spfacs,efacs,lfacs,sgfacs,tgfacs,kxfacs)
  end subroutine nprocs_yxf

  subroutine nprocs_e(nmesh)
    use species, only : nspec
    use kt_grids, only: gridopt_switch, gridopt_single, gridopt_range, gridopt_specified, gridopt_box
    use kt_grids, only: naky, ntheta0
    use le_grids, only: negrid, nlambda
    use theta_grid, only: ntgrid

    use gs2_layouts, only: layout
    implicit none
    real :: fac
    integer, intent (in) :: nmesh
    integer :: nefacs, nlfacs, nkxfacs, nkyfacs, nsgfacs, nspfacs, ntgfacs
    integer, dimension(:,:), allocatable :: facs
    integer :: npe
    real :: time
    integer :: maxfacs
    integer, allocatable, dimension(:):: spfacs, lfacs, sgfacs, tgfacs, kxfacs, kyfacs
    integer, dimension(6) :: nfac, sdim
    character(3), dimension(6):: sym

    write (report_unit, fmt="('#proc sweetspots for e_lo, up to:',i8)") npmax
    maxfacs=max(nspec,nlambda,2,2*ntgrid+1,ntheta0,naky)/2+1
    allocate (spfacs(maxfacs),lfacs(maxfacs),sgfacs(maxfacs),tgfacs(maxfacs),kxfacs(maxfacs),kyfacs(maxfacs),facs(maxfacs,6))
    call factors (nspec, nspfacs, spfacs)
    call factors (nlambda, nlfacs, lfacs)
    call factors (2, nsgfacs, sgfacs)
    call factors (2*ntgrid+1, ntgfacs, tgfacs)
    call factors (ntheta0, nkxfacs, kxfacs)
    call factors (naky, nkyfacs, kyfacs)

    select case (layout)
       case ('lexys','lxyes')
          sym = (/ "s  ","y  ","x  ","l  ","sgn","tg " /)
          sdim = (/ nspec, naky, ntheta0, nlambda, 2, 2*ntgrid+1 /)
          nfac= (/ nspfacs, nkyfacs, nkxfacs, nlfacs, nsgfacs, ntgfacs /)
          facs(:,1)=spfacs; facs(:,2)=kyfacs; facs(:,3)=kxfacs
          facs(:,4)=lfacs; facs(:,5)=sgfacs; facs(:,6)=tgfacs
          call wsweetspots(sym,sdim,nfac,facs,npmax,LUN=report_unit)
       case ('lyxes')
          sym = (/ "s  ","x  ","y  ","l  ","sgn","tg " /)
          sdim = (/ nspec, ntheta0, naky, nlambda, 2, 2*ntgrid+1 /)
          nfac= (/ nspfacs, nkxfacs, nkyfacs, nlfacs, nsgfacs, ntgfacs /)
          facs(:,1)=spfacs; facs(:,2)=kxfacs; facs(:,3)=kyfacs
          facs(:,4)=lfacs; facs(:,5)=sgfacs; facs(:,6)=tgfacs
          call wsweetspots(sym,sdim,nfac,facs,npmax,LUN=report_unit)
       case ('xyles')
          sym = (/ "s  ","l  ","y  ","x  ","sgn","tg " /)
          sdim = (/ nspec, nlambda, naky, ntheta0, 2, 2*ntgrid+1 /)
          nfac= (/ nspfacs, nlfacs, nkyfacs, nkxfacs, nsgfacs, ntgfacs /)
          facs(:,1)=spfacs; facs(:,2)=lfacs; facs(:,3)=kyfacs
          facs(:,4)=kxfacs; facs(:,5)=sgfacs; facs(:,6)=tgfacs
          call wsweetspots(sym,sdim,nfac,facs,npmax,LUN=report_unit)
       case ('yxles','yxels')
          sym = (/ "s  ","l  ","x  ","y  ","sgn","tg " /)
          sdim = (/ nspec, nlambda, ntheta0, naky, 2, 2*ntgrid+1 /)
          nfac= (/ nspfacs, nlfacs, nkxfacs, nkyfacs, nsgfacs, ntgfacs /)
          facs(:,1)=spfacs; facs(:,2)=lfacs; facs(:,3)=kxfacs
          facs(:,4)=kyfacs; facs(:,5)=sgfacs; facs(:,6)=tgfacs
          call wsweetspots(sym,sdim,nfac,facs,npmax,LUN=report_unit)
     end select
     deallocate (facs,spfacs,lfacs,sgfacs,tgfacs,kxfacs,kyfacs)
  end subroutine nprocs_e


  subroutine nprocs_lz(nmesh)
    use species, only : nspec
    use kt_grids, only: gridopt_switch, gridopt_single, gridopt_range, gridopt_specified, gridopt_box
    use kt_grids, only: naky, ntheta0
    use le_grids, only: negrid, nlambda
    use theta_grid, only: ntgrid

    use gs2_layouts, only: layout
    implicit none
    real :: fac
    integer, intent (in) :: nmesh
    integer :: nefacs, nlfacs, nkxfacs, nkyfacs, nsgfacs, nspfacs, ntgfacs
    integer, dimension(:,:), allocatable :: facs
    integer :: npe
    real :: time
    integer :: maxfacs
    integer, allocatable, dimension(:):: spfacs, efacs, sgfacs, tgfacs, kxfacs, kyfacs
    integer, dimension(5) :: nfac, sdim
    character(3), dimension(5):: sym

    write (report_unit, fmt="('#proc sweetspots for lz_lo, up to:',i8)") npmax

    maxfacs=max(nspec,negrid,2*ntgrid+1,ntheta0,naky)/2+1
    allocate (spfacs(maxfacs),efacs(maxfacs),tgfacs(maxfacs),kxfacs(maxfacs),kyfacs(maxfacs),facs(maxfacs,6))
    call factors (nspec, nspfacs, spfacs)
    call factors (negrid, nefacs, efacs)
    call factors (2*ntgrid+1, ntgfacs, tgfacs)
    call factors (ntheta0, nkxfacs, kxfacs)
    call factors (naky, nkyfacs, kyfacs)
 
    select case (layout)
       case ('lexys')
          sym = (/ "s  ","y  ","x  ","e  ","tg " /)
          sdim = (/ nspec, naky, ntheta0, negrid, 2*ntgrid+1 /)
          nfac= (/ nspfacs, nkyfacs, nkxfacs, nefacs, ntgfacs /)
          facs(:,1)=spfacs; facs(:,2)=kyfacs; facs(:,3)=kxfacs
          facs(:,4)=efacs; facs(:,5)=tgfacs
          call wsweetspots(sym,sdim,nfac,facs,npmax,LUN=report_unit)
       case ('lyxes','yxles','yxels')
          sym = (/ "s  ","e  ","x  ","y  ","tg " /)
          sdim = (/ nspec, negrid, ntheta0, naky, 2*ntgrid+1 /)
          nfac= (/ nspfacs, nefacs, nkxfacs, nkyfacs, ntgfacs /)
          facs(:,1)=spfacs; facs(:,2)=efacs; facs(:,3)=kxfacs
          facs(:,4)=kyfacs; facs(:,5)=tgfacs
          call wsweetspots(sym,sdim,nfac,facs,npmax,LUN=report_unit)
       case ('xyles','lxyes')
          sym = (/ "s  ","e  ","y  ","x  ","tg " /)
          sdim = (/ nspec, negrid, naky, ntheta0, 2*ntgrid+1 /)
          nfac= (/ nspfacs, nefacs, nkyfacs, nkxfacs, ntgfacs /)
          facs(:,1)=spfacs; facs(:,2)=efacs; facs(:,3)=kyfacs
          facs(:,4)=kxfacs; facs(:,5)=tgfacs
          call wsweetspots(sym,sdim,nfac,facs,npmax,LUN=report_unit)
    end select
    deallocate (spfacs,efacs,tgfacs,kxfacs,kyfacs,facs)
  end subroutine nprocs_lz


  subroutine nprocs_le(nmesh)
    use species, only : nspec
    use kt_grids, only: gridopt_switch, gridopt_single, gridopt_range, gridopt_specified, gridopt_box
    use kt_grids, only: naky, ntheta0
    use le_grids, only: negrid, nlambda
    use theta_grid, only: ntgrid

    use gs2_layouts, only: layout
    implicit none
    real :: fac
    integer, intent (in) :: nmesh
    integer :: nefacs, nlfacs, nkxfacs, nkyfacs, nsgfacs, nspfacs, ntgfacs
    integer, dimension(:,:), allocatable :: facs
    integer :: npe
    real :: time
    integer :: maxfacs
    integer, allocatable, dimension(:):: spfacs, efacs, sgfacs, tgfacs, kxfacs, kyfacs
    integer, dimension(4) :: nfac, sdim
    character(3), dimension(4):: sym

    write (report_unit, fmt="('#proc sweetspots for le_lo, up to:',i8)") npmax

    maxfacs=max(nspec,2*ntgrid+1,ntheta0,naky)/2+1
    allocate (spfacs(maxfacs),tgfacs(maxfacs),kxfacs(maxfacs),kyfacs(maxfacs),facs(maxfacs,4))
    call factors (nspec, nspfacs, spfacs)
    call factors (2*ntgrid+1, ntgfacs, tgfacs)
    call factors (ntheta0, nkxfacs, kxfacs)
    call factors (naky, nkyfacs, kyfacs)

    select case (layout)
       case ('lexys','xyles','lxyes')
          sym = (/ "s  ","y  ","x  ","tg " /)
          sdim = (/ nspec, naky, ntheta0, 2*ntgrid+1 /)
          nfac= (/ nspfacs, nkyfacs, nkxfacs, ntgfacs /)
          facs(:,1)=spfacs; facs(:,2)=kyfacs
          facs(:,3)=kxfacs; facs(:,4)=tgfacs
          call wsweetspots(sym,sdim,nfac,facs,npmax,LUN=report_unit)
       case ('lyxes','yxles','yxels')
          sym = (/ "s  ","x  ","y  ","tg " /)
          sdim = (/ nspec, ntheta0, naky, 2*ntgrid+1 /)
          nfac= (/ nspfacs, nkxfacs, nkyfacs, ntgfacs /)
          facs(:,1)=spfacs; facs(:,2)=kxfacs
          facs(:,3)=kyfacs; facs(:,4)=tgfacs
          call wsweetspots(sym,sdim,nfac,facs,npmax,LUN=report_unit)
    end select
    deallocate (spfacs,tgfacs,kxfacs,kyfacs,facs)
  end subroutine nprocs_le


  subroutine nprocs (nmesh)
    use nonlinear_terms, only : nonlin
    use species, only : nspec
    use kt_grids, only: gridopt_switch, gridopt_single, gridopt_range, gridopt_specified, gridopt_box
    use kt_grids, only: naky, ntheta0, nx, ny
    use le_grids, only: negrid, nlambda
    use theta_grid, only: ntgrid
    use gs2_layouts, only: layout, init_x_transform_layouts, init_y_transform_layouts
    implicit none
    real :: fac
    integer, intent (in) :: nmesh
    integer :: nefacs, nlfacs, nkyfacs, nkxfacs, nspfacs, nkxkyfacs
    integer, dimension(:,:), allocatable :: facs
    integer :: npe, checknpe, i, j
    logical :: onlyxoryfac
    real :: time

    write (report_unit, fmt="('Layout = ',a5,/)") layout 
    write (report_unit, fmt="('Recommended #proc up to:',i8)") npmax 
    if (nonlin) then

       write (report_unit, *) 
       write (report_unit, fmt="('----------------------------------------------------------------------------------------------------------------------------------------------')")
       write (report_unit, fmt="('|   Number of    | Dimension | Percentage data moved | Recommend unbalanced_xxf | xxf unbalanced | Recommend unbalanced_yxf | yxf unbalanced |')")
       write (report_unit, fmt="('|   processes    |   split   | by MPI from xxf->yxf  |    T = true F = false    |   amount (%)   |    T = true F = false    |   amount (%)   |')")
       write (report_unit, fmt="('----------------------------------------------------------------------------------------------------------------------------------------------')")
       
       
       call init_x_transform_layouts(ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx)
       call init_y_transform_layouts(ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny)

       select case (layout)
       case ('lexys')

!          write (report_unit, fmt="('Recommended numbers of processors, time on T3E')") 
          allocate (facs(max(nspec,naky,ntheta0)/2+1,5))
          call factors (nspec, nspfacs, facs(:,1))
          call factors (naky, nkyfacs, facs(:,2))
          call factors (ntheta0, nkxfacs, facs(:,3))
          call factors (negrid, nefacs, facs(:,4))
          call factors (nlambda, nlfacs, facs(:,5))
          do i=1,nspfacs
             npe = facs(i,1)
             if (npe .gt. npmax) exit
             call report_idle_procs(npe, 's')
          end do
          do i=2,nkyfacs
             npe = facs(i,2)*nspec
             if (npe .gt. npmax) exit
             call report_idle_procs(npe, 'y')
          end do
          do i=2,nkxfacs
             npe = facs(i,3)*naky*nspec
             if (npe .gt. npmax) exit
             call report_idle_procs(npe, 'x')
          end do
          do i=2,nefacs
             npe = facs(i,4)*ntheta0*naky*nspec
             if (npe .gt. npmax) exit
             call report_idle_procs(npe, 'e')
          end do
          do i=2,nlfacs
             npe = facs(i,5)*negrid*ntheta0*naky*nspec
             if (npe .gt. npmax) exit
             call report_idle_procs(npe, 'l')
          end do
          deallocate (facs)

       case ('lxyes')

!          write (report_unit, fmt="('Recommended numbers of processors, time on SP2')") 
          allocate (facs(max(nspec,negrid,naky,ntheta0)/2+1,5))
          call factors (nspec, nspfacs, facs(:,1))
          call factors (negrid, nefacs, facs(:,2))
          call factors (naky, nkyfacs, facs(:,3))
          call factors (ntheta0, nkxfacs, facs(:,4))
          call factors (nlambda, nlfacs, facs(:,5))
!          faclin = 3.5*(real(nmesh))**1.1/1.e7/5.
          fac = 3.5*(real(nmesh))**1.1/1.e7
          do i=1,nspfacs
             npe = facs(i,1)
             if (npe .gt. npmax) exit
             call report_idle_procs(npe, 's')
          end do
          do i=2,nefacs
             npe = facs(i,2)*nspec
             if (npe .gt. npmax) exit
             call report_idle_procs(npe, 'e')
          end do
          do i=2,nkyfacs
             npe = facs(i,3)*negrid*nspec
             if (npe .gt. npmax) exit
             call report_idle_procs(npe, 'y')
          end do
          do i=2,nkxfacs
             npe = facs(i,4)*naky*negrid*nspec
             if (npe .gt. npmax) exit
             call report_idle_procs(npe, 'x')
          end do
          do i=2,nlfacs
             npe = facs(i,5)*naky*ntheta0*negrid*nspec 
             if (npe .gt. npmax) exit
             call report_idle_procs(npe, 'l')
          end do
          deallocate (facs)

       case ('yxels')

          allocate (facs(max(naky,ntheta0,nspec,negrid,nlambda,ntheta0*naky)/2+1,6))
          call factors (nspec, nspfacs, facs(:,1))
          call factors (nlambda, nlfacs, facs(:,2))
          call factors (negrid, nefacs, facs(:,3))
          call factors (ntheta0, nkxfacs, facs(:,4))
          call factors (naky, nkyfacs, facs(:,5))
          call factors (naky*ntheta0, nkxkyfacs, facs(:,6))
          fac = 3.5*(real(nmesh))**1.1/1.e7
          do i=1,nspfacs
             npe = facs(i,1)
             if (npe .gt. npmax) exit
             call report_idle_procs(npe, 's')
          end do
          do i=2,nlfacs
             npe = facs(i,2)*nspec
             if (npe .gt. npmax) exit
             call report_idle_procs(npe, 'l')
          end do
          do i=2,nefacs
             npe = facs(i,3)*nlambda*nspec
             if (npe .gt. npmax) exit
             call report_idle_procs(npe, 'e')
          end do
          do i=2,nkxkyfacs
             npe = facs(i,6)*nlambda*negrid*nspec
             if (npe .gt. npmax) exit
! Check whether this process count would have been generated by the 
! plain x or y factorisation or if it is new from the combined x*y
! functionality             
             onlyxoryfac = .false.
             do j=2,nkxfacs
                checknpe = facs(j,4)*negrid*nlambda*nspec
                if (npe .eq. checknpe) then
                   onlyxoryfac = .true.
                   exit
                end if
             end do
             do j=2,nkyfacs
                checknpe = facs(j,5)*ntheta0*negrid*nlambda*nspec
                if (npe .gt. checknpe) then
                   onlyxoryfac = .true.
                   exit
                end if
             end do
             call report_idle_procs(npe, 'x*y', onlyxoryfac)
          end do
!          do i=2,nkxfacs
!             npe = facs(i,4)*negrid*nlambda*nspec
!             if (npe .gt. npmax) exit
!             call report_idle_procs(npe, 'x')
!          end do
!          do i=2,nkyfacs
!             npe = facs(i,5)*ntheta0*negrid*nlambda*nspec
!             if (npe .gt. npmax) exit
!             call report_idle_procs(npe, 'y')
!          end do
          deallocate (facs)

       case ('yxles')

          allocate (facs(max(ntheta0,naky,nspec,negrid,nlambda,ntheta0*naky)/2+1,6))
          call factors (nspec, nspfacs, facs(:,1))
          call factors (negrid, nefacs, facs(:,2))
          call factors (nlambda, nlfacs, facs(:,3))
          call factors (ntheta0, nkxfacs, facs(:,4))
          call factors (naky, nkyfacs, facs(:,5))
          call factors (naky*ntheta0, nkxkyfacs, facs(:,6))
          fac = 3.5*(real(nmesh))**1.1/1.e7
          do i=1,nspfacs
             npe = facs(i,1)
             if (npe .gt. npmax) exit
             call report_idle_procs(npe, 's')
          end do
          do i=2,nefacs
             npe = facs(i,2)*nspec
             if (npe .gt. npmax) exit
             call report_idle_procs(npe, 'e')
          end do
          do i=2,nlfacs
             npe = facs(i,3)*negrid*nspec
             if (npe .gt. npmax) exit
             call report_idle_procs(npe, 'l')
          end do
          do i=2,nkxkyfacs
             npe = facs(i,6)*nlambda*negrid*nspec
             if (npe .gt. npmax) exit
! Check whether this process count would have been generated by the 
! plain x or y factorisation or if it is new from the combined x*y
! functionality  
             onlyxoryfac = .false.
             do j=2,nkxfacs
                checknpe = facs(j,4)*nlambda*negrid*nspec
                if (npe .eq. checknpe) then
                   onlyxoryfac = .true.
                   exit
                end if
             end do
             do j=2,nkyfacs
                checknpe = facs(j,5)*ntheta0*nlambda*negrid*nspec
                if (npe .eq. checknpe) then
                   onlyxoryfac = .true.
                   exit
                end if
             end do
             call report_idle_procs(npe, 'x*y', onlyxoryfac)
          end do
!          do i=2,nkxfacs
!             npe = facs(i,4)*nlambda*negrid*nspec
!             if (npe .gt. npmax) exit
!             call report_idle_procs(npe, 'x')
!          end do
!          do i=2,nkyfacs
!             npe = facs(i,5)*ntheta0*nlambda*negrid*nspec
!             if (npe .gt. npmax) exit
!             call report_idle_procs(npe, 'y')
!          end do
          deallocate (facs)

       case ('xyles')
!CMR, 11/11/2009: add processor recommendations for xyles layout
!            NB added recommendations that also parallelise in y and x, 
!                          which may be unwise!

          allocate (facs(max(nspec,negrid,nlambda,naky,ntheta0,ntheta0*naky)/2+1,6))
          call factors (nspec, nspfacs, facs(:,1))
          call factors (negrid, nefacs, facs(:,2))
          call factors (nlambda, nlfacs, facs(:,3))
          call factors (naky, nkyfacs, facs(:,4))
          call factors (ntheta0, nkxfacs, facs(:,5))
          call factors (naky*ntheta0, nkxkyfacs, facs(:,6))
          fac = 3.5*(real(nmesh))**1.1/1.e7
          do i=1,nspfacs
             npe = facs(i,1)
             if (npe .gt. npmax) exit
             call report_idle_procs(npe, 's')
          end do
          do i=2,nefacs
             npe = facs(i,2)*nspec
             if (npe .gt. npmax) exit
             call report_idle_procs(npe, 'e')
          end do
          do i=2,nlfacs
             npe = facs(i,3)*negrid*nspec
             if (npe .gt. npmax) exit
             call report_idle_procs(npe, 'l')
          end do
          do i=2,nkxkyfacs
             npe = facs(i,6)*nlambda*negrid*nspec
             if (npe .gt. npmax) exit
! Check whether this process count would have been generated by the 
! plain x or y factorisation or if it is new from the combined x*y
! functionality  
             onlyxoryfac = .false.
             do j=2,nkyfacs
                checknpe = facs(j,4)*nlambda*negrid*nspec
                if (npe .eq. checknpe) then
                   onlyxoryfac = .true.
                   exit
                end if
             end do
             do j=2,nkxfacs
                checknpe = facs(j,5)*naky*nlambda*negrid*nspec
                if (npe .eq. checknpe) then
                   onlyxoryfac = .true.
                   exit
                end if
             end do
             call report_idle_procs(npe, 'x*y', onlyxoryfac)
          end do
!          do i=2,nkyfacs
!             npe = facs(i,4)*nlambda*negrid*nspec
!             if (npe .gt. npmax) exit
!             call report_idle_procs(npe, 'y')
!          end do
!          do i=2,nkxfacs
!             npe = facs(i,5)*naky*nlambda*negrid*nspec
!             if (npe .gt. npmax) exit
!             call report_idle_procs(npe, 'x')
!          end do
          deallocate (facs)

       end select

       write (report_unit, fmt="('----------------------------------------------------------------------------------------------------------------------------------------------')")
       write (report_unit, *)
       write (report_unit, fmt="('(*) denotes process counts that are from factors of the combined kx*ky index rather than the ordered kx or ky indices separately.')")
       write (report_unit, fmt="('To use the unbalanced functionality set unbalanced_xxf = .true. or unbalanced_yxf = .true. in the &layouts_knobs namelist in your GS2 ')")
       write (report_unit, fmt="('input file. You can also set the max_unbalanced_xxf and max_unbalanced_yxf flags in the same namelist in the input file to specify the')")
       write (report_unit, fmt="('maximum amount of computational imbalance allowed. These flags specify the maximum imbalance as 1 with no imbalance as 0, so to allow')")
       write (report_unit, fmt="('50% imbalanced on the xxf decomposition using the following flags in the input file: ')")
       write (report_unit, fmt="('                                                                                     unbalanced_xxf = .true. ')")
       write (report_unit, fmt="('                                                                                     max_unbalanced_xxf = 0.5 ')")
       write (report_unit, fmt="('And likewise for yxf.')")

    else
       write (report_unit, *) 
       write (report_unit, fmt="('Recommended numbers of processors:')")
       select case (gridopt_switch)

       case (gridopt_single)

          select case (layout)
          case ('lexys', 'yxles', 'lxyes', 'xyles')
             allocate (facs(max(negrid,nspec,nlambda)/2+1,3))
             call factors (nspec, nspfacs, facs(:,1))
             call factors (negrid, nefacs, facs(:,2))
             call factors (nlambda, nlfacs, facs(:,3))
             do i=1,nspfacs-1
                npe = facs(i,1)
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nefacs-1
                npe = facs(i,2)*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nlfacs
                npe = facs(i,3)*negrid*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             deallocate (facs)

          case ('yxels')
             allocate (facs(max(negrid,nspec,nlambda)/2+1,3))
             call factors (nspec, nspfacs, facs(:,1))
             call factors (nlambda, nlfacs, facs(:,2))
             call factors (negrid, nefacs, facs(:,3))
             do i=1,nspfacs-1
                npe = facs(i,1)
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nlfacs-1
                npe = facs(i,2)*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nefacs
                npe = facs(i,3)*nlambda*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             deallocate (facs)
          end select

       case (gridopt_range,gridopt_specified,gridopt_box)

          select case (layout)
          case ('lexys')

             allocate (facs(max(nspec,naky,ntheta0,negrid,nlambda)/2+1,5))
             call factors (nspec, nspfacs, facs(:,1))
             call factors (naky, nkyfacs, facs(:,2))
             call factors (ntheta0, nkxfacs, facs(:,3))
             call factors (negrid, nefacs, facs(:,4))
             call factors (nlambda, nlfacs, facs(:,5))
             do i=1,nspfacs-1
                npe = facs(i,1)
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nkyfacs-1
                npe = facs(i,2)*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nkxfacs-1
                npe = facs(i,3)*naky*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nefacs-1
                npe = facs(i,4)*ntheta0*naky*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nlfacs
                npe = facs(i,5)*negrid*ntheta0*naky*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             deallocate (facs)

          case ('lxyes')

             write (report_unit, *) 
             allocate (facs(max(nspec,negrid,naky,ntheta0)/2+1,4))
             call factors (nspec, nspfacs, facs(:,1))
             call factors (negrid, nefacs, facs(:,2))
             call factors (naky, nkyfacs, facs(:,3))
             call factors (ntheta0, nkxfacs, facs(:,4))
!             fac = 3.5*(real(nmesh))**1.1/1.e7/5.  ! okay for large runs, not small
             do i=1,nspfacs-1
                npe = facs(i,1)
!                if (nmesh/npe > ncut) write &
!                     & (report_unit, fmt="('  npe = ',i4,'    time = ',f8.2,'  seconds/time step')") &
!                     & npe, fac/npe**0.95
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nefacs-1
                npe = facs(i,2)*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nkyfacs-1
                npe = facs(i,3)*negrid*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nkxfacs
                npe = facs(i,4)*naky*negrid*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             deallocate (facs)

          case ('yxles')

             allocate (facs(max(nspec,naky,ntheta0,negrid,nlambda)/2+1,5))
             call factors (nspec, nspfacs, facs(:,1))
             call factors (negrid, nefacs, facs(:,2))
             call factors (nlambda, nlfacs, facs(:,3))
             call factors (ntheta0, nkxfacs, facs(:,4))
             call factors (naky, nkyfacs, facs(:,5))
             do i=1,nspfacs-1
                npe = facs(i,1)
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nefacs-1
                npe = facs(i,2)*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nlfacs-1
                npe = facs(i,3)*negrid*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nkxfacs-1
                npe = facs(i,4)*nlambda*negrid*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nkyfacs
                npe = facs(i,5)*ntheta0*nlambda*negrid*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             deallocate (facs)

          case ('xyles')

             allocate (facs(max(nspec,naky,ntheta0,negrid,nlambda)/2+1,5))
             call factors (nspec, nspfacs, facs(:,1))
             call factors (negrid, nefacs, facs(:,2))
             call factors (nlambda, nlfacs, facs(:,3))
             call factors (naky, nkyfacs, facs(:,4))
             call factors (ntheta0, nkxfacs, facs(:,5))
             do i=1,nspfacs-1
                npe = facs(i,1)
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i8)") npe
             end do
             do i=1,nefacs-1
                npe = facs(i,2)*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i8)") npe
             end do
             do i=1,nlfacs-1
                npe = facs(i,3)*negrid*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i8)") npe
             end do
             do i=1,nkyfacs-1
                npe = facs(i,4)*nlambda*negrid*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i8)") npe
             end do
             do i=1,nkxfacs
                npe = facs(i,5)*ntheta0*nlambda*negrid*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i8)") npe
             end do
             deallocate (facs)

          case ('yxels')

             allocate (facs(max(nspec,naky,ntheta0,negrid,nlambda)/2+1,5))
             call factors (nspec, nspfacs, facs(:,1))
             call factors (nlambda, nlfacs, facs(:,2))
             call factors (negrid, nefacs, facs(:,3))
             call factors (ntheta0, nkxfacs, facs(:,4))
             call factors (naky, nkyfacs, facs(:,5))
             do i=1,nspfacs-1
                npe = facs(i,1)
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nlfacs-1
                npe = facs(i,2)*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nefacs-1
                npe = facs(i,3)*nlambda*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nkxfacs-1
                npe = facs(i,4)*negrid*nlambda*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             do i=1,nkyfacs
                npe = facs(i,5)*ntheta0*negrid*nlambda*nspec
                if (nmesh/npe > ncut) write (report_unit, fmt="('  npe = ',i4)") npe
             end do
             deallocate (facs)

          end select
       end select
    end if

  end subroutine nprocs

  subroutine get_unbalanced_suggestions(npe, percentage_xxf_unbalanced_amount, percentage_yxf_unbalanced_amount, &
     use_unbalanced_xxf, use_unbalanced_yxf)
  !====================================================================================
  ! AJ June 2012
  ! This subroutine is used to return the unbalanced suggestions for a given process 
  ! count (npe).  The calculation is performed using the calculate_unbalanced_x 
  ! and calculate_unbalanced_y subroutines from gs2_layouts.  These return the 
  ! unbalanced decomposition size (difference between the suggested small and large
  ! block size for the decomposition) and from this a logical is set to recommend
  ! whether the unbalanced decomposition should be used or not.
  !====================================================================================
    use gs2_layouts, only : calculate_unbalanced_x, calculate_unbalanced_y
    implicit none

    integer, intent(in) :: npe
    integer, intent(out) :: percentage_xxf_unbalanced_amount, percentage_yxf_unbalanced_amount
    logical, intent(out) :: use_unbalanced_xxf, use_unbalanced_yxf

    real  :: xxf_unbalanced_amount, yxf_unbalanced_amount

    call calculate_unbalanced_x(npe, 0, xxf_unbalanced_amount)
    if(xxf_unbalanced_amount .gt. 0) then
       use_unbalanced_xxf = .true.
       percentage_xxf_unbalanced_amount = xxf_unbalanced_amount * 100
    else
       use_unbalanced_xxf = .false.  
       percentage_xxf_unbalanced_amount = 0
    end if
    call calculate_unbalanced_y(npe, 0, yxf_unbalanced_amount)
    if(yxf_unbalanced_amount .gt. 0) then
       use_unbalanced_yxf = .true.
       percentage_yxf_unbalanced_amount = yxf_unbalanced_amount * 100
    else
       use_unbalanced_yxf = .false. 
       percentage_yxf_unbalanced_amount = 0
    end if
    
  end subroutine get_unbalanced_suggestions

  subroutine get_idle_processes(npe, idle_percentage, use_unbalanced)
  !====================================================================================
  ! AJ June 2012
  ! This subroutine is used to return the idle processes from the xxf and yxf layouts.
  ! Idle processes sometimes occur, dependent on the process count used, because 
  ! the data domain does not evenly divide by the total number of processes available. 
  ! These can cause high communications overheads for the non-linear calculations
  ! for large numbers of processes so it is useful to print this data out in ingen
  ! to let users know which process counts this can happen at for a given input file
  ! We currently use an arbitrary cutoff of 10% difference in idle processes to 
  ! suggest that the unbalanced decomposition functionality should be used to 
  ! mitigate the impact of this difference. 
  !====================================================================================
    use gs2_layouts, only : calculate_idle_processes
    implicit none

    integer, intent(in) :: npe
    real, intent(out) :: idle_percentage
    logical, intent(out) :: use_unbalanced

    call calculate_idle_processes(npe, idle_percentage)
! This 0.1 represents the arbitrary 10% difference threshold discussed above.
    if(idle_percentage .gt. 0.1) then
       use_unbalanced = .true.
    else
       use_unbalanced = .false.
    end if
    idle_percentage = idle_percentage * 100

  end subroutine get_idle_processes
  
  subroutine report_idle_procs(npe, distchar, onlyxoryfac)
!====================================================================================
  ! AJ June 2012
  ! This subroutine wraps up the output functionality for the code that creates the 
  ! list of suggested process counts for the linear computation and checks whether 
  ! those process counts work well for the non-linear calculations as well.
  !
  ! The optional parameter onlyxoryfac (which is a logical) is used to highlight 
  ! process counts that have been generated using code which looks for the factors 
  ! of kx*ky rather than kx and ky separately.  This is only currently used for 
  ! layouts that start with x and y (for instance xyles or yxles).
  !====================================================================================
  
    implicit none

    integer, intent(in) :: npe
    character(*), intent(in) :: distchar
    logical, optional, intent(in) :: onlyxoryfac
    integer :: percentage_xxf_unbalanced_amount, percentage_yxf_unbalanced_amount
    logical :: use_unbalanced_xxf, use_unbalanced_yxf, use_unbalanced
    real :: idle_percentage

    call get_idle_processes(npe, idle_percentage, use_unbalanced)
    if(use_unbalanced) then
       call get_unbalanced_suggestions(npe, percentage_xxf_unbalanced_amount, percentage_yxf_unbalanced_amount, use_unbalanced_xxf, use_unbalanced_yxf)
       if(present(onlyxoryfac)) then
          if(onlyxoryfac) then
             write (report_unit, fmt="('|     ',i8,'   |    ',a3,'    |          ',i3,'          |            ',L1,'            |      ',i3,'       |            ',L1,'            |      ',i3,'       |')") & 
                  npe, distchar, INT(idle_percentage), use_unbalanced_xxf, percentage_xxf_unbalanced_amount, use_unbalanced_yxf, percentage_yxf_unbalanced_amount
          else
             write (report_unit, fmt="('|     ',i8,'(*)|    ',a3,'    |          ',i3,'          |            ',L1,'            |      ',i3,'       |            ',L1,'            |      ',i3,'       |')") & 
                  npe, distchar, INT(idle_percentage), use_unbalanced_xxf, percentage_xxf_unbalanced_amount, use_unbalanced_yxf, percentage_yxf_unbalanced_amount
          end if
       else
          write (report_unit, fmt="('|     ',i8,'   |    ',a3,'    |          ',i3,'          |            ',L1,'            |      ',i3,'       |            ',L1,'            |      ',i3,'       |')") & 
               npe, distchar, INT(idle_percentage), use_unbalanced_xxf, percentage_xxf_unbalanced_amount, use_unbalanced_yxf, percentage_yxf_unbalanced_amount
       end if
    else
       if(present(onlyxoryfac)) then
          if(onlyxoryfac) then
             write (report_unit, fmt="('|     ',i8,'   |    ',a3,'    |          ',i3,'          |                          |                |                          |                |')") & 
                  npe, distchar, INT(idle_percentage)
          else
             write (report_unit, fmt="('|     ',i8,'(*)|    ',a3,'    |          ',i3,'          |                          |                |                          |                |')") & 
                  npe, distchar, INT(idle_percentage)
          end if
       else
          write (report_unit, fmt="('|     ',i8,'   |    ',a3,'    |          ',i3,'          |                          |                |                          |                |')") & 
               npe, distchar, INT(idle_percentage)
       end if
    
    end if
    
  end subroutine report_idle_procs

  subroutine tell (a, b, c)
    
    integer :: j
    character (*), intent (in) :: a
    character (*), intent (in), optional :: b, c
    character (500) :: hack

    j = len_trim(a)

    if (present(b)) then
       j = max(j, len_trim(b))
       if (present(c)) j = max(j, len_trim(c))
    end if

    write (6, *)
    write (6, *)
    hack = repeat('*', j+2)
    write (6, 10) '  '//trim(hack)
    write (6, 10) '   '//trim(a)
    write (interactive_record, 10) trim(a)
    if (present(b)) then
       write (6, 10) '   '//trim(b)
       write (interactive_record, 10) trim(b)
       if (present(c)) then
          write (6, 10) '   '//trim(c)
          write (interactive_record, 10) trim(c)
       end if
    end if
    write (6, 10) '  '//trim(hack)
    write (6, *)
    write (6, *)

10 format (a)

  end subroutine tell

  subroutine text (a)
    
    character (*), intent (in), optional :: a

    if (present(a)) then
       write (6, 10) trim(a)
    else
       write (*, *)
    end if

10 format (a)

  end subroutine text

  subroutine run_number (sel, nbeta)
    
    integer, intent (out) :: sel
    integer, intent (in) :: nbeta
    
100 continue
    do 
       call text
       call text ('Lowest run number for this scan (0-999):')
       read (interactive_input, *, err = 100) sel
       if (sel + nbeta - 1 > 999) then
          write (6,*) 'The lowest run number for this scan must be less than ',1000-nbeta
          call try_again
       else if (sel < 0) then
          write (6,*) 'The lowest run number for this scan must be greater than zero.'
          call try_again
       else
          return
       end if
    end do
    
  end subroutine run_number
  
  subroutine try_again 
    
    write (6, *) '*****************'
    write (6, *) 'Please try again.'
    write (6, *) '*****************'

  end subroutine try_again

  subroutine get_choice (in, out)

    integer, intent (in) :: in
    integer, intent (out) :: out

    read (interactive_input, *, err=999) out
    if (out <= 0) out = 0
    if (out > in) out = 0

    return
999 out = 0
    
  end subroutine get_choice
  
  subroutine beta_range_low (x, a, x0)
    
    real, intent (out) :: x
    character (*), intent (in) :: a
    real, intent (in) :: x0
    
    do
100    continue
       call text
       if (trim(a) == 'le') then
          call text ('Lower limit of beta for this scan (must be > 0):')
       else if(trim(a) == 'lt') then
          call text ('Lower limit of beta for this scan (must be >= 0):')
       end if
       read (interactive_input, *, err=100) x
       if (trim(a) == 'le') then
          if (x <= x0) then
             call try_again
          else
             return
          end if
       else if (trim(a) == 'lt') then
          if (x < x0) then
             call try_again
          else
             return
          end if
       end if
    end do
    
  end subroutine beta_range_low

  subroutine beta_range_high (x, a, x0)
    
    real, intent (out) :: x
    character (*), intent (in) :: a
    real, intent (in) :: x0
    
    do
100    continue
       call text
       call text ('Upper limit of beta for this scan (must be > lower limit):')
       read (interactive_input, *, err=100) x
       if (trim(a) == 'le') then
          if (x <= x0) then
             call try_again
          else
             return
          end if
       end if
    end do
    
  end subroutine beta_range_high

  subroutine beta_prime_range_low (x)

    real, intent (out) :: x
    
    do
100    continue
       call tell ('Note: You will be asked to enter:         - d beta/d rho', &
            'Since d beta /d rho <= 0, you should enter a number >= 0')
       
       call text ('Weakest beta gradient in scan (zero or positive number):')
       read (interactive_input, *, err=100) x
       x = -x
       if (x > 0.) then
          call try_again
       else
          return
       end if
    end do

  end subroutine beta_prime_range_low
  
  subroutine beta_prime_range_high (x, x0)
    
    real, intent (out) :: x
    real, intent (in) :: x0
    
    do
100    continue
       call text
       call text ('Strongest gradient in scan: (positive number)')
       read (interactive_input, *, err=100) x
       x = -x
       if (x >= x0) then
          call try_again
       else
          return
       end if
    end do

  end subroutine beta_prime_range_high

  subroutine num_runs (n)

    integer, intent (out) :: n
    
    do 
100    continue
       call text
       call text ('How many runs should be done? (n > 1)')
       read (interactive_input, *, err=100) n
       if (n < 2) then
          call try_again
       else
          return
       end if
    end do
  end subroutine num_runs

end program ingen

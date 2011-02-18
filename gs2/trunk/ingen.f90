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
  logical:: debug=.false.

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
    use theta_grid, only: init_theta_grid, nbset, shat_real => shat
    use theta_grid_params, only: shat
    use constants, only: pi
    implicit none
    logical :: list, accelx=.false., accelv=.false.
    logical:: debug=.false.
    integer :: is
    call init_file_utils (list, name="template")
if (debug) write(6,*) 'get_namelists: called init_file_utils'
    ncut= 100000
    npmax=10000
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
if (debug) write(6,*) 'get_namelists: call init_theta_grid'
    call init_theta_grid

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
     integer, intent (in) :: n
     integer, intent (out) :: j
     integer, dimension (:), intent (out) :: div
     integer :: i, imax

     do i=2,n
        if (mod(n,i)==0) exit
     end do
     imax = n/i
     j=1
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
     use le_grids, only: leok_le_grids, check_le_grids, negrid, nlambda
     use nonlinear_terms, only: nonlin, cfl, check_nonlinear_terms
     use run_parameters, only: check_run_parameters
     use run_parameters, only: beta, tite, margin, code_delt_max
     use run_parameters, only: nstep, wstar_units
     use species, only: check_species, spec, nspec, has_electron_species
     use theta_grid, only: check_theta_grid
     use theta_grid, only: gb_to_cv, nbset, ntgrid_real => ntgrid
     use theta_grid_params, only: nperiod, ntheta, eps, epsl, rmaj, r_geo
     use theta_grid_params, only: pk, qinp, rhoc, shift, shat
     use theta_grid_params, only: akappa, akappri, tri, tripri

     implicit none
     real :: alne, dbetadrho_spec
     real :: kxfac, drhodpsi
     character (20) :: datestamp, timestamp, zone
     character (200) :: line
     logical :: le_ok = .true.
     integer :: ntgrid, j, nmesh
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

     le_ok=leok_le_grids(report_unit)
     
     if (le_ok) then

        ntgrid = ntgrid_real
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
        end if

        write (report_unit, fmt="('ntgrid :    ',i12)") ntgrid
        write (report_unit, fmt="('nlambda:    ',i12)") nlambda
        write (report_unit, fmt="('negrid :    ',i12)") negrid
        write (report_unit, fmt="('ntheta0:    ',i12)") ntheta0
        write (report_unit, fmt="('naky   :    ',i12)") naky
        write (report_unit, fmt="('ny     :    ',i12)") ny
        write (report_unit, fmt="('nx     :    ',i12)") nx
        write (report_unit, fmt="('nspec  :    ',i12)") nspec
        write (report_unit, fmt="('Number of meshpoints:    ',i12)") nmesh

        call nprocs (nmesh)

     end if

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

    write (report_unit, fmt="(/'------------------------------------------------------------'/)")
    
    call check_le_grids(report_unit,le_ok)

    write (report_unit, fmt="(/'------------------------------------------------------------'/)")
    call check_run_parameters(report_unit)

! diagnostic controls:
    call check_gs2_diagnostics(report_unit)

    call close_output_file (report_unit)

  end subroutine report

  subroutine nprocs (nmesh)
    use nonlinear_terms, only : nonlin
    use species, only : nspec
    use kt_grids, only: gridopt_switch, gridopt_single, gridopt_range, gridopt_specified, gridopt_box, gridopt_xbox
    use kt_grids, only: naky, ntheta0
    use le_grids, only: negrid, nlambda
    use gs2_layouts, only: layout
    implicit none
    real :: fac
    integer, intent (in) :: nmesh
    integer :: nefacs, nlfacs, nkyfacs, nkxfacs, nspfacs
    integer, dimension(:,:), allocatable :: facs
    integer :: npe
    real :: time

    write (report_unit, fmt="('Layout = ',a5)") layout 
    if (nonlin) then
       select case (layout)
       case ('lexys')

          write (report_unit, *) 
          write (report_unit, fmt="('Recommended numbers of processors, time on T3E')") 
          allocate (facs(max(nspec,naky,ntheta0)/2+1,5))
          call factors (nspec, nspfacs, facs(:,1))
          call factors (naky, nkyfacs, facs(:,2))
          call factors (ntheta0, nkxfacs, facs(:,3))
          call factors (negrid, nefacs, facs(:,4))
          call factors (nlambda, nlfacs, facs(:,5))
          do i=1,nspfacs
             npe = facs(i,1)
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= 40.*nmesh/1.e7/npe**0.85
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'s'
          end do
          do i=2,nkyfacs
             npe = facs(i,2)*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= 40.*nmesh/1.e7/npe**0.85
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'y'
          end do
          do i=2,nkxfacs
             npe = facs(i,3)*naky*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= 40.*nmesh/1.e7/npe**0.85
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step  (',a,')')") npe, time,'x'
          end do
          do i=2,nefacs
             npe = facs(i,4)*ntheta0*naky*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= 40.*nmesh/1.e7/npe**0.85
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step  (',a,')')") npe, time,'e'
          end do
          do i=2,nlfacs
             npe = facs(i,5)*negrid*ntheta0*naky*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= 40.*nmesh/1.e7/npe**0.85
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step  (',a,')')") npe, time,'l'
          end do
          deallocate (facs)

       case ('lxyes')

          write (report_unit, *) 
          write (report_unit, fmt="('Recommended numbers of processors, time on SP2')") 
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
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'s'
          end do
          do i=2,nefacs
             npe = facs(i,2)*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'e'
          end do
          do i=2,nkyfacs
             npe = facs(i,3)*negrid*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'y'
          end do
          do i=2,nkxfacs
             npe = facs(i,4)*naky*negrid*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'x'
          end do
          do i=2,nlfacs
             npe = facs(i,5)*naky*ntheta0*negrid*nspec 
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'l'
          end do
          deallocate (facs)

       case ('yxels')

          write (report_unit, *) 
          write (report_unit, fmt="('Recommended numbers of processors:')") 

          allocate (facs(max(nspec,negrid,nlambda)/2+1,5))
          call factors (nspec, nspfacs, facs(:,1))
          call factors (nlambda, nlfacs, facs(:,2))
          call factors (negrid, nefacs, facs(:,3))
          call factors (ntheta0, nkxfacs, facs(:,4))
          call factors (naky, nkyfacs, facs(:,5))
          fac = 3.5*(real(nmesh))**1.1/1.e7
          do i=1,nspfacs
             npe = facs(i,1)
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'s'
          end do
          do i=2,nlfacs
             npe = facs(i,2)*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'l'
          end do
          do i=2,nefacs
             npe = facs(i,3)*nlambda*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'e'
          end do
          do i=2,nkxfacs
             npe = facs(i,4)*negrid*nlambda*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'x'
          end do
          do i=2,nkyfacs
             npe = facs(i,5)*nkxfacs*negrid*nlambda*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'y'
          end do
          deallocate (facs)

       case ('yxles')

          write (report_unit, *) 
          write (report_unit, fmt="('Recommended numbers of processors:')") 

          allocate (facs(max(nspec,negrid,nlambda)/2+1,5))
          call factors (nspec, nspfacs, facs(:,1))
          call factors (negrid, nefacs, facs(:,2))
          call factors (nlambda, nlfacs, facs(:,3))
          call factors (ntheta0, nkxfacs, facs(:,4))
          call factors (naky, nkyfacs, facs(:,5))
          fac = 3.5*(real(nmesh))**1.1/1.e7
          do i=1,nspfacs
             npe = facs(i,1)
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'s'
          end do
          do i=2,nefacs
             npe = facs(i,2)*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'e'
          end do
          do i=2,nlfacs
             npe = facs(i,3)*negrid*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'l'
          end do
          do i=2,nkxfacs
             npe = facs(i,4)*nlambda*negrid*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'x'
          end do
          do i=2,nkyfacs
             npe = facs(i,5)*nkxfacs*nlambda*negrid*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'y'
          end do
          deallocate (facs)

       case ('xyles')
!CMR, 11/11/2009: add processor recommendations for xyles layout
!            NB added recommendations that also parallelise in y and x, 
!                          which may be unwise!
          write (report_unit, *) 
          write (report_unit, fmt="('Recommended numbers of processors:')") 

          allocate (facs(max(nspec,negrid,nlambda,naky,ntheta0)/2+1,5))
          call factors (nspec, nspfacs, facs(:,1))
          call factors (negrid, nefacs, facs(:,2))
          call factors (nlambda, nlfacs, facs(:,3))
          call factors (naky, nkyfacs, facs(:,4))
          call factors (ntheta0, nkxfacs, facs(:,5))
          fac = 3.5*(real(nmesh))**1.1/1.e7
          do i=1,nspfacs
             npe = facs(i,1)
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'s'
          end do
          do i=2,nefacs
             npe = facs(i,2)*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'e'
          end do
          do i=2,nlfacs
             npe = facs(i,3)*negrid*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'l'
          end do
          do i=2,nkyfacs
             npe = facs(i,4)*nlambda*negrid*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'y'
          end do
          do i=2,nkxfacs
             npe = facs(i,5)*naky*nlambda*negrid*nspec
             if (npe .gt. npmax) exit
             time=-9.9e9 ; if (nmesh/npe > ncut) time= fac/npe**0.95
             write (report_unit, fmt="('  npe = ',i8,'    time = ',1pe10.2,'  seconds/time step (',a,')')") npe, time,'x'
          end do
          deallocate (facs)

       end select
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

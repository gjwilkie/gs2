!> A module for calculating and writing gs2 outputs. It can write
!! these outputs both to a netcdf file <run_name>.cdf or to ascii text
!! files. It is controlled via the namelist diagnostics_config. 
!! This module is intended to replace the old gs2_diagnostics module with a 
!! simpler and more structured interface. 
module gs2_diagnostics_new
  use diagnostics_config, only: diagnostics_type
  implicit none

  private
  
  public :: init_gs2_diagnostics_new
  public :: finish_gs2_diagnostics_new
  public :: reset_averages_and_counters
  public :: run_diagnostics
  public :: gnostics
  
  !> Options passed to init_gs2_diagnostics_new 
  public :: diagnostics_init_options_type
  
  type diagnostics_init_options_type
     logical :: parallel_io_capable
     logical :: default_double
     logical :: initialized
     logical :: is_trinity_run
  end type diagnostics_init_options_type

  type(diagnostics_type) :: gnostics 

  logical, parameter :: debug=.false.

contains
  !> Read namelist diagnostics_config, initialise submodules,
  !! open output file 'run_name.cdf' and create dimensions.
  subroutine init_gs2_diagnostics_new(init_options)
    use diagnostics_config, only: init_diagnostics_config
    use volume_averages, only: init_volume_averages
    use diagnostics_fluxes, only: init_diagnostics_fluxes
    use diagnostics_omega, only: init_diagnostics_omega
    use diagnostics_velocity_space, only: init_diagnostics_velocity_space
    use diagnostics_heating, only: init_diagnostics_heating
    use diagnostics_ascii, only: init_diagnostics_ascii
    use diagnostics_antenna, only: init_diagnostics_antenna
    use diagnostics_nonlinear_convergence, only: init_nonlinear_convergence
    use nonlinear_terms, only: nonlin
    use gs2_save, only: save_many
    use file_utils, only: run_name, error_unit
    use mp, only: mp_comm, proc0
    use kt_grids, only: naky, aky
    !This reference needs removing to allow compilation on systems without netcdf
    !It should probably be moved to simpledataio as the rest of the diagnostics
    !doesn't know what backend simpledataio is using (i.e. if it is netcdf or not).
    use netcdf, only: NF90_CLOBBER
    use simpledataio, only: SDATIO_DOUBLE, SDATIO_FLOAT, SDATIO_INT
    use simpledataio, only: sdatio_init, simpledataio_functional
    use simpledataio, only: set_parallel, create_file
    use diagnostics_metadata, only: write_metadata
    use diagnostics_metadata, only: read_input_file_and_get_size
    implicit none
    type(diagnostics_init_options_type), intent(in) :: init_options
    
    if(proc0.and.debug) write (*,*) 'initializing new diagnostics'
    call init_diagnostics_config(gnostics)
    call check_parameters
    call check_restart_file_writeable
    
    call init_volume_averages
    
    ! Eventually we will remove this as we want the new diagnostics
    ! module to be built even if netCDF is not available.
    ! Or do we??
    if (.not. simpledataio_functional()) then
       if (proc0) then
          write (*,*) "WARNING: simpledataio is non-functional. &
               & Setting write_any to false in gs2_diagnostics_new"
       end if
       gnostics%write_any = .false.
    end if

    !!!!!!!!!!!!!!!!!!!!!!!
    !! Adjust other modules
    !!!!!!!!!!!!!!!!!!!!!!!
    save_many = gnostics%save_many

    if (.not. gnostics%write_any) return

    ! For the moment, hardwire these so as not to 
    ! conflict with the old module. 
    !gnostics%save_for_restart = .false.
    !gnostics%save_distfn = .false.
    
    ! Set whether this is a Trinity run.. enforces certain 
    ! calculations
    gnostics%is_trinity_run = init_options%is_trinity_run
    
    gnostics%parallel = .false.
    if (gnostics%enable_parallel) then
       if (init_options%parallel_io_capable) then 
          gnostics%parallel = .true.
       else
          if (proc0) write (*,*) "WARNING: you have selected &
               & enable_parallel but this build does not have &
               & parallel capability."
       end if
    end if
    
    !write (*,*) 'parallel', gnostics%parallel
    if (init_options%default_double) then
       gnostics%rtype = SDATIO_DOUBLE
    else
       gnostics%rtype = SDATIO_FLOAT
    end if
    gnostics%itype = SDATIO_INT
    
    gnostics%user_time_old = 0.0
    
    ! fluxfac is used for summing fields, fluxes etc over ky
    ! Mostly this is not needed, since the average_ky routine in 
    ! volume_averages takes care of the factor... you only need it
    ! if you are manually summing something over ky
    allocate(gnostics%fluxfac(naky))
    gnostics%fluxfac = 0.5
    !<DD>This is only correct if running in box mode surely?
    !    I think this should be if(aky(1)==0.0) fluxfac(1)=1.0 but I may be wrong
    if(aky(1)==0.0) gnostics%fluxfac(1) = 1.0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Open Text Files (if required)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (proc0) call set_ascii_file_switches
    if (proc0) call init_diagnostics_ascii(gnostics%ascii_files)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialise submodules
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    call init_diagnostics_fluxes(gnostics)
    call init_diagnostics_omega(gnostics)
    !if (gnostics%write_max_verr) gnostics%write_verr = .true.
    call init_diagnostics_velocity_space(gnostics)
    call init_diagnostics_antenna(gnostics)
    if (gnostics%write_heating) call init_diagnostics_heating(gnostics)
    call read_input_file_and_get_size(gnostics)
    
    call sdatio_init(gnostics%sfile, trim(run_name)//'.out.nc')
    if (gnostics%parallel) then 
       call set_parallel(gnostics%sfile, mp_comm)
    else
       if (.not. gnostics%serial_netcdf4) then 
          gnostics%sfile%mode = NF90_CLOBBER
       end if
    end if
    if (gnostics%parallel.or.proc0) then 
       call create_file(gnostics%sfile)
       call write_metadata(gnostics)
    end if

    !All procs initialise dimension data but if not parallel IO
    !only proc0 has to add them to file.
    call create_dimensions(gnostics%parallel.or.proc0)

    if (nonlin.and.gnostics%use_nonlin_convergence) call init_nonlinear_convergence(gnostics)

    
    if(proc0.and.debug) write (*,*) 'finished initializing new diagnostics'
  end subroutine init_gs2_diagnostics_new
  
  subroutine check_parameters
    use run_parameters, only: fapar
    use file_utils, only: error_unit
    implicit none
    if ((gnostics%print_line .or. gnostics%write_line) .and. .not.(gnostics%write_fields.and.gnostics%write_omega)) then 
       write (error_unit(), *) 'print_line and write_line require both write_fields and write_omega... enabling'
       gnostics%write_fields = .true.
       gnostics%write_omega = .true.
    end if
    if ((gnostics%print_flux_line .or. gnostics%write_flux_line) .and. .not.gnostics%write_fields) then 
       write (error_unit(), *) 'print_flux_line and write_flux_line require both write_fields ... enabling'
       gnostics%write_fields = .true.
    end if
    if (gnostics%write_jext .and. .not. fapar .gt. epsilon(0.0)) then
       write (*,*) "ERROR: it doesn't make sense to switch on write_jext without apar"
       !SHOULD THIS BE MP_ABORT? COULD WE NOT JUST DISABLE THE DIAGNOSTIC?
       stop 1
    end if

  end subroutine check_parameters


  subroutine check_restart_file_writeable
    use gs2_save, only: restart_writable
    use mp, only: proc0, mp_abort
    logical :: writable
    !Verify restart file can be written
    if((gnostics%save_for_restart.or.gnostics%save_distfn).and.(gnostics%file_safety_check))then
       !Can we write file?
       writable=restart_writable()

       !If we can't write the restart file then we should probably quit
       if((.not.writable).and.gnostics%save_for_restart) &
         call mp_abort("Cannot write to test file, maybe restart_dir &
         & doesn't exist --> Aborting.",to_screen=.true.)

       !If it's just a case of save_distfn then we can carry on but print a useful mesasge
       if((.not.writable).and.gnostics%save_distfn)then
          if(proc0)write(6,'("Warning: Cannot write to test restart_file --> Setting save_distfn=F.")')
          gnostics%save_distfn=.false.
       endif
    endif
  end subroutine check_restart_file_writeable

  !> This subroutine determines which ascii output files are enabled
  !! (i.e., opened, flushed at each write, and then closed).
  !! If an ascii file is not enabled here, writing to it will 
  !! cause some indeterminate unpleasant behaviour
  !!
  !! Note that the .out file is always enabled
  subroutine set_ascii_file_switches
    implicit none
    gnostics%ascii_files%write_to_out   = .true.
    !gnostics%ascii_files%write_to_fields = gnostics%write_fields  .and.  gnostics%write_ascii
    gnostics%ascii_files%write_to_heat   = gnostics%write_heating .and.  gnostics%write_ascii
    gnostics%ascii_files%write_to_heat2  = gnostics%write_heating .and.  gnostics%write_ascii
    gnostics%ascii_files%write_to_lpc    = gnostics%write_verr    .and.  gnostics%write_ascii
    gnostics%ascii_files%write_to_vres   = gnostics%write_verr    .and.  gnostics%write_ascii
    gnostics%ascii_files%write_to_vres2  = gnostics%write_verr    .and.  gnostics%write_ascii
    gnostics%ascii_files%write_to_cres   = gnostics%write_cerr    .and.  gnostics%write_ascii
    gnostics%ascii_files%write_to_phase  = gnostics%write_cross_phase .and.  gnostics%write_ascii
    gnostics%ascii_files%write_to_jext   = gnostics%write_jext    .and.  gnostics%write_ascii
    gnostics%ascii_files%write_to_parity = gnostics%write_parity  .and.  gnostics%write_ascii
    !write (*,*) 'gnostics%ascii_files%write_to_heat2', gnostics%ascii_files%write_to_heat2
  end subroutine set_ascii_file_switches

  !> Close the output file and deallocate arrays
  subroutine finish_gs2_diagnostics_new
    use diagnostics_fluxes, only: finish_diagnostics_fluxes
    use diagnostics_omega, only: finish_diagnostics_omega
    use diagnostics_heating, only: finish_diagnostics_heating
    use diagnostics_ascii, only: finish_diagnostics_ascii
    use diagnostics_config, only: finish_diagnostics_config
    use diagnostics_antenna, only: finish_diagnostics_antenna
    use diagnostics_nonlinear_convergence, only: finish_nonlinear_convergence
    use diagnostics_velocity_space, only: finish_diagnostics_velocity_space
    use nonlinear_terms, only: nonlin
    use dist_fn, only: write_fyx, write_f, write_poly, collision_error
    use mp, only: proc0
    use fields_arrays, only: phinew, bparnew
    use simpledataio, only: closefile
    use unit_tests, only: debug_message
    
    implicit none
    integer, parameter :: verb=3
    if (.not. gnostics%write_any) return
    
    call debug_message(verb, 'gs2_diagnostics_new::finish_gs2_diagnostics_new &
      & calling save_restart_dist_fn')
    call save_restart_dist_fn
    
    call debug_message(verb, 'gs2_diagnostics_new::finish_gs2_diagnostics_new &
      & calling run_old_final_routines')
    call run_old_final_routines
    
    deallocate(gnostics%fluxfac)
    call debug_message(verb, 'gs2_diagnostics_new::finish_gs2_diagnostics_new &
      & finishing submodules')
    call finish_diagnostics_fluxes
    call finish_diagnostics_omega
    call finish_diagnostics_antenna(gnostics)
    call finish_diagnostics_velocity_space(gnostics)
    if (nonlin.and.gnostics%use_nonlin_convergence) call finish_nonlinear_convergence(gnostics)
    if (gnostics%write_heating) call finish_diagnostics_heating(gnostics)
    if (gnostics%parallel .or. proc0) then
       if(proc0.and.debug) write(*,*) "Closing new diagnostics"
       call closefile(gnostics%sfile)
       !if (gnostics%write_movie) call closefile(gnostics%sfilemovie)
    end if
    !if (gnostics%write_ascii .and. proc0) call finish_diagnostics_ascii(gnostics%ascii_files)
    if (proc0) call finish_diagnostics_ascii(gnostics%ascii_files)
    
    ! Random stuff that needs to be put in properly or removed
    if (gnostics%write_gyx) call write_fyx (phinew,bparnew,.true.)
    if (gnostics%write_g) call write_f (.true.)
    if (gnostics%write_lpoly) call write_poly (phinew,bparnew,.true.,gnostics%istep)
    !if (gnostics%write_cerr) call collision_error(phinew,bparnew,.true.)
    
    call finish_diagnostics_config(gnostics)
  end subroutine finish_gs2_diagnostics_new

  subroutine save_restart_dist_fn
    use run_parameters, only: fphi, fapar, fbpar
    use collisions, only: vnmult
    use gs2_save, only: gs2_save_for_restart
    use fields_arrays, only: phinew, bparnew
    use gs2_time, only: user_dt
    use dist_fn_arrays, only: gnew
    use dist_fn_arrays, only: g_adjust
    integer :: istatus
    if (gnostics%save_for_restart) then
       call gs2_save_for_restart (gnew, gnostics%user_time, user_dt, vnmult, istatus, &
            fphi, fapar, fbpar, exit_in=.true.)
    end if
    
    !<DD> Added for saving distribution function
    if (gnostics%save_distfn) then
       !Convert h to distribution function
       call g_adjust(gnew,phinew,bparnew,fphi,fbpar)
       
       !Save dfn, fields and velocity grids to file
       call gs2_save_for_restart (gnew, gnostics%user_time, user_dt, vnmult, istatus, &
            fphi, fapar, fbpar, exit_in=.true.,distfn=.true.)
       
       !Convert distribution function back to h
       call g_adjust(gnew,phinew,bparnew,-fphi,-fbpar)
    end if
    !</DD> Added for saving distribution function
  end subroutine save_restart_dist_fn
  
  subroutine run_diagnostics_to_be_updated
    use fields_arrays, only: phinew, bparnew
    use dist_fn, only: write_fyx, write_f, write_poly, collision_error
    implicit none
    integer :: nwrite_large

    nwrite_large = gnostics%nwrite*gnostics%nwrite_mult
    ! Random stuff that needs to be put in properly or removed
    if (gnostics%write_gyx .and. &
      mod(gnostics%istep,nwrite_large) == 0) call write_fyx (phinew,bparnew,.false.)
    if (gnostics%write_g   .and. mod(gnostics%istep,nwrite_large) == 0) call write_f (.false.)
    if (gnostics%write_lpoly) call write_poly (phinew,bparnew,.false.,gnostics%istep)
    !if (gnostics%write_cerr) call collision_error(phinew,bparnew,.false.)
  end subroutine run_diagnostics_to_be_updated
  
  !> Create or write all variables according to the value of istep:
  !! istep=-1 --> Create all netcdf variables
  !! istep=0 --> Write constant arrays/parameters (e.g. aky) and initial values
  !! istep>0 --> Write variables
  subroutine run_diagnostics(istep, exit)
    use gs2_time, only: user_time
    use mp, only: proc0
    use diagnostics_printout, only: print_flux_line, print_line
    use diagnostics_printout, only: write_flux_line, write_line
    use diagnostics_fluxes, only: calculate_fluxes, write_symmetry, write_parity
    use diagnostics_fields, only: write_fields, write_movie
    use diagnostics_fields, only: dump_fields_periodically
    use diagnostics_moments, only: write_moments, write_full_moments_notgc
    use diagnostics_omega, only: calculate_omega, write_omega
    use diagnostics_velocity_space, only: write_velocity_space_checks
    use diagnostics_velocity_space, only: write_collision_error
    use diagnostics_heating, only: calculate_heating, write_heating
    use diagnostics_geometry, only: write_geometry
    use diagnostics_normalisations, only: write_normalisations
    use diagnostics_nonlinear_convergence, only: check_nonlin_convergence
    use diagnostics_turbulence, only: write_cross_phase, write_correlation
    use diagnostics_turbulence, only: write_correlation_extend
    use diagnostics_antenna, only: write_jext, write_lorentzian
    use diagnostics_ascii, only: flush_output_files
    use diagnostics_dimensions, only: dim_string
    use diagnostics_create_and_write, only: create_and_write_variable
    use collisions, only: vary_vnew
    use nonlinear_terms, only: nonlin
    use species, only: spec, has_electron_species
    use simpledataio, only: increment_start, syncfile, closefile, set_parallel, create_file
    use diagnostics_metadata, only: write_input_file
    implicit none
    integer, intent(in) :: istep
    logical, intent(inout) :: exit
    
    if (.not. gnostics%write_any) return
    
    gnostics%istep = istep
    gnostics%exit = exit
    
    ! If parallel, then everybody writes to netcdf,
    ! otherwise, only proc0. Also, creation of variables
    ! happens when istep == -1
    if (gnostics%parallel .or. proc0) then
       gnostics%create = (istep==-1)
       gnostics%wryte = (istep>-1)
    else
       gnostics%create=.false.
       gnostics%wryte=.false.
    end if
    
    !gnostics%wryte = .false.
    
    ! Sets whether field-like arrays are assumed
    ! to be distributed across processes
    ! This line is a temporary placeholder
    ! till distributed fields are up and running
    gnostics%distributed = gnostics%parallel
    !gnostics%distributed = .false.
    
    gnostics%calculate_fluxes = (gnostics%write_fluxes &
         .or.  gnostics%print_flux_line &
         .or.  gnostics%write_flux_line &
         .or.  gnostics%is_trinity_run)
    
    gnostics%user_time = user_time
    if (istep .eq. 0) gnostics%start_time = user_time
    
    ! Write constants/parameters
    if (istep < 1) then
       call write_dimensions
       call write_geometry(gnostics)
       call write_normalisations(gnostics)
       call write_input_file(gnostics)
    end if
    
    if (istep > 0) then
       call calculate_omega(gnostics)
       if (gnostics%write_heating) call calculate_heating (gnostics)
    end if

    ! EGH to be reinstated when old diagnostics removed.
    !if (mod(istep, gnostics%nsave).eq.0.or.exit) then
      !call save_restart_dist_fn
    !end if

    if (istep==-1.or.mod(istep, gnostics%nwrite).eq.0.or.exit) then
       gnostics%vary_vnew_only = .false.
       if (gnostics%write_omega)  call write_omega (gnostics)
       if (gnostics%write_fields) call write_fields(gnostics)
       if (gnostics%dump_fields_periodically) call dump_fields_periodically(gnostics)
       if (gnostics%calculate_fluxes) call calculate_fluxes(gnostics) ! NB  also writes fluxes if on
       if (gnostics%write_symmetry) call write_symmetry(gnostics)
       if (gnostics%write_parity) call write_parity(gnostics)
       if (gnostics%write_verr) call write_velocity_space_checks(gnostics)
       if (gnostics%write_cerr) call write_collision_error(gnostics) ! NB only ascii atm
       if (gnostics%write_moments) call write_moments(gnostics)
       if (gnostics%write_full_moments_notgc) call write_full_moments_notgc(gnostics)
       if (gnostics%make_movie) call write_movie(gnostics)
       if (gnostics%write_heating) call write_heating(gnostics)
       if (nonlin.and.gnostics%use_nonlin_convergence) call check_nonlin_convergence(gnostics)
       if (gnostics%write_cross_phase.and.has_electron_species(spec)) call write_cross_phase(gnostics)
       if (gnostics%write_jext) call write_jext(gnostics)
       if (gnostics%write_correlation) call write_correlation(gnostics)
       if (gnostics%write_correlation_extend) call write_correlation_extend(gnostics)
       if (gnostics%write_lorentzian) call write_lorentzian(gnostics)
       
       if (gnostics%print_line) call print_line(gnostics)
       if (gnostics%write_line) call write_line(gnostics)
       if (proc0) then
          if (gnostics%print_flux_line) call print_flux_line(gnostics)
          if (gnostics%write_flux_line) call write_flux_line(gnostics)
       end if
       call run_diagnostics_to_be_updated

       ! Finally, write time value and update time index
       call create_and_write_variable(gnostics, gnostics%rtype, "t", &
            trim(dim_string(gnostics%dims%time)), &
            "Values of the time coordinate", "a/v_thr", user_time) 
       if (gnostics%wryte) call increment_start(gnostics%sfile, &
         trim(dim_string(gnostics%dims%time)))
       if (gnostics%wryte) call syncfile(gnostics%sfile)
       if (proc0 .and. gnostics%write_ascii) call flush_output_files(gnostics%ascii_files)
       
       ! Update time used for time averages
       gnostics%user_time_old = gnostics%user_time
    else if (mod(istep, gnostics%ncheck).eq.0) then
       ! These lines cause the automated checking of velocity space resolution
       ! and correction by varying collisionality
       gnostics%vary_vnew_only = .true.
       if (gnostics%write_verr .and. vary_vnew) call write_velocity_space_checks(gnostics)
    end if
    
    exit = gnostics%exit
  end subroutine run_diagnostics

  !> Reset cumulative flux and heating averages
  !! that are used, e.g. for Trinity.
  !! Does not at the moment apply to average
  !! growth rates.
  subroutine reset_averages_and_counters
    use nonlinear_terms, only: nonlin
    use diagnostics_nonlinear_convergence, only: &
      dnc_reset => reset_averages_and_counters
    use diagnostics_fluxes, only: &
      fluxes_reset => reset_averages_and_counters
    use diagnostics_heating, only: &
      heating_reset => reset_averages_and_counters
    use gs2_time, only: user_time
    
! HJL > reset values for convergence condition
    if (nonlin.and.gnostics%use_nonlin_convergence) &
      call dnc_reset(gnostics)
    call fluxes_reset(gnostics)
    call heating_reset(gnostics)

    gnostics%start_time = user_time

  end subroutine reset_averages_and_counters
  
  subroutine run_old_final_routines
    use diagnostics_final_routines,only: do_write_eigenfunc
    use diagnostics_final_routines,only: do_write_final_fields
    use diagnostics_final_routines,only: do_write_kpar
    use diagnostics_final_routines,only: do_write_final_epar
    use diagnostics_final_routines,only: do_write_final_db
    use diagnostics_final_routines,only: do_write_final_moments
    use diagnostics_final_routines,only: do_write_final_antot
    use diagnostics_final_routines,only: do_write_gs
    use diagnostics_final_routines,only: init_par_filter
    use diagnostics_final_routines,only: ntg_out
    use nonlinear_terms, only: nonlin
    use antenna, only: dump_ant_amp
    use mp, only: proc0
    use kt_grids, only: ntheta0, naky
    use theta_grid, only: ntgrid
    use unit_tests, only: debug_message
    implicit none
    complex, dimension (ntheta0, naky) :: phi0
    integer, parameter :: verb=3
    
    if(gnostics%write_kpar.or.gnostics%write_gs) call init_par_filter
    
    ! ntg_out is imported from diagnostics_final_routines
    ntg_out = ntgrid
    if (proc0) then
       call debug_message(verb, 'gs2_diagnostics_new::run_old_final_routines &
         & calling do_write_eigenfunc')
       if (gnostics%write_eigenfunc) call do_write_eigenfunc(gnostics,phi0)
       call debug_message(verb, 'gs2_diagnostics_new::run_old_final_routines &
         & calling do_write_final_fields')
       if (gnostics%write_final_fields) call do_write_final_fields(gnostics)
       call debug_message(verb, 'gs2_diagnostics_new::run_old_final_routines &
         & calling do_write_kpar')
       if (gnostics%write_kpar) call do_write_kpar(gnostics)
       if (gnostics%write_final_epar) call do_write_final_epar(gnostics)
       
       ! definition here assumes we are not using wstar_units
       if (gnostics%write_final_db) call do_write_final_db(gnostics)
    end if
    !Note pass in phase factor phi0 which may not be initialised
    !this is ok as phi0 will be set in routine if not already set
    if (gnostics%write_final_moments) call do_write_final_moments(gnostics,phi0)
    
    if (gnostics%write_final_antot) call do_write_final_antot(gnostics)
    
    if (proc0) call dump_ant_amp
    
    if (nonlin.and.gnostics%write_gs) call do_write_gs(gnostics)
    
  end subroutine run_old_final_routines
  
  subroutine create_dimensions(add_to_file)
    use simpledataio_write, only: real_imaginary_dimension_name
    use kt_grids, only: naky, ntheta0, nx, ny, jtwist_out, box
    use theta_grid, only: ntgrid
    use le_grids, only: negrid, nlambda
    use species, only: nspec
    use diagnostics_metadata, only: inputfile_length

    implicit none
    logical, intent(in) :: add_to_file

    ! Please stick to the convention of CAPITALS for spectral
    ! dimensions and lower case for non-spectral
    
    ! The final two arguments of the add_dimension function call currently
    ! have no effect, but may in the future
    ! Their values are given for documentation purposes
    
    !Initialise the dimension instances and add to file
    !/Spectral space
    call gnostics%dims%kx%init("kx", ntheta0, "The kx dimension", "")
    if(add_to_file) call gnostics%dims%kx%add_to_file(gnostics%sfile)
    call gnostics%dims%ky%init("ky", naky, "The ky dimension", "")
    if(add_to_file) call gnostics%dims%ky%add_to_file(gnostics%sfile)
    !/Parallel space
    call gnostics%dims%theta%init("theta", 2*ntgrid+1, "The theta (parallel) dimension", "")
    if(add_to_file) call gnostics%dims%theta%add_to_file(gnostics%sfile)
    if (box) then 
       !What is this used for?
       call gnostics%dims%theta_ext%init("theta_ext", (2*ntgrid+1)*((ntheta0-1)/jtwist_out+1), &
            "The theta (parallel) dimension along the extended domain", "")
       if(add_to_file) call gnostics%dims%theta_ext%add_to_file(gnostics%sfile)
    end if
    !/Real space
    ! Since these are generated by using the transform2
    ! routines they include the aliased gridpoints (i.e.
    ! there is redundancy)
    if (box .and. .not. (nx .eq. 0 .or. ny .eq. 0)) then 
       call gnostics%dims%xx%init('x',nx,"The x dimension (inc aliased gridpoints)", "")
       if(add_to_file) call gnostics%dims%xx%add_to_file(gnostics%sfile)
       call gnostics%dims%yy%init('y',ny,"The y dimension", "")
       if(add_to_file) call gnostics%dims%yy%add_to_file(gnostics%sfile)
    endif
    !/Velocity space
    call gnostics%dims%energy%init("energy", negrid, "The energy dimension", "")
    if(add_to_file) call gnostics%dims%energy%add_to_file(gnostics%sfile)
    call gnostics%dims%lambda%init("lambda", nlambda, "The pitch angle dimension", "")
    if(add_to_file) call gnostics%dims%lambda%add_to_file(gnostics%sfile)
    call gnostics%dims%vpar%init("vpa",negrid*nlambda,"For writing functions of vparallel", "")
    if(add_to_file) call gnostics%dims%vpar%add_to_file(gnostics%sfile)
    call gnostics%dims%species%init("species", nspec, "The species dimension", "")
    if(add_to_file) call gnostics%dims%species%add_to_file(gnostics%sfile)
    !/Time | Note this is an unlimited dimension so the length is ignored
    call gnostics%dims%time%init("t", 0, "The time dimension", "", is_unlimited_in=.true.)
    if(add_to_file) call gnostics%dims%time%add_to_file(gnostics%sfile)
    !/Numeric/generic
    ! This tells simpledataio that we are using "ri" for real/imaginary (default is "r")
    real_imaginary_dimension_name = 'ri'
    call gnostics%dims%ri%init("ri", 2, "Real and imaginary components", "")
    if(add_to_file) call gnostics%dims%ri%add_to_file(gnostics%sfile)
    call gnostics%dims%generic_2%init("2", 2, "Generic 2", "")
    if(add_to_file) call gnostics%dims%generic_2%add_to_file(gnostics%sfile)
    call gnostics%dims%generic_3%init("3", 3, "Generic 3", "")
    if(add_to_file) call gnostics%dims%generic_3%add_to_file(gnostics%sfile)
    call gnostics%dims%generic_4%init("4", 4, "Generic 4", "")
    if(add_to_file) call gnostics%dims%generic_4%add_to_file(gnostics%sfile)
    call gnostics%dims%generic_5%init("5", 5, "Generic 5", "")
    if(add_to_file) call gnostics%dims%generic_5%add_to_file(gnostics%sfile)
    ! Special dimension for the input file
    call gnostics%dims%input_file_dim%init("input_file_dim", inputfile_length, "Length of input file", "")
    if(add_to_file) call gnostics%dims%input_file_dim%add_to_file(gnostics%sfile)
  end subroutine create_dimensions

  !subroutine create_dimensions_movie
  !use gs2_layouts, only: yxf_lo
  !call add_dimension(gnostics%sfilemovie, "x", yxf_lo%nx, "The x dimension", "")
  !call add_dimension(gnostics%sfilemovie, "y", yxf_lo%ny, "The y dimension", "")
  !end subroutine create_dimensions_movie

  subroutine write_dimensions
    use kt_grids, only: aky, akx
    use theta_grid, only: theta
    use le_grids, only: al, energy
    use diagnostics_create_and_write, only: create_and_write_variable
    use diagnostics_dimensions, only: dim_string
    implicit none
    call create_and_write_variable(gnostics, gnostics%rtype, "kx", dim_string(gnostics%dims%kx),  &
         "Values of kx, the wavenumber perpendicular to the flux surface ", "1/rho_r", akx)
    call create_and_write_variable(gnostics, gnostics%rtype, "ky", dim_string(gnostics%dims%ky),  &
         "Values of ky, the wavenumber in the direction of grad alpha ", "1/rho_r", aky)
    call create_and_write_variable(gnostics, gnostics%rtype, "theta", dim_string(gnostics%dims%theta),  &
         "Values of theta, the poloidal angle. ", "rad", theta)
    call create_and_write_variable(gnostics, gnostics%rtype, "energy", dim_string(gnostics%dims%energy),  &
         "Values of the energy grid. ", "T_s", energy)
    call create_and_write_variable(gnostics, gnostics%rtype, "lambda", dim_string(gnostics%dims%lambda),  &
         "Values of lambda = energy/magnetic moment", "1/ (2 B_a)", al)
  end subroutine write_dimensions
end module gs2_diagnostics_new

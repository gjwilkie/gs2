!> This module is intended to replace the old gs2_diagnostics module with a 
!! simpler and more structured interface. 
module gs2_diagnostics_new
  use diagnostics_create_and_write
  use simpledataio
  use diagnostics_config, only: diagnostics_type

  implicit none

  !integer, parameter :: gnostics%rtype = SDATIO_DOUBLE

  public :: init_gs2_diagnostics_new
  
  public :: finish_gs2_diagnostics_new

  public :: run_diagnostics

  public :: gnostics

  !> Options passed to init_gs2_diagnostics_new 
  public :: diagnostics_init_options_type

  type diagnostics_init_options_type
    logical :: parallel_io
    logical :: default_double
    logical :: initialized
  end type diagnostics_init_options_type

  private 

  type(diagnostics_type) :: gnostics


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
    use file_utils, only: run_name, error_unit
    use mp, only: mp_comm, proc0
    use kt_grids, only: naky, aky
    type(diagnostics_init_options_type), intent(in) :: init_options

    if(proc0) write (*,*) 'initializing new diagnostics'
    call init_diagnostics_config(gnostics)
    call check_parameters

    call init_volume_averages

    if (.not. simpledataio_functional()) then
      if (proc0) then
        write (*,*) "WARNING: simpledataio is non-functional. &
         & Setting write_any to false in gs2_diagnostics_new"
     end if
     gnostics%write_any = .false.
    end if
    if (.not. gnostics%write_any) return

    gnostics%parallel = init_options%parallel_io
    if (init_options%default_double) then
      gnostics%rtype = SDATIO_DOUBLE
    else
      gnostics%rtype = SDATIO_FLOAT
    end if
    !write (*,*) 'gnostics%rtype', gnostics%rtype, 'doub', SDATIO_DOUBLE, 'float', SDATIO_FLOAT

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



    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialise submodules
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    call init_diagnostics_fluxes(gnostics)
    call init_diagnostics_omega(gnostics)
    !if (gnostics%write_max_verr) gnostics%write_verr = .true.
    call init_diagnostics_velocity_space(gnostics)
    call init_diagnostics_antenna(gnostics)
    if (gnostics%write_heating) call init_diagnostics_heating(gnostics)

    !gnostics%create = .true.
   ! Integer below gives the sdatio type 
   ! which corresponds to a gs2 real
    !gnostics%rtype = SDATIO_DOUBLE

    if (gnostics%parallel) then
      call createfile_parallel(gnostics%sfile, trim(run_name)//'.cdf', mp_comm)
      !if (gnostics%write_movie) &
        !call createfile_parallel(gnostics%sfilemovie, trim(run_name)//'.movie.cdf', mp_comm)

    else if (proc0) then
      call createfile(gnostics%sfile, trim(run_name)//'.cdf')
      !if (gnostics%write_movie) &
        !call createfile(gnostics%sfilemovie, trim(run_name)//'.movie.cdf')
    end if

    if (gnostics%parallel .or. proc0) then 
      call create_dimensions
      !if (gnostics%write_movie) call create_dimensions_movie
    end if
    !if (gnostics%write_ascii) then 
      if (proc0) call set_ascii_file_switches
      if (proc0) call init_diagnostics_ascii(gnostics%ascii_files)
    !end if

    if (nonlin.and.gnostics%use_nonlin_convergence) call init_nonlinear_convergence(gnostics)



  end subroutine init_gs2_diagnostics_new
  
  subroutine check_parameters
    use run_parameters, only: fapar
    use file_utils, only: error_unit
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
      stop 1
    end if
  end subroutine check_parameters

  !> This subroutine determines which ascii output files are enabled
  !! (i.e., opened, flushed at each write, and then closed).
  !! If an ascii file is not enabled here, writing to it will 
  !! cause some indeterminate unpleasant behaviour
  !!
  !! Note that the .out file is always enabled
  subroutine set_ascii_file_switches
    gnostics%ascii_files%write_to_out   = .true.
    gnostics%ascii_files%write_to_fields = gnostics%write_fields  .and.  gnostics%write_ascii
    gnostics%ascii_files%write_to_heat   = gnostics%write_heating .and.  gnostics%write_ascii
    gnostics%ascii_files%write_to_heat2  = gnostics%write_heating .and.  gnostics%write_ascii
    gnostics%ascii_files%write_to_lpc    = gnostics%write_verr    .and.  gnostics%write_ascii
    gnostics%ascii_files%write_to_vres   = gnostics%write_verr    .and.  gnostics%write_ascii
    gnostics%ascii_files%write_to_vres2  = gnostics%write_verr    .and.  gnostics%write_ascii
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
    use nonlinear_terms, only: nonlin
    use dist_fn, only: write_fyx, write_f, write_poly, collision_error
    use dist_fn_arrays, only: g_adjust
    use mp, only: proc0
    use fields_arrays, only: phinew, bparnew
    use gs2_time, only: user_dt
    use dist_fn_arrays, only: gnew
    use run_parameters, only: fphi, fapar, fbpar
    use collisions, only: vnmult
    use gs2_save, only: gs2_save_for_restart
    integer :: istatus
    if (.not. gnostics%write_any) return

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

    call run_old_final_routines

    deallocate(gnostics%fluxfac)
    call finish_diagnostics_fluxes
    call finish_diagnostics_omega
    call finish_diagnostics_antenna(gnostics)
    if (nonlin.and.gnostics%use_nonlin_convergence) call finish_nonlinear_convergence(gnostics)
    if (gnostics%write_heating) call finish_diagnostics_heating(gnostics)
    if (gnostics%parallel .or. proc0) then
      if(proc0)write (*,*) "Closing new diagnostics"
      call closefile(gnostics%sfile)
      !if (gnostics%write_movie) call closefile(gnostics%sfilemovie)
    end if
    !if (gnostics%write_ascii .and. proc0) call finish_diagnostics_ascii(gnostics%ascii_files)
    if (proc0) call finish_diagnostics_ascii(gnostics%ascii_files)

    ! Random stuff that needs to be put in properly or removed
    if (gnostics%write_gyx) call write_fyx (phinew,bparnew,.true.)
    if (gnostics%write_g) call write_f (.true.)
    if (gnostics%write_lpoly) call write_poly (phinew,bparnew,.true.,gnostics%istep)
    if (gnostics%write_cerr) call collision_error(phinew,bparnew,.true.)

    call finish_diagnostics_config(gnostics)
  end subroutine finish_gs2_diagnostics_new

  subroutine run_diagnostics_to_be_updated
    use fields_arrays, only: phinew, bparnew
    use dist_fn, only: write_fyx, write_f, write_poly, collision_error
    ! Random stuff that needs to be put in properly or removed
    if (gnostics%write_gyx .and. mod(gnostics%istep,gnostics%nwrite_large) == 0) call write_fyx (phinew,bparnew,.false.)
    if (gnostics%write_g   .and. mod(gnostics%istep,gnostics%nwrite_large) == 0) call write_f (.false.)
    if (gnostics%write_lpoly) call write_poly (phinew,bparnew,.false.,gnostics%istep)
    if (gnostics%write_cerr) call collision_error(phinew,bparnew,.false.)
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
    use diagnostics_heating, only: calculate_heating, write_heating
    use diagnostics_geometry, only: write_geometry
    use diagnostics_nonlinear_convergence, only: check_nonlin_convergence
    use diagnostics_turbulence, only: write_cross_phase, write_correlation
    use diagnostics_turbulence, only: write_correlation_extend
    use diagnostics_antenna, only: write_jext, write_lorentzian
    use diagnostics_ascii, only: flush_output_files
    use collisions, only: vary_vnew
    use nonlinear_terms, only: nonlin
    use species, only: spec, has_electron_species, nspec
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

    ! Need to also add 'if Trinity run' condition to this
    gnostics%calculate_fluxes = (gnostics%write_fluxes &
                           .or.  gnostics%print_flux_line &
                           .or.  gnostics%write_flux_line)

    gnostics%user_time = user_time
    
    ! Write constants/parameters
    if (istep < 1) then
      call write_dimensions
      call write_geometry(gnostics)
    end if
    
    if (istep > 0) then
      call calculate_omega(gnostics)
      if (gnostics%write_heating) call calculate_heating (gnostics)
    end if




    if (istep==-1.or.mod(istep, gnostics%nwrite).eq.0.or.exit) then
      gnostics%vary_vnew_only = .false.
      if (gnostics%write_omega)  call write_omega (gnostics)
      if (gnostics%write_fields) call write_fields(gnostics)
      if (gnostics%dump_fields_periodically) call dump_fields_periodically(gnostics)
      if (gnostics%calculate_fluxes) call calculate_fluxes(gnostics)
      if (gnostics%write_symmetry) call write_symmetry(gnostics)
      if (gnostics%write_parity) call write_parity(gnostics)
      if (gnostics%write_verr) call write_velocity_space_checks(gnostics)
      if (gnostics%write_moments) call write_moments(gnostics)
      if (gnostics%write_full_moments_notgc) call write_full_moments_notgc(gnostics)
      if (gnostics%write_movie) call write_movie(gnostics)
      if (gnostics%write_heating) call write_heating(gnostics)
      if (nonlin.and.gnostics%use_nonlin_convergence) call check_nonlin_convergence(gnostics)
      if (gnostics%write_cross_phase.and.has_electron_species(spec)) call write_cross_phase(gnostics)
      if (gnostics%write_jext) call write_jext(gnostics)
      if (gnostics%write_correlation) call write_correlation(gnostics)
      if (gnostics%write_correlation_extend) call write_correlation_extend(gnostics)
      if (gnostics%write_lorentzian) call write_lorentzian(gnostics)
!
      if (proc0) then
        if (gnostics%print_line) call print_line(gnostics)
        if (gnostics%write_line) call write_line(gnostics)
        if (gnostics%print_flux_line) call print_flux_line(gnostics)
        if (gnostics%write_flux_line) call write_flux_line(gnostics)
      end if
      call run_diagnostics_to_be_updated

      ! Finally, write time value and update time index
      call create_and_write_variable(gnostics, gnostics%rtype, "t", "t", &
         "Values of the time coordinate", "a/v_thr", user_time) 
      if (gnostics%wryte) call increment_start(gnostics%sfile, "t")
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

  subroutine run_old_final_routines
    use diagnostics_final_routines,only: do_write_eigenfunc
    use diagnostics_final_routines,only: do_write_final_fields
    use diagnostics_final_routines,only: do_write_kpar
    use diagnostics_final_routines,only: do_write_final_epar
    use diagnostics_final_routines,only: do_write_final_db
    use diagnostics_final_routines,only: do_write_final_moments
    use diagnostics_final_routines,only: do_write_final_antot
    use diagnostics_final_routines,only: do_write_gs
    use diagnostics_final_routines,only: do_write_geom
    use diagnostics_final_routines,only: init_par_filter
    use diagnostics_final_routines,only: ntg_out
    use nonlinear_terms, only: nonlin
    use antenna, only: dump_ant_amp
    use mp, only: proc0
    use kt_grids, only: ntheta0, naky
    use theta_grid, only: ntgrid
    complex, dimension (ntheta0, naky) :: phi0

    if(gnostics%write_kpar.or.gnostics%write_gs) call init_par_filter

    ! ntg_out is imported from diagnostics_final_routines
    ntg_out = ntgrid
    if (proc0) then
       if (gnostics%write_eigenfunc) call do_write_eigenfunc(gnostics,phi0)
       if (gnostics%write_final_fields) call do_write_final_fields(gnostics)
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
    
    if (proc0) call do_write_geom

  end subroutine run_old_final_routines


  subroutine create_dimensions
    use kt_grids, only: naky, ntheta0, nx, ny, jtwist_out, box
    use theta_grid, only: ntgrid
    use le_grids, only: negrid, nlambda
    use species, only: nspec

    ! Please stick to the convention of CAPITALS for spectral
    ! dimensions and lower case for non-spectral

    ! The final two arguments of the add_dimension function call currently
    ! have no effect, but may in the future
    ! Their values are given for documentation purposes

    call add_dimension(gnostics%sfile, "X", ntheta0, "The kx dimension", "")
    call add_dimension(gnostics%sfile, "Y", naky, "The ky dimension", "")

    call add_dimension(gnostics%sfile, "z", 2*ntgrid+1, "The theta (parallel) dimension", "")
    call add_dimension(gnostics%sfile, "e", negrid, "", "")
    call add_dimension(gnostics%sfile, "l", nlambda, "", "")
    call add_dimension(gnostics%sfile, "s", nspec, "", "")
    call add_dimension(gnostics%sfile, "r", 2, "Real and imaginary parts", "")
    call add_dimension(gnostics%sfile, "t", SDATIO_UNLIMITED, "", "")


    ! Some specialised dimensions for specific cases
    call add_dimension(gnostics%sfile, "v", negrid*nlambda, "For writing functions of vparallel", "")
    if (box) then 
      call add_dimension(gnostics%sfile, "j", (2*ntgrid+1)*((ntheta0-1)/jtwist_out+1), &
        "The theta (parallel) dimension along the extended domain", "")
    end if

    ! A set of generic dimensions for writing arrays of data 
    ! which are grouped together for convenience like the 
    ! velocity space diagnostics
    call add_dimension(gnostics%sfile, "2", 2, "", "")
    call add_dimension(gnostics%sfile, "3", 3, "", "")
    call add_dimension(gnostics%sfile, "4", 4, "", "")
    call add_dimension(gnostics%sfile, "5", 5, "", "")

    ! Dimensions for writing quantities in real space
    ! Since these are generated by using the transform2
    ! routines they include the aliased gridpoints (i.e.
    ! there is redundancy)
    call add_dimension(gnostics%sfile, "x", nx, "The x dimension (inc aliased gridpoints)", "")
    call add_dimension(gnostics%sfile, "y", ny, "The y dimension", "")

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

    call create_and_write_variable(gnostics, gnostics%rtype, "kx",     "X",  &
        "Values of kx, the wavenumber perpendicular to the flux surface ", "1/rho_r", akx)
    call create_and_write_variable(gnostics, gnostics%rtype, "ky",     "Y",  &
        "Values of ky, the wavenumber in the direction of grad alpha ", "1/rho_r", aky)
    call create_and_write_variable(gnostics, gnostics%rtype, "theta",  "z",  &
        "Values of theta, the poloidal angle. ", "rad", theta)
    call create_and_write_variable(gnostics, gnostics%rtype, "energy", "e",  &
        "Values of the energy grid. ", "T_s", energy)
    call create_and_write_variable(gnostics, gnostics%rtype, "lambda", "l",  &
        "Values of lambda = energy/magnetic moment", "1/ (2 B_a)", al)
  end subroutine write_dimensions




end module gs2_diagnostics_new

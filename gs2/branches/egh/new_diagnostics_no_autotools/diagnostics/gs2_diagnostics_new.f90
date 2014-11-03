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
    use diagnostics_write_fluxes, only: init_diagnostics_write_fluxes
    use diagnostics_write_omega, only: init_diagnostics_write_omega
    use diagnostics_write_velocity_space_checks, only: init_diagnostics_write_velocity_space_checks
    use diagnostics_heating, only: init_diagnostics_heating
    use diagnostics_ascii, only: init_diagnostics_ascii
    use file_utils, only: run_name
    use mp, only: mp_comm, proc0
    type(diagnostics_init_options_type), intent(in) :: init_options
    if(proc0) write (*,*) 'initializing new diagnostics'
    call init_diagnostics_config(gnostics)
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




    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialise submodules
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    call init_diagnostics_write_fluxes
    call init_diagnostics_write_omega(gnostics)
    !if (gnostics%write_max_verr) gnostics%write_verr = .true.
    call init_diagnostics_write_velocity_space_checks(gnostics)
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
    if (gnostics%write_ascii) then 
      call set_ascii_file_switches
      call init_diagnostics_ascii(gnostics%ascii_files)
    end if
  end subroutine init_gs2_diagnostics_new
  
  subroutine set_ascii_file_switches
    if (gnostics%write_fields) gnostics%ascii_files%write_to_fields=.true.
    if (gnostics%write_heating) gnostics%ascii_files%write_to_heat=.true.
    if (gnostics%write_heating) gnostics%ascii_files%write_to_heat2=.true.
  end subroutine set_ascii_file_switches


  !> Close the output file and deallocate arrays
  subroutine finish_gs2_diagnostics_new
    use diagnostics_write_fluxes, only: finish_diagnostics_write_fluxes
    use diagnostics_write_omega, only: finish_diagnostics_write_omega
    use diagnostics_heating, only: finish_diagnostics_heating
    use mp, only: proc0
    if (.not. gnostics%write_any) return
    call finish_diagnostics_write_fluxes
    call finish_diagnostics_write_omega
    if (gnostics%write_heating) call finish_diagnostics_heating(gnostics)
    if (gnostics%parallel .or. proc0) then
      if(proc0)write (*,*) "Closing new diagnostics"
      call closefile(gnostics%sfile)
      !if (gnostics%write_movie) call closefile(gnostics%sfilemovie)
    end if
  end subroutine finish_gs2_diagnostics_new

  !> Create or write all variables according to the value of istep:
  !! istep=-1 --> Create all netcdf variables
  !! istep=0 --> Write constant arrays/parameters (e.g. aky) and initial values
  !! istep>0 --> Write variables
  subroutine run_diagnostics(istep, exit)
    use gs2_time, only: user_time
    use mp, only: proc0
    use diagnostics_write_fluxes, only: write_fluxes
    use diagnostics_write_fields, only: write_fields, write_movie
    use diagnostics_write_moments, only: write_moments
    use diagnostics_write_omega, only: calculate_omega, write_omega
    use diagnostics_write_velocity_space_checks, only: write_velocity_space_checks
    use diagnostics_heating, only: calculate_heating, write_heating
    use diagnostics_geometry, only: write_geometry
    integer, intent(in) :: istep
    logical, intent(inout) :: exit

    if (.not. gnostics%write_any) return

    gnostics%istep = istep
    gnostics%exit = exit

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
      if (gnostics%write_omega)  call write_omega (gnostics)
      if (gnostics%write_fields) call write_fields(gnostics)
      if (gnostics%write_fluxes) call write_fluxes(gnostics)
      if (gnostics%write_verr) call write_velocity_space_checks(gnostics)
      if (gnostics%write_moments) call write_moments(gnostics)
      if (gnostics%write_movie) call write_movie(gnostics)
      if (gnostics%write_heating) call write_heating(gnostics)
!
      ! Finally, write time value and update time index
      call create_and_write_variable(gnostics, gnostics%rtype, "t", "t", &
         "Values of the time coordinate", "a/v_thr", user_time) 
      if (gnostics%wryte) call increment_start(gnostics%sfile, "t")
      if (gnostics%wryte) call syncfile(gnostics%sfile)
    end if

    exit = gnostics%exit

  end subroutine run_diagnostics

  subroutine create_dimensions
    use kt_grids, only: naky, ntheta0, nx, ny
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

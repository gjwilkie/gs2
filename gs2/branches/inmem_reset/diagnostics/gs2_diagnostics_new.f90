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




    call init_diagnostics_write_fluxes
    call init_diagnostics_write_omega(gnostics)

    !gnostics%create = .true.
   ! Integer below gives the sdatio type 
   ! which corresponds to a gs2 real
    !gnostics%rtype = SDATIO_DOUBLE

    if (gnostics%parallel) then
      call createfile_parallel(gnostics%sfile, trim(run_name)//'.cdf', mp_comm)
    else if (proc0) then
      call createfile(gnostics%sfile, trim(run_name)//'.cdf')
    end if

    if (gnostics%parallel .or. proc0) call create_dimensions
  end subroutine init_gs2_diagnostics_new

  !> Close the output file and deallocate arrays
  subroutine finish_gs2_diagnostics_new
    use diagnostics_write_fluxes, only: finish_diagnostics_write_fluxes
    use diagnostics_write_omega, only: finish_diagnostics_write_omega
    use mp, only: proc0
    if (.not. gnostics%write_any) return
    call finish_diagnostics_write_fluxes
    call finish_diagnostics_write_omega
    if (gnostics%parallel .or. proc0) then
      if(proc0)write (*,*) "Closing new diagnostics"
      call closefile(gnostics%sfile)
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
    use diagnostics_write_fields, only: write_fields
    use diagnostics_write_omega, only: calculate_omega, write_omega
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
    gnostics%distributed = gnostics%parallel

    ! Write constants/parameters
    if (istep < 1) then
      call write_dimensions
      call write_geometry
    end if
    
    if (istep > 0) then
      call calculate_omega(gnostics)
    end if



    if (istep==-1.or.mod(istep, gnostics%nwrite).eq.0.or.exit) then
      if (gnostics%write_omega)  call write_omega (gnostics)
      if (gnostics%write_fields) call write_fields(gnostics)
      if (gnostics%write_fluxes) call write_fluxes(gnostics)
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
    use kt_grids, only: naky, ntheta0
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
    call add_dimension(gnostics%sfile, "r", 2, "", "")
    call add_dimension(gnostics%sfile, "t", SDATIO_UNLIMITED, "", "")

  end subroutine create_dimensions

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

  subroutine write_geometry
    use theta_grid, only: bmag, gradpar, gbdrift, gbdrift0, &
         cvdrift, cvdrift0, gds2, gds21, gds22, grho, jacob, &
         shat, drhodpsi, eps, cdrift, cdrift0, qval
    call create_and_write_variable(gnostics, gnostics%rtype, "bmag", "z", &
      "Values of bmag, the magnitude of the magnetic field ", "B_a", bmag)
    call create_and_write_variable(gnostics, gnostics%rtype, "gradpar", "z", &
      "Values of gradpar, which multiplies the parallel derivative", "a", gradpar)
    call create_and_write_variable(gnostics, gnostics%rtype, "gbdrift", "z", &
      "Values of gbdrift, the magnetic gradient drift ", "TBC", gbdrift)
    call create_and_write_variable(gnostics, gnostics%rtype, "gbdrift0", "z", &
      "Values of gbdrift0, ", "TBC", gbdrift0)
    call create_and_write_variable(gnostics, gnostics%rtype, "cvdrift", "z", &
      "Values of cvdrift, the magnetic curvature drift", "TBC", cvdrift)
    call create_and_write_variable(gnostics, gnostics%rtype, "cvdrift0", "z", &
      "Values of cvdrift0, ", "TBC", cvdrift0)
    call create_and_write_variable(gnostics, gnostics%rtype, "gds2", "z", &
      "Values of gds2, ", "TBC", gds2)
    call create_and_write_variable(gnostics, gnostics%rtype, "gds21", "z", &
      "Values of gds21, ", "TBC", gds21)
    call create_and_write_variable(gnostics, gnostics%rtype, "gds22", "z", &
      "Values of gds22, ", "TBC", gds22)
    call create_and_write_variable(gnostics, gnostics%rtype, "grho", "z", &
      "Values of grho, ", "TBC", grho)
    call create_and_write_variable(gnostics, gnostics%rtype, "jacob", "z", &
      "Values of jacob, ", "TBC", jacob)
    call create_and_write_variable(gnostics, gnostics%rtype, "shat", "", &
      "Values of shat, the magnetic shear", "TBC", shat)
    call create_and_write_variable(gnostics, gnostics%rtype, "drhodpsi", "", &
      "Values of drhodpsi, ", "TBC", drhodpsi)
    call create_and_write_variable(gnostics, gnostics%rtype, "eps", "", &
      "Values of eps, ", "TBC", eps)
    call create_and_write_variable(gnostics, gnostics%rtype, "cdrift", "z", &
      "Values of cdrift, ", "TBC", cdrift)
    call create_and_write_variable(gnostics, gnostics%rtype, "cdrift0", "z", &
      "Values of cdrift0, ", "TBC", cdrift0)
    call create_and_write_variable(gnostics, gnostics%rtype, "qval", "", &
      "Values of qval, ", "TBC", qval)
  end subroutine write_geometry



end module gs2_diagnostics_new

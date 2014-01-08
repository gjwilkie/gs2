!> This module is intended to replace the old gs2_diagnostics module with a 
!! simpler and more structured interface. 
module gs2_diagnostics_new
  use diagnostics_create_and_write
  use simpledataio
  use diagnostics_config, only: diagnostics_type

  implicit none

  type(diagnostics_type) :: gnostics
  !integer, parameter :: gnostics%rtype = SDATIO_DOUBLE

  !> Read namelist diagnostics_config, initialise submodules,
  !! open output file 'run_name.cdf' and create dimensions.
  public :: init_gs2_diagnostics_new
  
  !> Close output file and finish submodules
  public :: finish_gs2_diagnostics_new

  !> Create or write all variables according to the value of istep:
  !!   istep=-1 --> Create all netcdf variables
  !!   istep=0 --> Write constant arrays/parameters (e.g. aky) and initial values
  !!   istep>0 --> Write variables
  public :: run_diagnostics

  private 


contains
  subroutine init_gs2_diagnostics_new(parallel_io)
    use diagnostics_config, only: init_diagnostics_config
    use volume_averages, only: init_volume_averages
    use diagnostics_write_fluxes, only: init_diagnostics_write_fluxes
    use file_utils, only: run_name
    use mp, only: mp_comm, proc0
    logical, intent(in) :: parallel_io

    call init_diagnostics_config(gnostics)
    call init_volume_averages

    gnostics%parallel = parallel_io


    if (.not. gnostics%write_any) return


    call init_diagnostics_write_fluxes

    !gnostics%create = .true.
   !> Integer below gives the sdatio type 
   !! which corresponds to a gs2 real
    gnostics%rtype = SDATIO_DOUBLE

    if (gnostics%parallel) then
      call createfile_parallel(gnostics%sfile, trim(run_name)//'.cdf', mp_comm)
    else if (proc0) then
      call createfile(gnostics%sfile, trim(run_name)//'.cdf')
    end if

    if (gnostics%parallel .or. proc0) call create_dimensions
  end subroutine init_gs2_diagnostics_new

  subroutine finish_gs2_diagnostics_new
    use diagnostics_write_fluxes, only: finish_diagnostics_write_fluxes
    use mp, only: proc0
    if (.not. gnostics%write_any) return
    call finish_diagnostics_write_fluxes
    if (gnostics%parallel .or. proc0) call closefile(gnostics%sfile)
  end subroutine finish_gs2_diagnostics_new

  subroutine run_diagnostics(istep, exit)
    use gs2_diagnostics, only: nwrite
    use gs2_time, only: user_time
    use mp, only: proc0
    use diagnostics_write_fluxes, only: write_fluxes
    use diagnostics_write_fields, only: write_fields
    integer, intent(in) :: istep
    logical, intent(inout) :: exit

    if (.not. gnostics%write_any) return

    gnostics%istep = istep

    if (gnostics%parallel .or. proc0) then
      gnostics%create = (istep==-1)
      gnostics%wryte = (istep>-1)
    else
      gnostics%create=.false.
      gnostics%wryte=.false.
    end if

    ! Sets whether field-like arrays are assumed
    ! to be distributed across processes
    gnostics%distributed = gnostics%parallel

    ! Write constants/parameters
    if (istep < 1) then
      call write_dimensions
      call write_geometry
    end if


    if (istep==-1.or.mod(istep, gnostics%nwrite).eq.0.or.exit) then
      if (gnostics%write_fields) call write_fields(gnostics)
      if (gnostics%write_fluxes) call write_fluxes(gnostics)
!
      ! Finally, write time value and update time index
      call create_and_write_variable(gnostics, SDATIO_DOUBLE, "t", "t", &
         "Values of the time coordinate", "a/v_thr", user_time) 
      if (gnostics%wryte) call increment_start(gnostics%sfile, "t")
    end if
  end subroutine run_diagnostics

  subroutine create_dimensions
    use kt_grids, only: naky, ntheta0
    use theta_grid, only: ntgrid
    use le_grids, only: negrid, nlambda
    use species, only: nspec

    call add_dimension(gnostics%sfile, "x", ntheta0, "", "")
    call add_dimension(gnostics%sfile, "y", naky, "", "")
    call add_dimension(gnostics%sfile, "z", 2*ntgrid+1, "", "")
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

    call create_and_write_variable(gnostics, gnostics%rtype, "kx",     "x",  &
        "Values of kx, the wavenumber perpendicular to the flux surface ", "1/rho_r", akx)
    call create_and_write_variable(gnostics, gnostics%rtype, "ky",     "y",  &
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
    call create_and_write_variable(gnostics, gnostics%rtype, "shat", "z", &
      "Values of shat, the magnetic shear", "TBC", shat)
    call create_and_write_variable(gnostics, gnostics%rtype, "drhodpsi", "z", &
      "Values of drhodpsi, ", "TBC", drhodpsi)
    call create_and_write_variable(gnostics, gnostics%rtype, "eps", "z", &
      "Values of eps, ", "TBC", eps)
    call create_and_write_variable(gnostics, gnostics%rtype, "cdrift", "z", &
      "Values of cdrift, ", "TBC", cdrift)
    call create_and_write_variable(gnostics, gnostics%rtype, "cdrift0", "z", &
      "Values of cdrift0, ", "TBC", cdrift0)
    call create_and_write_variable(gnostics, gnostics%rtype, "qval", "z", &
      "Values of qval, ", "TBC", qval)
  end subroutine write_geometry



end module gs2_diagnostics_new

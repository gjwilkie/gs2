!> This module is intended to replace the old gs2_diagnostics module with a 
!! simpler and more structured interface. All of the technicalities of dealing
!! with the netcdf file are in diagnostics_write and
!! diagnostics_create_and_write, and all the calculations 
!! required are in diagnostics_calculate. The naming convention of the
!! subroutines in diagnostics_write follows the naming conventions of variables
!! within the netcdf file itself.
module gs2_diagnostics_new
  use diagnostics_create_and_write
  use simpledataio

  implicit none

  type(sdatio_file) :: netcdf_file
  integer, parameter :: REAL_TYPE = SDATIO_DOUBLE



contains
  subroutine init_gs2_diagnostics_new
    use volume_averages, only: init_volume_averages
    use diagnostics_write_fluxes, only: init_diagnostics_write_fluxes
    use file_utils, only: run_name
    use mp, only: mp_comm
    call createfile_parallel(netcdf_file, trim(run_name)//'.cdf', mp_comm)
    call create_dimensions
    call init_volume_averages
    call init_diagnostics_write_fluxes
  end subroutine init_gs2_diagnostics_new

  subroutine finish_gs2_diagnostics_new
    call closefile(netcdf_file)
  end subroutine finish_gs2_diagnostics_new

  subroutine run_diagnostics(istep, exit)
    use gs2_diagnostics, only: nwrite
    use gs2_time, only: user_time
    use mp, only: proc0
    use diagnostics_write_fluxes, only: write_fluxes
    integer, intent(in) :: istep
    logical, intent(inout) :: exit

    if (mod(istep, nwrite).eq.0.or.exit) then
      call write_fields(istep)
      call write_fluxes(netcdf_file, istep)

      ! Finally, write time value and update time index
      call create_and_write_variable(netcdf_file, SDATIO_DOUBLE, "t", "t", &
         "Values of the time coordinate", "a/v_thr", user_time) 
      call increment_start(netcdf_file, "t")
    end if
  end subroutine run_diagnostics

  subroutine create_dimensions
    use kt_grids, only: naky, ntheta0, aky, akx
    use theta_grid, only: ntgrid, theta
    use le_grids, only: negrid, nlambda, al, energy
    use species, only: nspec

    call add_dimension(netcdf_file, "x", ntheta0, "", "")
    call add_dimension(netcdf_file, "y", naky, "", "")
    call add_dimension(netcdf_file, "z", 2*ntgrid+1, "", "")
    call add_dimension(netcdf_file, "e", negrid, "", "")
    call add_dimension(netcdf_file, "l", nlambda, "", "")
    call add_dimension(netcdf_file, "s", nspec, "", "")
    call add_dimension(netcdf_file, "r", 2, "", "")
    call add_dimension(netcdf_file, "t", SDATIO_UNLIMITED, "", "")

    call create_and_write_variable(netcdf_file, REAL_TYPE, "kx",     "x",  &
        "Values of kx, the wavenumber perpendicular to the flux surface ", "1/rho_r", akx)
    call create_and_write_variable(netcdf_file, REAL_TYPE, "ky",     "y",  &
        "Values of ky, the wavenumber in the direction of grad alpha ", "1/rho_r", aky)
    call create_and_write_variable(netcdf_file, REAL_TYPE, "theta",  "z",  &
        "Values of theta, the poloidal angle. ", "rad", theta)
    call create_and_write_variable(netcdf_file, REAL_TYPE, "energy", "e",  &
        "Values of the energy grid. ", "T_s", energy)
    call create_and_write_variable(netcdf_file, REAL_TYPE, "lambda", "l",  &
        "Values of lambda = energy/magnetic moment", "1/ (2 B_a)", al)
  end subroutine create_dimensions

  subroutine write_fields(istep)
    use diagnostics_write_fields, only: fields_local
    use diagnostics_write_fields, only: write_standard_field_properties
    use fields_arrays, only: phinew, aparnew, bparnew
    use run_parameters, only: fphi, fapar, fbpar
    integer, intent(in) :: istep

    fields_local = .false.
    
    if (fphi >epsilon(0.0)) call write_standard_field_properties(netcdf_file, &
      'phi',  'The electrostatic potential', 'T_r/e', phinew)
    if (fapar>epsilon(0.0)) call write_standard_field_properties(netcdf_file, &
      'apar', 'The parallel magnetic potential', '...', aparnew)
    if (fbpar>epsilon(0.0)) call write_standard_field_properties(netcdf_file, &
      'bpar', 'The parallel magnetic potential', '...', bparnew)
  end subroutine write_fields


end module gs2_diagnostics_new

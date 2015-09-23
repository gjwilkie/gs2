!> A module for writing out geometric coefficients and information
!! such as bpar, drifts and plotting coefficients
module diagnostics_geometry
  implicit none
  private
  public :: write_geometry
contains
  subroutine write_geometry(gnostics)
    use theta_grid, only: bmag, gradpar, gbdrift, gbdrift0, &
         cvdrift, cvdrift0, gds2, gds21, gds22, grho, jacob, &
         shat, drhodpsi, eps, cdrift, cdrift0, qval, shape, theta, ntgrid, Bpol
    use theta_grid, only: Rplot, Zplot, aplot, Rprime, Zprime, aprime, drhodpsi
    use diagnostics_create_and_write, only: create_and_write_variable
    use diagnostics_dimensions, only: dim_string
    use diagnostics_config, only: diagnostics_type
    use file_utils, only: open_output_file, close_output_file
    use mp, only: proc0
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    integer :: i, unit

    call create_and_write_variable(gnostics, gnostics%rtype, "bmag", &
         dim_string([gnostics%dims%theta]), &
         "Values of bmag, the magnitude of the magnetic field", "B_a", bmag)
    call create_and_write_variable(gnostics, gnostics%rtype, "bpol", &
         dim_string([gnostics%dims%theta]), &
         "Values of Bpol, the poloidal magnetic field", "B_a", Bpol)
    call create_and_write_variable(gnostics, gnostics%rtype, "gradpar", &
         dim_string([gnostics%dims%theta]), &
         "Values of gradpar, which multiplies the parallel derivative", "a", gradpar)
    call create_and_write_variable(gnostics, gnostics%rtype, "gbdrift", &
         dim_string([gnostics%dims%theta]), &
         "Values of gbdrift, the magnetic gradient drift ", "TBC", gbdrift)
    call create_and_write_variable(gnostics, gnostics%rtype, "gbdrift0", &
         dim_string([gnostics%dims%theta]), &
      "Values of gbdrift0, ", "TBC", gbdrift0)
    call create_and_write_variable(gnostics, gnostics%rtype, "cvdrift", &
         dim_string([gnostics%dims%theta]), &
         "Values of cvdrift, the magnetic curvature drift", "TBC", cvdrift)
    call create_and_write_variable(gnostics, gnostics%rtype, "cvdrift0", &
         dim_string([gnostics%dims%theta]), &
         "Values of cvdrift0, ", "TBC", cvdrift0)
    call create_and_write_variable(gnostics, gnostics%rtype, "gds2", &
         dim_string([gnostics%dims%theta]), &
         "Values of gds2, ", "TBC", gds2)
    call create_and_write_variable(gnostics, gnostics%rtype, "gds21", &
         "theta", &
         "Values of gds21, ", "TBC", gds21)
    call create_and_write_variable(gnostics, gnostics%rtype, "gds22", &
         dim_string([gnostics%dims%theta]), &
         "Values of gds22, ", "TBC", gds22)
    call create_and_write_variable(gnostics, gnostics%rtype, "grho", &
         dim_string([gnostics%dims%theta]), &
         "Values of grho, ", "TBC", grho)
    call create_and_write_variable(gnostics, gnostics%rtype, "jacob", &
         dim_string([gnostics%dims%theta]), &
         "Values of jacob, ", "TBC", jacob)
    call create_and_write_variable(gnostics, gnostics%rtype, "shat", "", &
         "Values of shat, the magnetic shear", "TBC", shat)
    call create_and_write_variable(gnostics, gnostics%rtype, "drhodpsi", "", &
         "Values of drhodpsi, ", "TBC", drhodpsi)
    call create_and_write_variable(gnostics, gnostics%rtype, "eps", "", &
         "Values of eps, ", "TBC", eps)
    call create_and_write_variable(gnostics, gnostics%rtype, "cdrift", &
         dim_string([gnostics%dims%theta]), &
         "Values of cdrift, ", "TBC", cdrift)
    call create_and_write_variable(gnostics, gnostics%rtype, "cdrift0", &
         dim_string([gnostics%dims%theta]), &
         "Values of cdrift0, ", "TBC", cdrift0)
    call create_and_write_variable(gnostics, gnostics%rtype, "qval", "", &
         "Values of qval, ", "TBC", qval)

    call create_and_write_variable(gnostics, gnostics%rtype, "Rplot", &
         dim_string([gnostics%dims%theta]), &
         "The major radius at the centre of the flux tube ", "a", Rplot)
    call create_and_write_variable(gnostics, gnostics%rtype, "Zplot", &
         dim_string([gnostics%dims%theta]), &
         "The height above the midplane at the centre of the flux tube ", "a", Zplot)
    call create_and_write_variable(gnostics, gnostics%rtype, "aplot", &
         dim_string([gnostics%dims%theta]), &
         "The toroidal angle at the centre of the flux tube ", "rad", aplot)
    call create_and_write_variable(gnostics, gnostics%rtype, "Rprime", &
         dim_string([gnostics%dims%theta]), &
         "d/dx of the major radius at the centre of the flux tube ", "a/rho_r", Rprime)
    call create_and_write_variable(gnostics, gnostics%rtype, "Zprime", &
         dim_string([gnostics%dims%theta]), &
         "d/dx of the height above the midplane at the centre of the flux tube ", "a/rho_r", Zprime)
    call create_and_write_variable(gnostics, gnostics%rtype, "aprime", &
         dim_string([gnostics%dims%theta]), &
         "d/dx of the toroidal angle at the centre of the flux tube ", "rad/rho_r", aprime)


    !Write data to ascii file. 
    !Should probably disable this if not write_ascii
    !Netcdf missing "shape" data
    if(proc0) then 
       call open_output_file (unit, ".g")
       write (unit,fmt="('# shape: ',a)") trim(shape)
       write (unit,fmt="('# q = ',e11.4,' drhodpsi = ',e11.4)") qval, drhodpsi
       write (unit,fmt="('# theta1             R2                  Z3               alpha4      ', &
            '       Rprime5              Zprime6           alpha_prime7     ', & 
            '      bpol8')")
       do i=-ntgrid,ntgrid
          write (unit,'(20(1x,1pg18.11))') theta(i),Rplot(i),Zplot(i),aplot(i), &
               Rprime(i),Zprime(i),aprime(i),Bpol(i)
       enddo
       call close_output_file (unit)
    endif
  end subroutine write_geometry
end module diagnostics_geometry

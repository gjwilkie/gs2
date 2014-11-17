
!> A module for writing out geometric coefficients and information
!! such as bpar, drifts and plotting coefficients
module diagnostics_geometry
  use diagnostics_create_and_write
  use simpledataio
  use diagnostics_config, only: diagnostics_type
  contains
  subroutine write_geometry(gnostics)
    use theta_grid, only: bmag, gradpar, gbdrift, gbdrift0, &
         cvdrift, cvdrift0, gds2, gds21, gds22, grho, jacob, &
         shat, drhodpsi, eps, cdrift, cdrift0, qval
    use theta_grid, only: Rplot, Zplot, aplot, Rprime, Zprime, aprime, drhodpsi, shape
    type(diagnostics_type), intent(in) :: gnostics
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

    call create_and_write_variable(gnostics, gnostics%rtype, "Rplot", "z", &
      "The major radius at the centre of the flux tube ", "a", Rplot)
    call create_and_write_variable(gnostics, gnostics%rtype, "Zplot", "z", &
      "The height above the midplane at the centre of the flux tube ", "a", Zplot)
    call create_and_write_variable(gnostics, gnostics%rtype, "aplot", "z", &
      "The toroidal angle at the centre of the flux tube ", "rad", aplot)
    call create_and_write_variable(gnostics, gnostics%rtype, "Rprime", "z", &
      "d/dx of the major radius at the centre of the flux tube ", "a/rho_r", Rprime)
    call create_and_write_variable(gnostics, gnostics%rtype, "Zprime", "z", &
      "d/dx of the height above the midplane at the centre of the flux tube ", "a/rho_r", Zprime)
    call create_and_write_variable(gnostics, gnostics%rtype, "aprime", "z", &
      "d/dx of the toroidal angle at the centre of the flux tube ", "rad/rho_r", aprime)
  end subroutine write_geometry

end module diagnostics_geometry

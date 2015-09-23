!> This module writes the parameters in the
!! normalisations namelist to the netcdf output file.
!! These parameters play no part in the simulation,
!! but will typically have been used in calculating
!! the dimensionless inputs, and will typically be used
!! to dimensionalise the outputs.
module diagnostics_normalisations
contains
  subroutine write_normalisations(gnostics)
    use normalisations, only: norms
    use diagnostics_config, only: diagnostics_type
    use diagnostics_create_and_write, only: create_and_write_variable_noread
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    call create_and_write_variable_noread(gnostics, gnostics%rtype, "mref", "", &
         "Reference mass in atomic mass units ", "a.m.u.", norms%get_value("mref"))
    call create_and_write_variable_noread(gnostics, gnostics%rtype, "zref", "", &
         "Reference charge", "proton charge", norms%get_value("zref"))
    call create_and_write_variable_noread(gnostics, gnostics%rtype, "nref", "", &
         "The reference density ", "m^-3", norms%get_value("nref"))
    call create_and_write_variable_noread(gnostics, gnostics%rtype, "tref", "", &
         "The reference temperature ", "eV", norms%get_value("tref"))
    call create_and_write_variable_noread(gnostics, gnostics%rtype, "aref", "", &
         "The reference length ", "m", norms%get_value("aref"))
    call create_and_write_variable_noread(gnostics, gnostics%rtype, "vref", "", &
         "The reference (thermal) velocity ", "m/s", norms%get_value("vref"))
    call create_and_write_variable_noread(gnostics, gnostics%rtype, "bref", "", &
         "The reference magnetic field ", "Tesla", norms%get_value("bref"))
    call create_and_write_variable_noread(gnostics, gnostics%rtype, "rhoref", "", &
         "The reference larmour radius ", "m", norms%get_value("rhoref"))
  end subroutine write_normalisations

end module diagnostics_normalisations

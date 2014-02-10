

module parameter_scan_arrays

  public :: current_scan_parameter_value
  real :: current_scan_parameter_value

  public :: write_scan_parameter
  logical :: write_scan_parameter

  public :: run_scan
  logical :: run_scan

  integer :: scan_parameter_switch
  public :: scan_parameter_switch

  integer, parameter :: scan_parameter_tprim = 1, &
                        scan_parameter_g_exb = 2

  public :: scan_parameter_tprim, &
            scan_parameter_g_exb

  integer :: scan_spec
  public :: scan_spec

  public :: hflux_tot, momflux_tot, phi2_tot
  real, dimension(:), allocatable :: hflux_tot, momflux_tot, phi2_tot

  public :: growth_rate 
  real, dimension(:), allocatable :: growth_rate 

  public :: nout
  integer :: nout

end module parameter_scan_arrays




module parameter_scan_arrays

  public :: current_scan_parameter_value
  real :: current_scan_parameter_value

  public :: write_scan_parameter
  logical :: write_scan_parameter

  public :: hflux_tot, momflux_tot, phi2_tot
  real, dimension(:), allocatable :: hflux_tot, momflux_tot, phi2_tot

  public :: growth_rate 
  real, dimension(:), allocatable :: growth_rate 

  public :: nout
  integer :: nout

end module parameter_scan_arrays


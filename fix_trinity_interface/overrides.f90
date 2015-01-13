module overrides
end module overrides

module profile_overrides
  logical :: profile_overrides_set = .false.
  integer, parameter :: otprim = 1
  integer, parameter :: ofprim = 2
  integer, parameter :: otemp = 3
  integer, parameter :: odens = 4
  integer, parameter :: ovnewk = 5
  integer, parameter :: og_exb = 6
  integer, parameter :: omach = 7
end module profile_overrides

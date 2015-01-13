module overrides
end module overrides

module gs2_profile_overrides
  integer, parameter :: otprim = 1
  integer, parameter :: ofprim = 2
  integer, parameter :: otemp = 3
  integer, parameter :: odens = 4
  integer, parameter :: ovnewk = 5
  integer, parameter :: og_exb = 6
  integer, parameter :: omach = 7
end module gs2_profile_overrides

module gs2_miller_geometry_overrides
    !if (.not. use_gs2_geo) call init_trin_geo (rhoc, qval, shat, &
         !rgeo_lcfs, rgeo_local, kap, kappri, tri, tripri, shift, betaprim)
  integer, parameter :: orhoc = 100 
  integer, parameter :: oqval = 101 
  integer, parameter :: oshat = 102 
  integer, parameter :: orgeo_lcfs = 104
  integer, parameter :: orgeo_local = 105
  integer, parameter :: okap = 106
  integer, parameter :: okappri = 107
  integer, parameter :: otri = 108
  integer, parameter :: otripri = 109
  integer, parameter :: oshift = 110
  integer, parameter :: obetaprim = 111
  
end module gs2_miller_geometry_overrides


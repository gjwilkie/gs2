module dist_fn_arrays

  ! dist fn
  complex, dimension (:,:,:), allocatable :: g, gnew, gold
  ! (-ntgrid:ntgrid,2, -g-layout-)

  real, dimension (:,:,:), allocatable :: vpa, vpac
  ! (-ntgrid:ntgrid,2, -g-layout-)

  real, dimension (:,:), allocatable :: vperp2
  ! (-ntgrid:ntgrid, -g-layout-)

  real, dimension (:,:,:), allocatable :: vpar
  ! (-ntgrid:ntgrid,2, -g-layout-)

  integer, dimension (:), allocatable :: ittp
  ! (-ntgrid:ntgrid)

  real, dimension (:,:), allocatable :: aj0, aj1, aj0f, aj1f
  ! (-ntgrid:ntgrid, -g-layout-)

  ! fieldeq
  complex, dimension (:,:,:), allocatable :: apar_ext
  real, dimension (:,:,:), allocatable :: kperp2
  ! (-ntgrid:ntgrid,ntheta0,naky) replicated

  ! collisional diagnostic of heating rate
  complex, dimension (:,:,:,:), allocatable :: c_rate
  ! (-ntgrid:ntgrid,ntheta0,naky) replicated

end module dist_fn_arrays

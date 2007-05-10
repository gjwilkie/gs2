module dist_fn_arrays

  ! dist fn
  complex, dimension (:,:,:), allocatable :: g, gnew
  ! (-ntgrid:ntgrid,2, -g-layout-)

  real, dimension (:,:,:), allocatable :: vpa, vpac
  ! (-ntgrid:ntgrid,2, -g-layout-)

  real, dimension (:,:), allocatable :: vperp2
  ! (-ntgrid:ntgrid, -g-layout-)

  real, dimension (:,:,:), allocatable :: vpar
  ! (-ntgrid:ntgrid,2, -g-layout-)

  real, dimension (:,:), allocatable :: aj0, aj1, aj0f, aj1f
  ! (-ntgrid:ntgrid, -g-layout-)

  ! fieldeq
  real, dimension (:,:,:), allocatable :: kperp2
  ! (-ntgrid:ntgrid,naky,ntheta0) replicated

end module dist_fn_arrays

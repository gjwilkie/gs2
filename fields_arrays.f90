module fields_arrays
  implicit none

  complex, dimension (:,:,:), allocatable :: phi,    apar,    aperp
  complex, dimension (:,:,:), allocatable :: phinew, aparnew, aperpnew
  ! (-ntgrid:ntgrid,naky,ntheta0) replicated

end module fields_arrays

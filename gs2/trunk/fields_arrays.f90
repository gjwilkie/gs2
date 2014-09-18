module fields_arrays
  implicit none

  !Main fields
  complex, dimension (:,:,:), allocatable :: phi,    apar,    bpar
  complex, dimension (:,:,:), allocatable :: phinew, aparnew, bparnew
  !For antenna
  complex, dimension (:,:,:), allocatable :: apar_ext
  ! (-ntgrid:ntgrid,ntheta0,naky) replicated

  !Timing data
  real, save :: time_field(2)=0.
end module fields_arrays

module fields_arrays
  implicit none

  private

  public :: phi, apar, bpar, phinew, aparnew, bparnew
  public :: gf_phi, gf_apar, gf_bpar, gf_phinew, gf_aparnew, gf_bparnew
  public :: apar_ext, time_field, response_file

  !Main fields
  complex, dimension (:,:,:), allocatable :: phi,    apar,    bpar
  complex, dimension (:,:,:), allocatable :: phinew, aparnew, bparnew
  complex, dimension (:,:,:), allocatable :: gf_phi, gf_apar, gf_bpar
  complex, dimension (:,:,:), allocatable :: gf_phinew, gf_aparnew, gf_bparnew
  !For antenna
  complex, dimension (:,:,:), allocatable :: apar_ext
  ! (-ntgrid:ntgrid,ntheta0,naky) replicated

  !Timing data
  real :: time_field(2)=0.
  
  !For response data
  character(len=300) :: response_file
end module fields_arrays

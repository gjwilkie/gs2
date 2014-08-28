module fields_arrays
  implicit none

  !Main fields
  complex, dimension (:,:,:), allocatable :: phi,    apar,    bpar
  complex, dimension (:,:,:), allocatable :: phinew, aparnew, bparnew
  !For antenna
  complex, dimension (:,:,:), allocatable :: apar_ext
  ! (-ntgrid:ntgrid,ntheta0,naky) replicated

  integer, save :: nidx

  real, save :: time_field(2)=0.

! Within each supercell, there are are N_class primary cells.  Each 
! has (2*ntgrid+1)*nfield points.
  type dcell_type
     complex, dimension(:), pointer :: supercell => null()
  end type dcell_type

! Within each class, there may be multiple supercells.

! The number of supercells in each class is M_class.

! When aminv is laid out over PEs, the supercells of each class 
! are distributed -- thus, "dcells"

  type :: field_matrix_type
     type(dcell_type), dimension (:), pointer :: dcell => null()
  end type field_matrix_type

! There may be supercells of different sizes or "classes".  
  type (field_matrix_type), dimension (:), allocatable :: aminv
end module fields_arrays

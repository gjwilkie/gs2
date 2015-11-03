module fields_arrays
  implicit none

  complex, dimension (:,:,:), allocatable :: phi,    apar,    bpar
  complex, dimension (:,:,:), allocatable :: phi_ext, apar_ext
  complex, dimension (:,:,:), allocatable :: phinew, aparnew, bparnew
  complex, dimension (:,:,:), allocatable :: phitmp, apartmp, bpartmp
  complex, dimension (:,:,:), allocatable :: phitmp1, apartmp1, bpartmp1
  ! (-ntgrid:ntgrid,ntheta0,naky) replicated

!!  complex, dimension (:,:), allocatable :: aminv
!!  ! (nidx, -*- f_lo -*-)

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

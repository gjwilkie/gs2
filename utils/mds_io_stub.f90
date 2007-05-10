!-----------------------------------------------------------------------
!     file mdsplus_io.f.
!     performs MDSplus tree io
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     code organization.
!-----------------------------------------------------------------------
!     0. mdsplus_io_mod.
!     1. GetMdsLogical.
!     2. GetMdsInteger.
!     3. GetMdsDouble.
!     4. GetMdsDouble2DArray.
!     5. GetMdsDouble3DArray.
!     6. GetMdsText.
!     7. CheckMds.
!     8. GetMdsErrorText.
!     9. GetMdsReal
!-----------------------------------------------------------------------
!     subprogram 0. mdsplus_io_mod.
!     module declarations.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.
!-----------------------------------------------------------------------
MODULE mdsio
  use prec
  IMPLICIT NONE
  Private
  public :: mds_read
  public :: mds_write
  public :: CheckMDS
  public :: GetMdsErrorText


  interface mds_read
     module procedure mds_r_log_0
     module procedure mds_r_char_1
     module procedure mds_r_int_0
     module procedure mds_r_int_1
     module procedure mds_r_int_2
     module procedure mds_r_int_3
     module procedure mds_r_int_4
     module procedure mds_r_int_5
     module procedure mds_r_int_6
     module procedure mds_r_int_7
     module procedure mds_r_real_0
     module procedure mds_r_real_1
     module procedure mds_r_real_2
     module procedure mds_r_real_3
     module procedure mds_r_real_4
     module procedure mds_r_real_5
     module procedure mds_r_real_6
     module procedure mds_r_real_7
     module procedure mds_r_doub_0
     module procedure mds_r_doub_1
     module procedure mds_r_doub_2
     module procedure mds_r_doub_3
     module procedure mds_r_doub_4
     module procedure mds_r_doub_5
     module procedure mds_r_doub_6
     module procedure mds_r_doub_7
     module procedure mds_r_cmplx_0
     module procedure mds_r_cmplx_1
     module procedure mds_r_cmplx_2
     module procedure mds_r_cmplx_3
     module procedure mds_r_cmplx_4
     module procedure mds_r_cmplx_5
     module procedure mds_r_cmplx_6
     module procedure mds_r_cmplx_7
     module procedure mds_r_dcmplx_0
     module procedure mds_r_dcmplx_1
     module procedure mds_r_dcmplx_2
     module procedure mds_r_dcmplx_3
     module procedure mds_r_dcmplx_4
     module procedure mds_r_dcmplx_5
     module procedure mds_r_dcmplx_6
     module procedure mds_r_dcmplx_7
  end interface
  
  interface mds_write
     module procedure mds_w_log_0
     module procedure mds_w_char_1
     module procedure mds_w_int_0
     module procedure mds_w_int_1
     module procedure mds_w_int_2
     module procedure mds_w_int_3
     module procedure mds_w_int_4
     module procedure mds_w_int_5
     module procedure mds_w_int_6
     module procedure mds_w_int_7
     module procedure mds_w_real_0
     module procedure mds_w_real_1
     module procedure mds_w_real_2
     module procedure mds_w_real_3
     module procedure mds_w_real_4
     module procedure mds_w_real_5
     module procedure mds_w_real_6
     module procedure mds_w_real_7
     module procedure mds_w_doub_0
     module procedure mds_w_doub_1
     module procedure mds_w_doub_2
     module procedure mds_w_doub_3
     module procedure mds_w_doub_4
     module procedure mds_w_doub_5
     module procedure mds_w_doub_6
     module procedure mds_w_doub_7
     module procedure mds_w_cmplx_0
     module procedure mds_w_cmplx_1
     module procedure mds_w_cmplx_2
     module procedure mds_w_cmplx_3
     module procedure mds_w_cmplx_4
     module procedure mds_w_cmplx_5
     module procedure mds_w_cmplx_6
     module procedure mds_w_cmplx_7
     module procedure mds_w_dcmplx_0
     module procedure mds_w_dcmplx_1
     module procedure mds_w_dcmplx_2
     module procedure mds_w_dcmplx_3
     module procedure mds_w_dcmplx_4
     module procedure mds_w_dcmplx_5
     module procedure mds_w_dcmplx_6
     module procedure mds_w_dcmplx_7
  end interface

CONTAINS
!---------------------------------------------
! mdsplus internal routines
!---------------------------------------------

  SUBROUTINE CheckMds(status,msg)
    LOGICAL, INTENT(IN) :: status
    CHARACTER(*), INTENT(IN) :: msg
  END SUBROUTINE CheckMds
  
  SUBROUTINE GetMdsErrorText(status,text)
    LOGICAL, INTENT(IN) :: status
    CHARACTER(*) :: text
  END SUBROUTINE GetMdsErrorText

!--------------------------------------------
! mds_read subroutines
!--------------------------------------------
 
  SUBROUTINE mds_r_log_0(name,value)
    CHARACTER(*), INTENT(IN) :: name
    LOGICAL, intent(out) :: value
  END SUBROUTINE mds_r_log_0
  
  SUBROUTINE mds_r_char_1(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Character*(*) , intent(out):: value
  END SUBROUTINE mds_r_char_1

  SUBROUTINE mds_r_int_0(name,value)
    CHARACTER(*), INTENT(IN) :: name
    INTEGER , intent(out) :: value
  END SUBROUTINE mds_r_int_0

  SUBROUTINE mds_r_int_1(name,value)
    CHARACTER(*), INTENT(IN) :: name 
    INTEGER , dimension(:), intent(out) :: value
  END SUBROUTINE mds_r_int_1

  SUBROUTINE mds_r_int_2(name,value)
    CHARACTER(*), INTENT(IN) :: name 
    INTEGER , dimension(:,:), intent(out) :: value
  END SUBROUTINE mds_r_int_2

  SUBROUTINE mds_r_int_3(name,value)
    CHARACTER(*), INTENT(IN) :: name 
    INTEGER , dimension(:,:,:), intent(out) :: value
  END SUBROUTINE mds_r_int_3

  SUBROUTINE mds_r_int_4(name,value)
    CHARACTER(*), INTENT(IN) :: name 
    INTEGER , dimension(:,:,:,:), intent(out) :: value
  END SUBROUTINE mds_r_int_4

  SUBROUTINE mds_r_int_5(name,value)
    CHARACTER(*), INTENT(IN) :: name 
    INTEGER , dimension(:,:,:,:,:), intent(out) :: value
  END SUBROUTINE mds_r_int_5

  SUBROUTINE mds_r_int_6(name,value)
    CHARACTER(*), INTENT(IN) :: name 
    INTEGER , dimension(:,:,:,:,:,:), intent(out) :: value
  END SUBROUTINE mds_r_int_6

  SUBROUTINE mds_r_int_7(name,value)
    CHARACTER(*), INTENT(IN) :: name 
    INTEGER , dimension(:,:,:,:,:,:,:), intent(out) :: value
  END SUBROUTINE mds_r_int_7

  SUBROUTINE mds_r_real_0(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=SP) , intent(out) :: value
  END SUBROUTINE mds_r_real_0

  SUBROUTINE mds_r_real_1(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=SP) , dimension(:), intent(out) :: value
  END SUBROUTINE mds_r_real_1

  SUBROUTINE mds_r_real_2(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=SP) , dimension(:,:), intent(out) :: value
  END SUBROUTINE mds_r_real_2

  SUBROUTINE mds_r_real_3(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=SP) , dimension(:,:,:), intent(out) :: value
  END SUBROUTINE mds_r_real_3

  SUBROUTINE mds_r_real_4(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=SP) , dimension(:,:,:,:), intent(out) :: value
  END SUBROUTINE mds_r_real_4

  SUBROUTINE mds_r_real_5(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=SP) , dimension(:,:,:,:,:), intent(out) :: value
  END SUBROUTINE mds_r_real_5

  SUBROUTINE mds_r_real_6(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=SP) , dimension(:,:,:,:,:,:), intent(out) :: value
  END SUBROUTINE mds_r_real_6

  SUBROUTINE mds_r_real_7(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=SP) , dimension(:,:,:,:,:,:,:), intent(out) :: value
  END SUBROUTINE mds_r_real_7

  SUBROUTINE mds_r_doub_0(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=DP) , intent(out):: value
  END SUBROUTINE mds_r_doub_0

  SUBROUTINE mds_r_doub_1(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=DP), DIMENSION(:), intent(out) :: value
  END SUBROUTINE mds_r_doub_1

  SUBROUTINE mds_r_doub_2(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=DP), DIMENSION(:,:), intent(out) :: value
  END SUBROUTINE mds_r_doub_2

  SUBROUTINE mds_r_doub_3(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=DP), DIMENSION(:,:,:), intent(out) :: value
  END SUBROUTINE mds_r_doub_3

  SUBROUTINE mds_r_doub_4(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=DP), DIMENSION(:,:,:,:), intent(out) :: value
  END SUBROUTINE mds_r_doub_4

  SUBROUTINE mds_r_doub_5(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=DP), DIMENSION(:,:,:,:,:), intent(out) :: value
  END SUBROUTINE mds_r_doub_5

  SUBROUTINE mds_r_doub_6(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=DP), DIMENSION(:,:,:,:,:,:), intent(out) :: value
  END SUBROUTINE mds_r_doub_6

  SUBROUTINE mds_r_doub_7(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=DP), DIMENSION(:,:,:,:,:,:,:), intent(out) :: value
  END SUBROUTINE mds_r_doub_7

  SUBROUTINE mds_r_cmplx_0(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=SP), intent(out) :: value
  END SUBROUTINE mds_r_cmplx_0

  SUBROUTINE mds_r_cmplx_1(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=SP), DIMENSION(:), intent(out) :: value
  END SUBROUTINE mds_r_cmplx_1

  SUBROUTINE mds_r_cmplx_2(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=SP), DIMENSION(:,:), intent(out) :: value
  END SUBROUTINE mds_r_cmplx_2

  SUBROUTINE mds_r_cmplx_3(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=SP), DIMENSION(:,:,:), intent(out) :: value
  END SUBROUTINE mds_r_cmplx_3

  SUBROUTINE mds_r_cmplx_4(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=SP), DIMENSION(:,:,:,:), intent(out) :: value
  END SUBROUTINE mds_r_cmplx_4

  SUBROUTINE mds_r_cmplx_5(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=SP), DIMENSION(:,:,:,:,:), intent(out) :: value
  END SUBROUTINE mds_r_cmplx_5

  SUBROUTINE mds_r_cmplx_6(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=SP), DIMENSION(:,:,:,:,:,:), intent(out) :: value
  END SUBROUTINE mds_r_cmplx_6

  SUBROUTINE mds_r_cmplx_7(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=SP), DIMENSION(:,:,:,:,:,:,:), intent(out) :: value
  END SUBROUTINE mds_r_cmplx_7

  SUBROUTINE mds_r_dcmplx_0(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=DP), intent(out) :: value
  END SUBROUTINE mds_r_dcmplx_0

  SUBROUTINE mds_r_dcmplx_1(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=DP), DIMENSION(:), intent(out) :: value
  END SUBROUTINE mds_r_dcmplx_1

  SUBROUTINE mds_r_dcmplx_2(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=DP), DIMENSION(:,:), intent(out) :: value
  END SUBROUTINE mds_r_dcmplx_2

  SUBROUTINE mds_r_dcmplx_3(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=DP), DIMENSION(:,:,:), intent(out) :: value
  END SUBROUTINE mds_r_dcmplx_3

  SUBROUTINE mds_r_dcmplx_4(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=DP), DIMENSION(:,:,:,:), intent(out) :: value
  END SUBROUTINE mds_r_dcmplx_4

  SUBROUTINE mds_r_dcmplx_5(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=DP), DIMENSION(:,:,:,:,:), intent(out) :: value
  END SUBROUTINE mds_r_dcmplx_5

  SUBROUTINE mds_r_dcmplx_6(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=DP), DIMENSION(:,:,:,:,:,:), intent(out) :: value
  END SUBROUTINE mds_r_dcmplx_6

  SUBROUTINE mds_r_dcmplx_7(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=DP), DIMENSION(:,:,:,:,:,:,:), intent(out) :: value
  END SUBROUTINE mds_r_dcmplx_7

  !--------------------------------------------
  ! mds_write subroutines
  !--------------------------------------------

  SUBROUTINE mds_w_log_0(name,value)
    CHARACTER(*), INTENT(IN) :: name
    LOGICAL :: value
  END SUBROUTINE mds_w_log_0

  SUBROUTINE mds_w_char_1(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Character*(*), intent(in) :: value
  END SUBROUTINE mds_w_char_1

  SUBROUTINE mds_w_int_0(name,value)
    CHARACTER(*), INTENT(IN) :: name
    INTEGER , intent(in) :: value
  END SUBROUTINE mds_w_int_0

  SUBROUTINE mds_w_int_1(name,value)
    CHARACTER(*), INTENT(IN) :: name
    INTEGER , DIMENSION(:), intent(in) :: value
  END SUBROUTINE mds_w_int_1

  SUBROUTINE mds_w_int_2(name,value)
    CHARACTER(*), INTENT(IN) :: name 
    INTEGER , dimension(:,:), intent(in) :: value
  END SUBROUTINE mds_w_int_2

  SUBROUTINE mds_w_int_3(name,value)
    CHARACTER(*), INTENT(IN) :: name 
    INTEGER , dimension(:,:,:), intent(in) :: value
  END SUBROUTINE mds_w_int_3

  SUBROUTINE mds_w_int_4(name,value)
    CHARACTER(*), INTENT(IN) :: name 
    INTEGER , dimension(:,:,:,:), intent(in) :: value
  END SUBROUTINE mds_w_int_4

  SUBROUTINE mds_w_int_5(name,value)
    CHARACTER(*), INTENT(IN) :: name 
    INTEGER , dimension(:,:,:,:,:), intent(in) :: value
  END SUBROUTINE mds_w_int_5

  SUBROUTINE mds_w_int_6(name,value)
    CHARACTER(*), INTENT(IN) :: name 
    INTEGER , dimension(:,:,:,:,:,:), intent(in) :: value
  END SUBROUTINE mds_w_int_6

  SUBROUTINE mds_w_int_7(name,value)
    CHARACTER(*), INTENT(IN) :: name 
    INTEGER , dimension(:,:,:,:,:,:,:), intent(in) :: value
  END SUBROUTINE mds_w_int_7

  SUBROUTINE mds_w_real_0(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=SP), intent(in) :: value
  END SUBROUTINE mds_w_real_0

  SUBROUTINE mds_w_real_1(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=SP), DIMENSION(:), intent(in) :: value
  END SUBROUTINE mds_w_real_1

  SUBROUTINE mds_w_real_2(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=SP) , dimension(:,:), intent(in) :: value
  END SUBROUTINE mds_w_real_2

  SUBROUTINE mds_w_real_3(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=SP) , dimension(:,:,:), intent(in) :: value
  END SUBROUTINE mds_w_real_3

  SUBROUTINE mds_w_real_4(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=SP) , dimension(:,:,:,:), intent(in) :: value
  END SUBROUTINE mds_w_real_4

  SUBROUTINE mds_w_real_5(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=SP) , dimension(:,:,:,:,:), intent(in) :: value
  END SUBROUTINE mds_w_real_5

  SUBROUTINE mds_w_real_6(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=SP) , dimension(:,:,:,:,:,:), intent(in) :: value
  END SUBROUTINE mds_w_real_6

  SUBROUTINE mds_w_real_7(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=SP) , dimension(:,:,:,:,:,:,:), intent(in) :: value
  END SUBROUTINE mds_w_real_7

  SUBROUTINE mds_w_doub_0(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=DP), intent(in) :: value
  END SUBROUTINE mds_w_doub_0

  SUBROUTINE mds_w_doub_1(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=DP), DIMENSION(:), intent(in) :: value
  END SUBROUTINE mds_w_doub_1

  SUBROUTINE mds_w_doub_2(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=DP), DIMENSION(:,:), intent(in) :: value
  END SUBROUTINE mds_w_doub_2

  SUBROUTINE mds_w_doub_3(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=DP), DIMENSION(:,:,:), intent(in) :: value
  END SUBROUTINE mds_w_doub_3

  SUBROUTINE mds_w_doub_4(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=DP), DIMENSION(:,:,:,:), intent(in) :: value
  END SUBROUTINE mds_w_doub_4

  SUBROUTINE mds_w_doub_5(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=DP), DIMENSION(:,:,:,:,:), intent(in) :: value
  END SUBROUTINE mds_w_doub_5

  SUBROUTINE mds_w_doub_6(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=DP), DIMENSION(:,:,:,:,:,:), intent(in) :: value
  END SUBROUTINE mds_w_doub_6

  SUBROUTINE mds_w_doub_7(name,value)
    CHARACTER(*), INTENT(IN) :: name
    REAL(KIND=DP), DIMENSION(:,:,:,:,:,:,:), intent(in) :: value
  END SUBROUTINE mds_w_doub_7

  SUBROUTINE mds_w_cmplx_0(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=SP), intent(in) :: value
  END SUBROUTINE mds_w_cmplx_0

  SUBROUTINE mds_w_cmplx_1(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=SP), DIMENSION(:), intent(in) :: value
  END SUBROUTINE mds_w_cmplx_1

  SUBROUTINE mds_w_cmplx_2(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=SP), DIMENSION(:,:), intent(in) :: value
  END SUBROUTINE mds_w_cmplx_2

  SUBROUTINE mds_w_cmplx_3(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=SP), DIMENSION(:,:,:), intent(in) :: value
  END SUBROUTINE mds_w_cmplx_3

  SUBROUTINE mds_w_cmplx_4(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=SP), DIMENSION(:,:,:,:), intent(in) :: value
  END SUBROUTINE mds_w_cmplx_4

  SUBROUTINE mds_w_cmplx_5(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=SP), DIMENSION(:,:,:,:,:), intent(in) :: value
  END SUBROUTINE mds_w_cmplx_5

  SUBROUTINE mds_w_cmplx_6(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=SP), DIMENSION(:,:,:,:,:,:), intent(in) :: value
  END SUBROUTINE mds_w_cmplx_6

  SUBROUTINE mds_w_cmplx_7(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=SP), DIMENSION(:,:,:,:,:,:,:), intent(in) :: value
  END SUBROUTINE mds_w_cmplx_7

  SUBROUTINE mds_w_dcmplx_0(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=DP), intent(in) :: value
  END SUBROUTINE mds_w_dcmplx_0

  SUBROUTINE mds_w_dcmplx_1(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=DP), DIMENSION(:), intent(in) :: value
  END SUBROUTINE mds_w_dcmplx_1

  SUBROUTINE mds_w_dcmplx_2(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=DP), DIMENSION(:,:), intent(in) :: value
  END SUBROUTINE mds_w_dcmplx_2

  SUBROUTINE mds_w_dcmplx_3(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=DP), DIMENSION(:,:,:), intent(in) :: value
  END SUBROUTINE mds_w_dcmplx_3

  SUBROUTINE mds_w_dcmplx_4(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=DP), DIMENSION(:,:,:,:), intent(in) :: value
  END SUBROUTINE mds_w_dcmplx_4

  SUBROUTINE mds_w_dcmplx_5(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=DP), DIMENSION(:,:,:,:,:), intent(in) :: value
  END SUBROUTINE mds_w_dcmplx_5

  SUBROUTINE mds_w_dcmplx_6(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=DP), DIMENSION(:,:,:,:,:,:), intent(in) :: value
  END SUBROUTINE mds_w_dcmplx_6

  SUBROUTINE mds_w_dcmplx_7(name,value)
    CHARACTER(*), INTENT(IN) :: name
    Complex(KIND=DP), DIMENSION(:,:,:,:,:,:,:), intent(in) :: value
  END SUBROUTINE mds_w_dcmplx_7

END MODULE mdsio


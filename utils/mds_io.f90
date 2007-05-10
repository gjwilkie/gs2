!-----------------------------------------------------------------------
!     file mdsplus_io.f.
!     performs mdsplus tree io
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     code organization.
!-----------------------------------------------------------------------
!     0. mdsplus_io_mod.
!     1. getmdslogical.
!     2. getmdsinteger.
!     3. getmdsdouble.
!     4. getmdsdouble2darray.
!     5. getmdsdouble3darray.
!     6. getmdstext.
!     7. checkmds.
!     8. getmdserrortext.
!     9. getmdsreal
!-----------------------------------------------------------------------
!     subprogram 0. mdsplus_io_mod.
!     module declarations.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.
!-----------------------------------------------------------------------
module mdsio
  use prec
  implicit none
  private
  external mdssetdefault
  logical  mdssetdefault
  external mdsvalue, mdsput
  logical  mdsvalue, mdsput
  public :: mds_read
  public :: mds_write
  public :: checkmds
  public :: getmdserrortext


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

contains
!---------------------------------------------
! mdsplus internal routines
!---------------------------------------------

  subroutine checkmds(status,msg)
    logical, intent(in) :: status
    character(*), intent(in) :: msg
    character(512) text
    integer istat
    logical lstat
    equivalence(istat,lstat)      
    lstat = status
    if (iand(istat,1) .eq. 0) then
       call getmdserrortext(status,text)
       write(*,*) text
       write(*,*) msg
       !       call program_stop(msg//", "//trim(text))
    endif
  end subroutine checkmds

  subroutine getmdserrortext(status,text)
    logical, intent(in) :: status
    character(*) :: text
    logical loc_status
    integer istat
    equivalence (loc_status,istat)
    integer length
    integer :: descr
    loc_status = mdsvalue("getmsg($)",descr(8,status,0),&
         &descr(14,text,0,len(text),0),0,length)
    if (iand(istat,1) .eq. 0) then
       loc_status = status
       write(text,*) "error status = ",istat
    endif
    return
  end subroutine getmdserrortext


!--------------------------------------------
! mds_read subroutines
!--------------------------------------------

  subroutine mds_r_log_0(name,value)
    character(*), intent(in) :: name
    logical, intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    status = mdsvalue(name//char(0),descr(8,value,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_log_0

  subroutine mds_r_char_1(name,value)
    character(*), intent(in) :: name
    character*(*) , intent(out):: value
    logical :: status
    integer :: clen 
    integer :: descr
    status = mdsvalue(name//char(0),&
         &descr(14,value,0,len(value)),0,clen)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_char_1

  subroutine mds_r_int_0(name,value)
    character(*), intent(in) :: name
    integer , intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    status = mdsvalue(name//char(0),descr(8,value,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_int_0

  subroutine mds_r_int_1(name,value)
    character(*), intent(in) :: name 
    integer , dimension(:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1
    dim1 = size(value,1)
    status = mdsvalue(name//char(0),descr(8,value,dim1,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_int_1

  subroutine mds_r_int_2(name,value)
    character(*), intent(in) :: name 
    integer , dimension(:,:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2
    dim1 = size(value,1)
    dim2 = size(value,2)
    status = mdsvalue(name//char(0),descr(8,value,dim1,dim2,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_int_2

  subroutine mds_r_int_3(name,value)
    character(*), intent(in) :: name 
    integer , dimension(:,:,:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    status = mdsvalue(name//char(0),descr(8,value,dim1,dim2,dim3,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_int_3

  subroutine mds_r_int_4(name,value)
    character(*), intent(in) :: name 
    integer , dimension(:,:,:,:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    status = mdsvalue(name//char(0),&
         &descr(8,value,dim1,dim2,dim3,dim4,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_int_4

  subroutine mds_r_int_5(name,value)
    character(*), intent(in) :: name 
    integer , dimension(:,:,:,:,:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    status = mdsvalue(name//char(0),&
         &descr(8,value,dim1,dim2,dim3,dim4,dim5,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_int_5

  subroutine mds_r_int_6(name,value)
    character(*), intent(in) :: name 
    integer , dimension(:,:,:,:,:,:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    status = mdsvalue(name//char(0),&
         &descr(8,value,dim1,dim2,dim3,dim4,dim5,dim6,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_int_6

  subroutine mds_r_int_7(name,value)
    character(*), intent(in) :: name 
    integer , dimension(:,:,:,:,:,:,:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6,dim7
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    dim7 = size(value,7)
    status = mdsvalue(name//char(0),&
         &descr(8,value,dim1,dim2,dim3,dim4,dim5,dim6,dim7,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_int_7

  subroutine mds_r_real_0(name,value)
    character(*), intent(in) :: name
    real(kind=sp) , intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    status = mdsvalue(name//char(0),descr(10,value,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_real_0

  subroutine mds_r_real_1(name,value)
    character(*), intent(in) :: name
    real(kind=sp) , dimension(:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1
    dim1 = size(value,1)
    status = mdsvalue(name//char(0),descr(10,value,dim1,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_real_1

  subroutine mds_r_real_2(name,value)
    character(*), intent(in) :: name
    real(kind=sp) , dimension(:,:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2
    dim1 = size(value,1)
    dim2 = size(value,2)
    status = mdsvalue(name//char(0),&
         &descr(10,value,dim1,dim2,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_real_2

  subroutine mds_r_real_3(name,value)
    character(*), intent(in) :: name
    real(kind=sp) , dimension(:,:,:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    status = mdsvalue(name//char(0),&
         &descr(10,value,dim1,dim2,dim3,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_real_3

  subroutine mds_r_real_4(name,value)
    character(*), intent(in) :: name
    real(kind=sp) , dimension(:,:,:,:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    status = mdsvalue(name//char(0),&
         &descr(10,value,dim1,dim2,dim3,dim4,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_real_4

  subroutine mds_r_real_5(name,value)
    character(*), intent(in) :: name
    real(kind=sp) , dimension(:,:,:,:,:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    status = mdsvalue(name//char(0),&
         &descr(10,value,dim1,dim2,dim3,dim4,dim5,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_real_5

  subroutine mds_r_real_6(name,value)
    character(*), intent(in) :: name
    real(kind=sp) , dimension(:,:,:,:,:,:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    status = mdsvalue(name//char(0),&
         &descr(10,value,dim1,dim2,dim3,dim4,dim5,dim6,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_real_6

  subroutine mds_r_real_7(name,value)
    character(*), intent(in) :: name
    real(kind=sp) , dimension(:,:,:,:,:,:,:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6,dim7
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    dim7 = size(value,7)
    status = mdsvalue(name//char(0),&
         &descr(10,value,dim1,dim2,dim3,dim4,dim5,dim6,dim7,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_real_7

  subroutine mds_r_doub_0(name,value)
    character(*), intent(in) :: name
    real(kind=dp) , intent(out):: value
    logical :: status
    integer :: len
    integer :: descr
    status = mdsvalue(name//char(0),descr(11,value,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_doub_0

  subroutine mds_r_doub_1(name,value)
    character(*), intent(in) :: name
    real(kind=dp), dimension(:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1
    dim1 = size(value,1)
    status = mdsvalue(name//char(0),descr(11,value,dim1,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_doub_1

  subroutine mds_r_doub_2(name,value)
    character(*), intent(in) :: name
    real(kind=dp), dimension(:,:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2
    dim1 = size(value,1)
    dim2 = size(value,2)
    status = mdsvalue(name//char(0),&
         &descr(11,value,dim1,dim2,0,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_doub_2

  subroutine mds_r_doub_3(name,value)
    character(*), intent(in) :: name
    real(kind=dp), dimension(:,:,:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    status = mdsvalue(name//char(0),&
         &descr(11,value,dim1,dim2,dim3,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_doub_3

  subroutine mds_r_doub_4(name,value)
    character(*), intent(in) :: name
    real(kind=dp), dimension(:,:,:,:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    status = mdsvalue(name//char(0),&
         &descr(11,value,dim1,dim2,dim3,dim4,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_doub_4

  subroutine mds_r_doub_5(name,value)
    character(*), intent(in) :: name
    real(kind=dp), dimension(:,:,:,:,:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    status = mdsvalue(name//char(0),&
         &descr(11,value,dim1,dim2,dim3,dim4,dim5,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_doub_5

  subroutine mds_r_doub_6(name,value)
    character(*), intent(in) :: name
    real(kind=dp), dimension(:,:,:,:,:,:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    status = mdsvalue(name//char(0),&
         &descr(11,value,dim1,dim2,dim3,dim4,dim5,dim6,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_doub_6

  subroutine mds_r_doub_7(name,value)
    character(*), intent(in) :: name
    real(kind=dp), dimension(:,:,:,:,:,:,:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6,dim7
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    dim7 = size(value,7)
    status = mdsvalue(name//char(0),&
         &descr(11,value,dim1,dim2,dim3,dim4,dim5,dim6,dim7,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_doub_7

  subroutine mds_r_cmplx_0(name,value)
    character(*), intent(in) :: name
    complex(kind=sp), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    status = mdsvalue(name//char(0),descr(12,value,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_cmplx_0

  subroutine mds_r_cmplx_1(name,value)
    character(*), intent(in) :: name
    complex(kind=sp), dimension(:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1
    dim1 = size(value,1)
    status = mdsvalue(name//char(0),descr(12,value,dim1,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_cmplx_1

  subroutine mds_r_cmplx_2(name,value)
    character(*), intent(in) :: name
    complex(kind=sp), dimension(:,:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2
    dim1 = size(value,1)
    dim2 = size(value,2)
    status = mdsvalue(name//char(0),&
         &descr(12,value,dim1,dim2,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_cmplx_2

  subroutine mds_r_cmplx_3(name,value)
    character(*), intent(in) :: name
    complex(kind=sp), dimension(:,:,:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    status = mdsvalue(name//char(0),&
         &descr(12,value,dim1,dim2,dim3,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_cmplx_3

  subroutine mds_r_cmplx_4(name,value)
    character(*), intent(in) :: name
    complex(kind=sp), dimension(:,:,:,:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    status = mdsvalue(name//char(0),&
         &descr(12,value,dim1,dim2,dim3,dim4,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_cmplx_4

  subroutine mds_r_cmplx_5(name,value)
    character(*), intent(in) :: name
    complex(kind=sp), dimension(:,:,:,:,:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    status = mdsvalue(name//char(0),&
         &descr(12,value,dim1,dim2,dim3,dim4,dim5,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_cmplx_5

  subroutine mds_r_cmplx_6(name,value)
    character(*), intent(in) :: name
    complex(kind=sp), dimension(:,:,:,:,:,:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    status = mdsvalue(name//char(0),&
         &descr(12,value,dim1,dim2,dim3,dim4,dim5,dim6,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_cmplx_6

  subroutine mds_r_cmplx_7(name,value)
    character(*), intent(in) :: name
    complex(kind=sp), dimension(:,:,:,:,:,:,:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6,dim7
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    dim7 = size(value,7)
    status = mdsvalue(name//char(0),&
         &descr(12,value,dim1,dim2,dim3,dim4,dim5,dim6,dim7,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_cmplx_7

  subroutine mds_r_dcmplx_0(name,value)
    character(*), intent(in) :: name
    complex(kind=dp), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    status = mdsvalue(name//char(0),descr(13,value,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_dcmplx_0

  subroutine mds_r_dcmplx_1(name,value)
    character(*), intent(in) :: name
    complex(kind=dp), dimension(:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1
    dim1 = size(value,1)
    status = mdsvalue(name//char(0),descr(13,value,dim1,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_dcmplx_1

  subroutine mds_r_dcmplx_2(name,value)
    character(*), intent(in) :: name
    complex(kind=dp), dimension(:,:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2
    dim1 = size(value,1)
    dim2 = size(value,2)
    status = mdsvalue(name//char(0),&
         &descr(13,value,dim1,dim2,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_dcmplx_2

  subroutine mds_r_dcmplx_3(name,value)
    character(*), intent(in) :: name
    complex(kind=dp), dimension(:,:,:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    status = mdsvalue(name//char(0),&
         &descr(13,value,dim1,dim2,dim3,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_dcmplx_3

  subroutine mds_r_dcmplx_4(name,value)
    character(*), intent(in) :: name
    complex(kind=dp), dimension(:,:,:,:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    status = mdsvalue(name//char(0),&
         &descr(13,value,dim1,dim2,dim3,dim4,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_dcmplx_4

  subroutine mds_r_dcmplx_5(name,value)
    character(*), intent(in) :: name
    complex(kind=dp), dimension(:,:,:,:,:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    status = mdsvalue(name//char(0),&
         &descr(13,value,dim1,dim2,dim3,dim4,dim5,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_dcmplx_5

  subroutine mds_r_dcmplx_6(name,value)
    character(*), intent(in) :: name
    complex(kind=dp), dimension(:,:,:,:,:,:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    status = mdsvalue(name//char(0),&
         &descr(13,value,dim1,dim2,dim3,dim4,dim5,dim6,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_dcmplx_6

  subroutine mds_r_dcmplx_7(name,value)
    character(*), intent(in) :: name
    complex(kind=dp), dimension(:,:,:,:,:,:,:), intent(out) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6,dim7
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    dim7 = size(value,7)
    status = mdsvalue(name//char(0),&
         &descr(13,value,dim1,dim2,dim3,dim4,dim5,dim6,dim7,0),0,len)
    call checkmds(status,"error reading "//name)
    return
  end subroutine mds_r_dcmplx_7

!--------------------------------------------
! mds_write subroutines
!--------------------------------------------

  subroutine mds_w_log_0(name,value)
    character(*), intent(in) :: name
    logical :: value
    logical :: status
    integer :: descr
    status = mdsput(name//char(0),"$",descr(8,value,0),0)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_log_0

  subroutine mds_w_char_1(name,value)
    character(*), intent(in) :: name
    character*(*), intent(in) :: value
    logical :: status
    integer :: descr
    status = mdsput(name//char(0),"$",&
         &descr(14,value,0,len(value)),0)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_char_1

  subroutine mds_w_int_0(name,value)
    character(*), intent(in) :: name
    integer , intent(in) :: value
    logical :: status
    integer :: descr
    status = mdsput(name//char(0),"$",descr(8 ,value, 0),0)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_int_0

  subroutine mds_w_int_1(name,value)
    character(*), intent(in) :: name
    integer , dimension(:), intent(in) :: value
    logical :: status
    integer :: descr
    integer :: dim1
    dim1 = size(value,1)
    status = mdsput(name//char(0),"$",descr(8 ,value, dim1, 0),0)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_int_1

  subroutine mds_w_int_2(name,value)
    character(*), intent(in) :: name 
    integer , dimension(:,:), intent(in) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2
    dim1 = size(value,1)
    dim2 = size(value,2)
    status = mdsput(name//char(0),&
         &"$",descr(8,value,dim1,dim2,0),0,len)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_int_2

  subroutine mds_w_int_3(name,value)
    character(*), intent(in) :: name 
    integer , dimension(:,:,:), intent(in) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    status = mdsput(name//char(0),&
         &"$",descr(8,value,dim1,dim2,dim3,0),0,len)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_int_3

  subroutine mds_w_int_4(name,value)
    character(*), intent(in) :: name 
    integer , dimension(:,:,:,:), intent(in) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    status = mdsput(name//char(0),&
         &"$",descr(8,value,dim1,dim2,dim3,dim4,0),0,len)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_int_4

  subroutine mds_w_int_5(name,value)
    character(*), intent(in) :: name 
    integer , dimension(:,:,:,:,:), intent(in) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    status = mdsput(name//char(0),&
         &"$",descr(8,value,dim1,dim2,dim3,dim4,dim5,0),0,len)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_int_5

  subroutine mds_w_int_6(name,value)
    character(*), intent(in) :: name 
    integer , dimension(:,:,:,:,:,:), intent(in) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    status = mdsput(name//char(0),&
         &"$",descr(8,value,dim1,dim2,dim3,dim4,dim5,dim6,0),0,len)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_int_6

  subroutine mds_w_int_7(name,value)
    character(*), intent(in) :: name 
    integer , dimension(:,:,:,:,:,:,:), intent(in) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6,dim7
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    dim7 = size(value,7)
    status = mdsput(name//char(0),&
         &"$",descr(8,value,dim1,dim2,dim3,dim4,dim5,dim6,dim7,0),0,len)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_int_7

  subroutine mds_w_real_0(name,value)
    character(*), intent(in) :: name
    real(kind=sp), intent(in) :: value
    logical :: status
    integer :: descr
    status = mdsput(name//char(0),"$",descr(10,value,0),0)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_real_0

  subroutine mds_w_real_1(name,value)
    character(*), intent(in) :: name
    real(kind=sp), dimension(:), intent(in) :: value
    logical :: status
    integer :: descr
    integer :: dim1
    dim1 = size(value,1)
    status = mdsput(name//char(0),"$",descr(10,value,dim1,0),0)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_real_1

  subroutine mds_w_real_2(name,value)
    character(*), intent(in) :: name
    real(kind=sp) , dimension(:,:), intent(in) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2
    dim1 = size(value,1)
    dim2 = size(value,2)
    status = mdsput(name//char(0),&
         &"$",descr(10,value,dim1,dim2,0),0,len)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_real_2

  subroutine mds_w_real_3(name,value)
    character(*), intent(in) :: name
    real(kind=sp) , dimension(:,:,:), intent(in) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    status = mdsput(name//char(0),&
         &"$",descr(10,value,dim1,dim2,dim3,0),0,len)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_real_3

  subroutine mds_w_real_4(name,value)
    character(*), intent(in) :: name
    real(kind=sp) , dimension(:,:,:,:), intent(in) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    status = mdsput(name//char(0),&
         &"$",descr(10,value,dim1,dim2,dim3,dim4,0),0,len)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_real_4

  subroutine mds_w_real_5(name,value)
    character(*), intent(in) :: name
    real(kind=sp) , dimension(:,:,:,:,:), intent(in) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    status = mdsput(name//char(0),&
         &"$",descr(10,value,dim1,dim2,dim3,dim4,dim5,0),0,len)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_real_5

  subroutine mds_w_real_6(name,value)
    character(*), intent(in) :: name
    real(kind=sp) , dimension(:,:,:,:,:,:), intent(in) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    status = mdsput(name//char(0),&
         &"$",descr(10,value,dim1,dim2,dim3,dim4,dim5,dim6,0),0,len)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_real_6

  subroutine mds_w_real_7(name,value)
    character(*), intent(in) :: name
    real(kind=sp) , dimension(:,:,:,:,:,:,:), intent(in) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6,dim7
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    dim7 = size(value,7)
    status = mdsput(name//char(0),&
         &"$",descr(10,value,dim1,dim2,dim3,dim4,dim5,dim6,dim7,0),0,len)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_real_7

  subroutine mds_w_doub_0(name,value)
    character(*), intent(in) :: name
    real(kind=dp), intent(in) :: value
    logical :: status
    integer :: descr
    status = mdsput(name//char(0),"$",descr(11,value,0),0)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_doub_0

  subroutine mds_w_doub_1(name,value)
    character(*), intent(in) :: name
    real(kind=dp), dimension(:), intent(in) :: value
    logical :: status
    integer :: descr
    integer :: dim1
    dim1 = size(value,1)
    status = mdsput(name//char(0),"$",descr(11,value,dim1,0),0)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_doub_1

  subroutine mds_w_doub_2(name,value)
    character(*), intent(in) :: name
    real(kind=dp), dimension(:,:), intent(in) :: value
    logical :: status
    integer :: descr
    integer :: dim1,dim2
    dim1 = size(value,1)
    dim2 = size(value,2)
    status = mdsput(name//char(0),"$",&
         &descr(11,value,dim1,dim2,0,0),0)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_doub_2

  subroutine mds_w_doub_3(name,value)
    character(*), intent(in) :: name
    real(kind=dp), dimension(:,:,:), intent(in) :: value
    logical :: status
    integer :: descr
    integer :: dim1,dim2,dim3
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    status = mdsput(name//char(0),"$",&
         &descr(11,value,dim1,dim2,dim3,0),0)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_doub_3

  subroutine mds_w_doub_4(name,value)
    character(*), intent(in) :: name
    real(kind=dp), dimension(:,:,:,:), intent(in) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    status = mdsput(name//char(0),&
         &"$",descr(11,value,dim1,dim2,dim3,dim4,0),0,len)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_doub_4

  subroutine mds_w_doub_5(name,value)
    character(*), intent(in) :: name
    real(kind=dp), dimension(:,:,:,:,:), intent(in) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    status = mdsput(name//char(0),&
         &"$",descr(11,value,dim1,dim2,dim3,dim4,dim5,0),0,len)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_doub_5

  subroutine mds_w_doub_6(name,value)
    character(*), intent(in) :: name
    real(kind=dp), dimension(:,:,:,:,:,:), intent(in) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    status = mdsput(name//char(0),&
         &"$",descr(11,value,dim1,dim2,dim3,dim4,dim5,dim6,0),0,len)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_doub_6

  subroutine mds_w_doub_7(name,value)
    character(*), intent(in) :: name
    real(kind=dp), dimension(:,:,:,:,:,:,:), intent(in) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6,dim7
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    dim7 = size(value,7)
    status = mdsput(name//char(0),&
         &"$",descr(11,value,dim1,dim2,dim3,dim4,dim5,dim6,dim7,0),0,len)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_doub_7

  subroutine mds_w_cmplx_0(name,value)
    character(*), intent(in) :: name
    complex(kind=sp), intent(in) :: value
    logical :: status
    integer :: descr
    status = mdsput(name//char(0),"$",descr(12,value,0),0)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_cmplx_0

  subroutine mds_w_cmplx_1(name,value)
    character(*), intent(in) :: name
    complex(kind=sp), dimension(:), intent(in) :: value
    logical :: status
    integer :: descr
    integer :: dim1
    dim1 = size(value,1)
    status = mdsput(name//char(0),"$",descr(12,value,dim1,0),0)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_cmplx_1

  subroutine mds_w_cmplx_2(name,value)
    character(*), intent(in) :: name
    complex(kind=sp), dimension(:,:), intent(in) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2
    dim1 = size(value,1)
    dim2 = size(value,2)
    status = mdsput(name//char(0),&
         &"$",descr(12,value,dim1,dim2,0),0,len)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_cmplx_2

  subroutine mds_w_cmplx_3(name,value)
    character(*), intent(in) :: name
    complex(kind=sp), dimension(:,:,:), intent(in) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    status = mdsput(name//char(0),&
         &"$",descr(12,value,dim1,dim2,dim3,0),0,len)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_cmplx_3

  subroutine mds_w_cmplx_4(name,value)
    character(*), intent(in) :: name
    complex(kind=sp), dimension(:,:,:,:), intent(in) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    status = mdsput(name//char(0),&
         &"$",descr(12,value,dim1,dim2,dim3,dim4,0),0,len)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_cmplx_4

  subroutine mds_w_cmplx_5(name,value)
    character(*), intent(in) :: name
    complex(kind=sp), dimension(:,:,:,:,:), intent(in) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    status = mdsput(name//char(0),&
         &"$",descr(12,value,dim1,dim2,dim3,dim4,dim5,0),0,len)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_cmplx_5

  subroutine mds_w_cmplx_6(name,value)
    character(*), intent(in) :: name
    complex(kind=sp), dimension(:,:,:,:,:,:), intent(in) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    status = mdsput(name//char(0),&
         &"$",descr(12,value,dim1,dim2,dim3,dim4,dim5,dim6,0),0,len)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_cmplx_6

  subroutine mds_w_cmplx_7(name,value)
    character(*), intent(in) :: name
    complex(kind=sp), dimension(:,:,:,:,:,:,:), intent(in) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6,dim7
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    dim7 = size(value,7)
    status = mdsput(name//char(0),&
         &"$",descr(12,value,dim1,dim2,dim3,dim4,dim5,dim6,dim7,0),0,len)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_cmplx_7

  subroutine mds_w_dcmplx_0(name,value)
    character(*), intent(in) :: name
    complex(kind=dp), intent(in) :: value
    logical :: status
    integer :: descr
    status = mdsput(name//char(0),"$",descr(13,value,0),0)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_dcmplx_0

  subroutine mds_w_dcmplx_1(name,value)
    character(*), intent(in) :: name
    complex(kind=dp), dimension(:), intent(in) :: value
    logical :: status
    integer :: descr
    integer :: dim1
    dim1 = size(value,1)
    status = mdsput(name//char(0),"$",descr(13,value,dim1,0),0)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_dcmplx_1

  subroutine mds_w_dcmplx_2(name,value)
    character(*), intent(in) :: name
    complex(kind=dp), dimension(:,:), intent(in) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2
    dim1 = size(value,1)
    dim2 = size(value,2)
    status = mdsput(name//char(0),&
         &"$",descr(13,value,dim1,dim2,0),0,len)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_dcmplx_2

  subroutine mds_w_dcmplx_3(name,value)
    character(*), intent(in) :: name
    complex(kind=dp), dimension(:,:,:), intent(in) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    status = mdsput(name//char(0),&
         &"$",descr(13,value,dim1,dim2,dim3,0),0,len)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_dcmplx_3

  subroutine mds_w_dcmplx_4(name,value)
    character(*), intent(in) :: name
    complex(kind=dp), dimension(:,:,:,:), intent(in) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    status = mdsput(name//char(0),&
         &"$",descr(13,value,dim1,dim2,dim3,dim4,0),0,len)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_dcmplx_4

  subroutine mds_w_dcmplx_5(name,value)
    character(*), intent(in) :: name
    complex(kind=dp), dimension(:,:,:,:,:), intent(in) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    status = mdsput(name//char(0),&
         &"$",descr(13,value,dim1,dim2,dim3,dim4,dim5,0),0,len)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_dcmplx_5

  subroutine mds_w_dcmplx_6(name,value)
    character(*), intent(in) :: name
    complex(kind=dp), dimension(:,:,:,:,:,:), intent(in) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    status = mdsput(name//char(0),&
         &"$",descr(13,value,dim1,dim2,dim3,dim4,dim5,dim6,0),0,len)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_dcmplx_6

  subroutine mds_w_dcmplx_7(name,value)
    character(*), intent(in) :: name
    complex(kind=dp), dimension(:,:,:,:,:,:,:), intent(in) :: value
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6,dim7
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    dim7 = size(value,7)
    status = mdsput(name//char(0),&
         &"$",descr(13,value,dim1,dim2,dim3,dim4,dim5,dim6,dim7,0),0,len)
    call checkmds(status,"error writing "//name)
    return
  end subroutine mds_w_dcmplx_7

end module mdsio


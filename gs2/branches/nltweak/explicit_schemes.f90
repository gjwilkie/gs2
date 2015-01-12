module explicit_schemes
  implicit none

  private

  public :: multistep_scheme

  type multistep_scheme
     private
     integer :: MaxSupportedOrder
     real, dimension(:,:), allocatable :: Coeffs
     character(len=128) :: name
   contains
     procedure :: init => multistep_scheme_init
     procedure :: set_coeffs => multistep_scheme_set_coeffs
     procedure :: get_coeff => multistep_scheme_get_coeff
     procedure :: finish => multistep_scheme_finish
     procedure :: print => multistep_print
  end type multistep_scheme

contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! TYPE BOUND PROCEDURES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!###########################
!# MULTISTEP SCHEMES
!###########################

  subroutine multistep_scheme_init(self,MaxOrder,name)
    implicit none
    class(multistep_scheme), intent(in out) :: self
    integer, intent(in) :: MaxOrder
    character(len=*), intent(in) :: name

    self%MaxSupportedOrder=MaxOrder
    self%name=name
    allocate(self%Coeffs(MaxOrder,MaxOrder))
    self%Coeffs=0
  end subroutine multistep_scheme_init

  subroutine multistep_scheme_finish(self)
    implicit none
    class(multistep_scheme), intent(in out) :: self
    if(allocated(self%Coeffs)) deallocate(self%Coeffs)
    self%MaxSupportedOrder=0
    self%name=''
  end subroutine multistep_scheme_finish

  subroutine multistep_scheme_set_coeffs(self,Order,CoeffIn)
    implicit none
    class(multistep_scheme), intent(in out) :: self
    integer, intent(in) :: Order !Order to set
    real, dimension(:), intent(in) :: CoeffIn !Coefficients for current Order
    !Could check here that 0<Order<=MaxOrder
    self%Coeffs(1:Order,Order)=CoeffIn !Expect to have as many coefficients as Order
  end subroutine multistep_scheme_set_coeffs

  function multistep_scheme_get_coeff(self,Order,Position)
    implicit none
    class(multistep_scheme), intent(in) :: self
    integer, intent(in) :: Order !Order to set
    integer, intent(in) :: Position !Which position do we want?
    real :: multistep_scheme_get_coeff
    !Could check here that 0<Order<=MaxOrder
    !Could check here that 0<Position<=Order

    multistep_scheme_get_coeff=self%Coeffs(Position,Order)
  end function multistep_scheme_get_coeff

  subroutine multistep_print(self)
    implicit none
    class(multistep_scheme), intent(in) :: self
    integer :: i, j

    write(6,'("Multistep scheme : ",A)') self%name
    write(6,'("Max order        : ",I0)') self%MaxSupportedOrder
    write(6,'()')
    do i=1,self%MaxSupportedOrder
       write(6,'("Coefficients for order : ",I0)') i
       do j=1,i
          write(6,'(F12.5,",")',advance='no') self%get_coeff(i,j)
       enddo
       write(6,'()')
    enddo
  end subroutine multistep_print
end module explicit_schemes

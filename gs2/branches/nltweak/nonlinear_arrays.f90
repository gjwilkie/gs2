module nonlinear_arrays
  use explicit_schemes, only: multistep_scheme

  implicit none
  private

  public :: get_exp_source, init_multistep, finish_multistep
  public :: nl_order, current_order, gexp, nonlin, store_history

  !Internal step counter, used to determine what order scheme to use
  integer :: current_order

  integer :: nl_order !Order of method to use

  !The source term history
  complex, dimension (:,:,:,:), allocatable :: gexp
  ! (explicit_order,-ntgrid:ntgrid,2, -g-layout-)

  logical :: nonlin = .false.
  logical :: store_history

  !!Advance scheme
  type(multistep_scheme) :: multistep !Currently only Adams-Bashforth supported
  integer, parameter :: max_adams_order=6

contains
  subroutine init_multistep
    implicit none

    !Setup the multistep scheme
    !Could have case select here to allow support for different schemes
    call multistep%init(max_adams_order,'Adams-Bashforth')
    !Note it is possible to calculate these coefficients for any arbritary order
    !if we ever decide we want to support this.
    call multistep%set_coeffs(1,[1]/1.0)
    call multistep%set_coeffs(2,[3,-1]/2.0)
    call multistep%set_coeffs(3,[23,-16,5]/12.0)
    call multistep%set_coeffs(4,[55,-59,37,-9]/24.0)
    call multistep%set_coeffs(5,[1901,-2774,2616,-1274,251]/720.0)
    call multistep%set_coeffs(6,[4277,-7923,9982,-7298,2877,-475]/1440.0)
    !if(proc0)call multistep%print
  end subroutine init_multistep

  subroutine finish_multistep
    implicit none
    call multistep%finish
  end subroutine finish_multistep

  !>Calculate the explicit source term at given location
  !!Note, doesn't include 0.5*dt needed to actually advance
  function get_exp_source(ig,isgn,iglo,orderIn)
    implicit none
    complex :: get_exp_source
    integer, intent(in) :: ig, isgn, iglo
    integer, intent(in), optional ::orderIn
    integer :: desiredOrder,order, i

    !Initialise return value
    get_exp_source=0.0

    !Exit early if no steps taken
    if(current_order.eq.0) return

    !Set the target order
    desiredOrder=nl_order
    if(present(orderIn)) desiredOrder=orderIn

    !Set the order to use 
    order=min(desiredOrder,current_order)

    !Build up term
    do i=1,order
       get_exp_source=get_exp_source+multistep%get_coeff(order,i)*gexp(i,ig,isgn,iglo)
    enddo    
  end function get_exp_source
end module nonlinear_arrays

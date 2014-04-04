
!> This module contains functions for averaging quantities
!! which have the dimensions of space, and potentially species.
!! For example, (y), (x), (x,y), (theta,x,y), (x,y,s) etc.
!! The names of the interfaces are self-explanatory: for example,
!! the function average_theta averages the theta dimension.
!! 
!! The quantities may or may not be distributed across processes.
!! If they are not local, they must be distributed in the same way
!! that the fields are distributed, i.e. in x and y only, and
!! according to the same pattern.
module volume_averages

  use theta_grid, only: ntheta, nperiod
  use fields_parallelization, only: field_k_local
  use species, only: nspec

  implicit none
  integer :: ntg_out
  interface average_all
    module procedure average_all_txy
    module procedure average_all_real_xys
  end interface average_all

  interface average_theta
    module procedure average_theta_txy
  end interface

  interface average_ky
    module procedure average_ky_xy
  end interface

  interface average_kx
    module procedure average_kx_xy
    module procedure average_kx_xys
  end interface


contains
  subroutine init_volume_averages
    use theta_grid, only: ntgrid
    !ntg_out =  ntheta/2 + (nperiod-1)*ntheta
    ntg_out = ntgrid
  end subroutine init_volume_averages

  !subroutine average_theta_distributed(a, b, axb_by_mode)
    !use theta_grid, only: ntgrid, delthet, jacob, nperiod, ntheta
    !use kt_grids, only: naky, ntheta0
    !implicit none
    !complex, dimension (-ntgrid:,:,:), intent (in) :: a, b
    !real, dimension (:,:), intent (out) :: axb_by_mode


    !integer :: ik, it
    !integer :: ng
    !real, dimension (-ntg_out:ntg_out) :: wgt
    !real :: anorm

    !ng = ntg_out
    !wgt = delthet(-ng:ng)*jacob(-ng:ng)
    !anorm = sum(wgt)

    !do ik = 1, naky
       !do it = 1, ntheta0
         !if (field_k_local(it,ik)) then
            !axb_by_mode(it,ik) &
                 != sum(real(conjg(a(-ng:ng,it,ik))*b(-ng:ng,it,ik))*wgt)/anorm
          !end if
       !end do
    !end do

  !end subroutine average_theta_distributed
  subroutine average_all_real_xys (a, avg, distributed)
    use theta_grid, only: ntgrid, delthet, jacob, nperiod, ntheta
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
    implicit none
    real, dimension (ntheta0, naky,nspec), intent (in) :: a
    real, dimension (ntheta0, nspec) :: axb_by_kx
    real, dimension(nspec), intent(out) :: avg
    logical,intent(in) :: distributed
    integer :: is

    do is = 1,nspec
      call average_ky(a(:,:,is), axb_by_kx(:,is), distributed)
      avg(is) = sum(axb_by_kx(:,is))
    end do
  end subroutine average_all_real_xys

  subroutine average_all_txy (a, b, avg, distributed)
    use theta_grid, only: ntgrid, delthet, jacob, nperiod, ntheta
    use kt_grids, only: naky, ntheta0
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: a, b
    real, intent(out) :: avg
    real, dimension (ntheta0,naky) :: axb_by_mode
    real, dimension (naky) :: axb_by_ky
    logical,intent(in) :: distributed

    call average_theta(a, b, axb_by_mode, distributed)
    call average_kx(axb_by_mode, axb_by_ky, distributed)
    avg = sum(axb_by_ky)
  end subroutine average_all_txy

  subroutine average_theta_txy (a, b, axb_by_mode, distributed)
    use theta_grid, only: ntgrid, delthet, jacob, nperiod, ntheta
    use kt_grids, only: naky, ntheta0
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: a, b
    real, dimension (:,:), intent (out) :: axb_by_mode
    logical,intent(in) :: distributed


    integer :: ik, it
    integer :: ng
    real, dimension (-ntg_out:ntg_out) :: wgt
    real :: anorm

    ng = ntg_out
    wgt = delthet(-ng:ng)*jacob(-ng:ng)
    anorm = sum(wgt)

    axb_by_mode = 0

    do ik = 1, naky
       do it = 1, ntheta0
         if (.not. distributed .or. field_k_local(it,ik)) then
          axb_by_mode(it,ik) &
               = sum(real(conjg(a(-ng:ng,it,ik))*b(-ng:ng,it,ik))*wgt)/anorm
         end if
       end do
    end do

  end subroutine average_theta_txy

  subroutine average_ky_xy (f, favg, distributed)
    use kt_grids, only: naky, ntheta0, aky
    use fields_parallelization, only: field_k_local
    use mp, only: sum_allreduce
    implicit none
    real, dimension (:,:), intent (in) :: f
    real, dimension (:), intent (out) :: favg
    logical,intent(in) :: distributed
    real :: fac
    integer :: ik, it

! ky=0 modes have correct amplitudes; rest must be scaled
! note contrast with scaling factors in FFT routines.

!CMR+GC: 2/9/2013
!  fac values here arise because gs2's Fourier coefficients, F_k^gs2, not standard form: 
!          i.e. f(x) = f_k e^(i k.x)
!  With standard Fourier coeffs in gs2, we would instead need:  fac=2.0 for ky > 0
!      (fac=2.0 would account ky<0 contributions, not stored due to reality condition)

    favg = 0.
    do ik = 1, naky
       fac = 0.5
       if (aky(ik) == 0.) fac = 1.0
       do it = 1, ntheta0
         if (.not. distributed .or. field_k_local(it,ik)) then

          favg(it) = favg(it) + f(it, ik) * fac
         end if
       end do
    end do
    if (distributed) call sum_allreduce(favg)

  end subroutine average_ky_xy
  subroutine average_kx_xys (f, favg, distributed)
    use kt_grids, only: naky, ntheta0, aky
    use fields_parallelization, only: field_k_local
    use mp, only: sum_allreduce
    implicit none
    real, dimension (:,:,:), intent (in) :: f
    real, dimension (:,:), intent (out) :: favg
    logical,intent(in) :: distributed
    integer :: is
    do is = 1,nspec
      call average_kx_xy(f(:,:,is), favg(:, is), distributed)
    end do
  end subroutine average_kx_xys
  subroutine average_kx_xy (f, favg, distributed)
    use kt_grids, only: naky, ntheta0, aky
    use fields_parallelization, only: field_k_local
    use mp, only: sum_allreduce
    implicit none
    real, dimension (:,:), intent (in) :: f
    real, dimension (:), intent (out) :: favg
    logical,intent(in) :: distributed
    real :: fac
    integer :: ik, it

! ky=0 modes have correct amplitudes; rest must be scaled
! note contrast with scaling factors in FFT routines.

!CMR+GC: 2/9/2013
!  fac values here arise because gs2's Fourier coefficients, F_k^gs2, not standard form: 
!          i.e. f(x) = f_k e^(i k.x)
!  With standard Fourier coeffs in gs2, we would instead need:  fac=2.0 for ky > 0
!      (fac=2.0 would account ky<0 contributions, not stored due to reality condition)

    favg = 0.
    do ik = 1, naky
       fac = 0.5
       if (aky(ik) == 0.) fac = 1.0
       do it = 1, ntheta0
         if (.not. distributed .or. field_k_local(it,ik)) then

          favg(ik) = favg(ik) + f(it, ik) * fac
         end if
       end do
    end do
    if (distributed) call sum_allreduce(favg)

  end subroutine average_kx_xy

end module volume_averages

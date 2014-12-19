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

  implicit none

  private

  public :: average_all, average_theta, average_ky, average_kx
  public :: init_volume_averages

  integer :: ntg_out

  interface average_all
    module procedure average_all_txy
    module procedure average_all_real_xys
    module procedure average_all_real_xy
    module procedure average_all_complex_xy
    module procedure average_all_2_complex_xy
    module procedure average_all_2_complex_txy
  end interface average_all

  interface average_theta
    module procedure average_theta_txy
    module procedure average_theta_1_real_txy
    module procedure average_theta_2_complex_txy
    module procedure average_theta_txys
    module procedure average_theta_complex_complex_txys

    module procedure average_theta_real_real_t
    module procedure average_theta_complex_complex_t
    module procedure average_theta_complex_complex_real_t
    module procedure average_theta_complex_complex_complex_t

  end interface

  interface average_ky
    module procedure average_ky_xy
    module procedure average_ky_complex_xy
    module procedure average_ky_2_complex_xy
    module procedure average_ky_xys
  end interface

  interface average_kx
    module procedure average_kx_xy
    module procedure average_kx_xys
  end interface

contains
  subroutine init_volume_averages
    use theta_grid, only: ntgrid
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

  subroutine average_all_2_complex_txy(a, b, avg, distributed)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    implicit none
    complex, dimension (-ntgrid:ntgrid, ntheta0, naky), intent (in) :: a, b
    complex, dimension (ntheta0, naky) :: axb_by_mode
    complex, dimension (ntheta0) :: axb_by_kx
    complex, intent(out) :: avg
    logical,intent(in) :: distributed

    call average_theta(a, b, axb_by_mode, distributed)
    call average_ky(axb_by_mode, axb_by_kx, distributed)
    avg = sum(axb_by_kx(:))
  end subroutine average_all_2_complex_txy

  subroutine average_all_2_complex_xy(a, b, avg, distributed)
    use kt_grids, only: naky, ntheta0
    implicit none
    complex, dimension (ntheta0, naky), intent (in) :: a, b
    complex, dimension (ntheta0) :: axb_by_kx
    complex, intent(out) :: avg
    logical,intent(in) :: distributed

    call average_ky(a(:,:), b(:,:), axb_by_kx(:), distributed)
    avg = sum(axb_by_kx(:))
  end subroutine average_all_2_complex_xy

  subroutine average_all_complex_xy(a, avg, distributed)
    use kt_grids, only: naky, ntheta0
    implicit none
    complex, dimension (ntheta0, naky), intent (in) :: a
    complex, dimension (ntheta0) :: axb_by_kx
    complex, intent(out) :: avg
    logical,intent(in) :: distributed

    call average_ky(a(:,:), axb_by_kx(:), distributed)
    avg = sum(axb_by_kx(:))
  end subroutine average_all_complex_xy

  subroutine average_all_real_xy(a, avg, distributed)
    use kt_grids, only: naky, ntheta0
    implicit none
    real, dimension (ntheta0, naky), intent (in) :: a
    real, dimension (ntheta0) :: axb_by_kx
    real, intent(out) :: avg
    logical,intent(in) :: distributed

    call average_ky(a(:,:), axb_by_kx(:), distributed)
    avg = sum(axb_by_kx(:))
  end subroutine average_all_real_xy

  subroutine average_all_real_xys (a, avg, distributed)
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
    implicit none
    real, dimension (ntheta0, naky, nspec), intent (in) :: a
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
    use theta_grid, only: ntgrid
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

  subroutine average_theta_complex_complex_txys (a, axb_by_mode, distributed)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use fields_parallelization, only: field_k_local
    use species, only: nspec
    implicit none
    complex, dimension (-ntgrid:,:,:,:), intent (in) :: a
    complex, dimension (:,:,:), intent (out) :: axb_by_mode
    logical,intent(in) :: distributed
    integer :: ik, it, is
    integer :: ng

    ng = ntg_out
    axb_by_mode = 0

    do is = 1,nspec
       do ik = 1, naky
          do it = 1, ntheta0
             if (.not. distributed .or. field_k_local(it,ik)) then
                call average_theta(a(-ng:ng,it,ik,is), axb_by_mode(it,ik,is))
             end if
          end do
       end do
    end do
  end subroutine average_theta_complex_complex_txys

  subroutine average_theta_txys (a, b, axb_by_mode, distributed)
    use theta_grid, only: ntgrid, delthet, jacob
    use kt_grids, only: naky, ntheta0
    use fields_parallelization, only: field_k_local
    use species, only: nspec
    implicit none
    complex, dimension (-ntgrid:,:,:,:), intent (in) :: a, b
    real, dimension (:,:,:), intent (out) :: axb_by_mode
    logical,intent(in) :: distributed
    integer :: ik, it, is
    integer :: ng
    real, dimension (-ntg_out:ntg_out) :: wgt
    real :: anorm

    ng = ntg_out
    wgt = delthet(-ng:ng)*jacob(-ng:ng)
    anorm = sum(wgt)
    axb_by_mode = 0

    do is = 1,nspec
       do ik = 1, naky
          do it = 1, ntheta0
             if (.not. distributed .or. field_k_local(it,ik)) then
                axb_by_mode(it,ik,is) &
                     = sum(real(conjg(a(-ng:ng,it,ik,is))*b(-ng:ng,it,ik,is))*wgt)/anorm
             end if
          end do
       end do
    end do
  end subroutine average_theta_txys

  subroutine average_theta_1_real_txy (a, axb_by_mode, distributed)
    use theta_grid, only: ntgrid, delthet, jacob
    use kt_grids, only: naky, ntheta0
    use fields_parallelization, only: field_k_local
    implicit none
    real, dimension (-ntgrid:,:,:), intent (in) :: a
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
                  = sum((a(-ng:ng,it,ik))*wgt)/anorm
          end if
       end do
    end do
  end subroutine average_theta_1_real_txy

  subroutine average_theta_real_real_t (a, axb)
    use theta_grid, only: ntgrid, delthet, jacob
    implicit none
    real, dimension (-ntgrid:), intent (in) :: a
    real, intent (out) :: axb
    integer :: ng
    real, dimension (-ntg_out:ntg_out) :: wgt
    real :: anorm

    ng = ntg_out
    wgt = delthet(-ng:ng)*jacob(-ng:ng)
    anorm = sum(wgt)

    axb = sum(a*wgt)/anorm
  end subroutine average_theta_real_real_t

  subroutine average_theta_complex_complex_complex_t (a, b, axb)
    use theta_grid, only: ntgrid, delthet, jacob
    implicit none
    complex, dimension (-ntgrid:), intent (in) :: a, b
    complex, intent (out) :: axb
    integer :: ng
    real, dimension (-ntg_out:ntg_out) :: wgt
    real :: anorm

    ng = ntg_out
    wgt = delthet(-ng:ng)*jacob(-ng:ng)
    anorm = sum(wgt)

    axb = sum(conjg(a(-ng:ng))*b(-ng:ng)*wgt)/anorm
  end subroutine average_theta_complex_complex_complex_t

  subroutine average_theta_complex_complex_t (a, axb)
    use theta_grid, only: ntgrid, delthet, jacob
    implicit none
    complex, dimension (-ntgrid:), intent (in) :: a
    complex, intent (out) :: axb
    integer :: ng
    real, dimension (-ntg_out:ntg_out) :: wgt
    real :: anorm

    ng = ntg_out
    wgt = delthet(-ng:ng)*jacob(-ng:ng)
    anorm = sum(wgt)

    axb = sum(a(-ng:ng)*wgt)/anorm
  end subroutine average_theta_complex_complex_t

  subroutine average_theta_complex_complex_real_t (a, b, axb)
    use theta_grid, only: ntgrid
    implicit none
    complex, dimension (-ntgrid:), intent (in) :: a, b
    real, intent (out) :: axb
    complex :: axbcomplex

    call average_theta(a, b, axbcomplex)
    axb = real(axbcomplex)
  end subroutine average_theta_complex_complex_real_t

  subroutine average_theta_txy (a, b, axb_by_mode, distributed)
    use theta_grid, only: ntgrid, delthet, jacob
    use kt_grids, only: naky, ntheta0
    use fields_parallelization, only: field_k_local
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

  subroutine average_theta_2_complex_txy (a, b, axb_by_mode, distributed)
    use theta_grid, only: ntgrid, delthet, jacob
    use kt_grids, only: naky, ntheta0
    use fields_parallelization, only: field_k_local
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: a, b
    complex, dimension (:,:), intent (out) :: axb_by_mode
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
                  = sum(conjg(a(-ng:ng,it,ik))*b(-ng:ng,it,ik)*wgt)/anorm
          end if
       end do
    end do
  end subroutine average_theta_2_complex_txy

  subroutine average_ky_complex_xy (f, favg, distributed)
    use kt_grids, only: naky, ntheta0, aky
    use fields_parallelization, only: field_k_local
    use mp, only: sum_allreduce
    implicit none
    complex, dimension (:,:), intent (in) :: f
    complex, dimension (:), intent (out) :: favg
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
  end subroutine average_ky_complex_xy

  subroutine average_ky_2_complex_xy (f, g, favg, distributed)
    use kt_grids, only: naky, ntheta0, aky
    use fields_parallelization, only: field_k_local
    use mp, only: sum_allreduce
    implicit none
    complex, dimension (:,:), intent (in) :: f, g
    complex, dimension (:), intent (out) :: favg
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
             favg(it) = favg(it) + conjg(f(it, ik)) * g(it,ik) * fac
          end if
       end do
    end do
    if (distributed) call sum_allreduce(favg)
  end subroutine average_ky_2_complex_xy

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

  subroutine average_ky_xys (f, favg, distributed)
    use kt_grids, only: naky, ntheta0, aky
    use species, only: nspec
    use fields_parallelization, only: field_k_local
    use mp, only: sum_allreduce
    implicit none
    real, dimension (:,:,:), intent (in) :: f
    real, dimension (:,:), intent (out) :: favg
    logical,intent(in) :: distributed
    real :: fac
    integer :: ik, it, is

! ky=0 modes have correct amplitudes; rest must be scaled
! note contrast with scaling factors in FFT routines.

!CMR+GC: 2/9/2013
!  fac values here arise because gs2's Fourier coefficients, F_k^gs2, not standard form: 
!          i.e. f(x) = f_k e^(i k.x)
!  With standard Fourier coeffs in gs2, we would instead need:  fac=2.0 for ky > 0
!      (fac=2.0 would account ky<0 contributions, not stored due to reality condition)

    favg = 0.
    do is = 1, nspec
       do ik = 1, naky
          fac = 0.5
          if (aky(ik) == 0.) fac = 1.0
          do it = 1, ntheta0
             if (.not. distributed .or. field_k_local(it,ik)) then
                favg(it,is) = favg(it,is) + f(it, ik,is) * fac
             end if
          end do
       end do
    end do
    if (distributed) call sum_allreduce(favg)
  end subroutine average_ky_xys

  subroutine average_kx_xys (f, favg, distributed)
    use species, only: nspec
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
    use mp, only: sum_allreduce
    use fields_parallelization, only: field_k_local
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

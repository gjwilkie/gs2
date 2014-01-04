
module volume_averages

  use theta_grid, only: ntheta, nperiod
  use fields_parallelization, only: field_k_local
    integer :: ntg_out


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

  subroutine average_theta (a, b, axb_by_mode, distributed)
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

  end subroutine average_theta

  subroutine average_ky (f, favg, distributed)
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

  end subroutine average_ky
  subroutine average_kx (f, favg, distributed)
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

  end subroutine average_kx
end module volume_averages

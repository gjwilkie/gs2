module gs2_time

  implicit none

  private

  real :: user_dt, code_dt 
  real :: user_tstart, code_tstart
  real :: user_dt_cfl = -1.
  real :: code_dt_cfl = -1.

  real :: code_dt_max, user_dt_max
  real :: code_dt_min, user_dt_min

  ! added May 18, 2009 to take care of problems
  ! in exb_shear calculation after change in time step size
  real :: code_dt_old = 0.

! GGH points out that this initialization is not necessary (we think)
  real :: user_time = 0.
  real :: code_time = 0.

!  real :: dt
  real :: vtfac

  public :: user_dt, code_dt, update_time, code_dt_old
  public :: user_time, code_time, code_dt_max
  public :: save_dt_min, save_dt, save_dt_cfl, write_dt
  public :: init_tstart, init_delt
  public :: code_dt_cfl, code_dt_min, user2code
  
contains

  subroutine init_tstart (tstart)

    real, intent (in) :: tstart

    user_time = tstart
    code_time = tstart * vtfac

  end subroutine init_tstart

  subroutine init_delt (delt, tnorm)
    real, intent (in) :: delt, tnorm

! The decision of how to define v_t, i.e., with or without
! a sqrt(2) factor, is made in kt_grids, with the "norm_option"
! variable.  It would be TOO CRAZY to allow v_t to be defined 
! in two different modules, and I do not want to introduce a 
! module solely to define v_t, so there is a little dance here.
!
! The point is:
!
! If norm_option -> v_t == sqrt(2 T/m ) then tnorm = 1.
! If norm_option -> v_t == sqrt(T/m ) then tnorm = sqrt(2)
!
! This comes up because time is normalized as follows:
!
! t_phys = t_N * a / v_t
!
! Internal to the code, v_t is always v_t == sqrt(2T/m)
! because that is how MTK coded everything.
!
! delt_in is a user input, from the run_parameters module.
! In a perfect world, we could have a gs2_time namelist. 
!
! So: we require the time in the code to be normalized by 
!
! t_N = t_phys * sqrt(2T/m) / a
! 
! When v_t is defined to be sqrt(T/m), the additional factor of 
! sqrt(2) is supplied by tnorm.
!
! 
    vtfac = tnorm
    
    user_dt = delt
    code_dt = delt * vtfac

  end subroutine init_delt

  subroutine update_time
! MAB+CMR, 21/5/09: set code_dt_old to code_dt BEFORE any changes in timestep
    code_dt_old = code_dt
    code_time = code_time + code_dt
    user_time = user_time + user_dt

  end subroutine update_time

  subroutine save_dt_cfl (delt_cfl)

    real, intent (in) :: delt_cfl

    code_dt_cfl = delt_cfl
    user_dt_cfl = delt_cfl / vtfac

  end subroutine save_dt_cfl

  subroutine save_dt_min (dt_min)

    real, intent (in) :: dt_min

    user_dt_min = dt_min
    code_dt_min = dt_min * vtfac

  end subroutine save_dt_min

  subroutine save_dt(delt)
    
    real, intent (in) :: delt

    code_dt = delt
    user_dt = delt / vtfac 
    
  end subroutine save_dt

  subroutine write_dt

    if (user_dt_cfl > 0. .and. user_dt_cfl < 1.e7) &
         write(*,*) user_dt_cfl,' : ',user_dt
       
  end subroutine write_dt

  subroutine user2code (usertime, codetime)

    real, intent (in) :: usertime
    real, intent (out) :: codetime

    codetime = usertime * vtfac

  end subroutine user2code

end module gs2_time

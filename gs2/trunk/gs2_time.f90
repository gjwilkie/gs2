module gs2_time

  implicit none

  real :: simtime = 0.
  real, private :: dt
  real, private :: dt_cfl = -1.

contains

  subroutine init_time (tstart, delt)
    real, intent (in) :: tstart, delt

    simtime = tstart
    dt = delt

  end subroutine init_time

  subroutine update (delt)
    real, intent (in) :: delt

    simtime = simtime + delt

  end subroutine update

  function stime() 
    real :: stime
    
    stime = simtime

  end function stime

  function simdt() 
    real :: simdt
    
    simdt = dt

  end function simdt

  subroutine save_dt(delt, delt_cfl)

    real :: delt, delt_cfl

    dt = delt
    dt_cfl = delt_cfl

  end subroutine save_dt

  subroutine write_dt(tnorm)

    real, intent (in) :: tnorm

    if (dt_cfl > 0.) then
       write(*,*) dt_cfl/tnorm,' : ',dt/tnorm
    end if
       
  end subroutine write_dt

end module gs2_time

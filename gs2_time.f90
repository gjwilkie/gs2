module gs2_time

  implicit none

  real :: simtime = 0.

contains

  subroutine init (tstart)
    real, intent (in) :: tstart

    simtime = tstart

  end subroutine init

  subroutine update (delt)
    real, intent (in) :: delt

    simtime = simtime + delt

  end subroutine update

  function stime() 
    real :: stime
    
    stime = simtime

  end function stime

end module gs2_time

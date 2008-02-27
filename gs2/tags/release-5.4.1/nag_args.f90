  function iargc()
    use f90_unix_env, only: nag_iargc => iargc
    integer :: iargc
    iargc = nag_iargc ()
  end function iargc

  subroutine getarg (k, arg)
    use f90_unix_env, only: nag_getarg => getarg
    integer,       intent (in)  :: k
    character (*), intent (out) :: arg

    call nag_getarg (k, arg)
  end subroutine getarg


# include "define.inc"

module ran
! D. Ernst: change Cray/SGI ranf to standard f90 calls to random_number
  implicit none
  private

# if FCOMPILER == _CRAY_
! make at least one public item to avoid a compiler warning  
  public :: ran_dummy_var
  logical :: ran_dummy_var
# else

  public :: ranf

contains

# ifdef USE_NR_RAN
  function ranf(idum_arg)
!CMR 6/4/2008: fix ranf by NOT resetting idum to 1 on every call
!            : and to allow resetting the seed
!

!.....Replaces cray routine ranf(). What follows is the
!     same as FUNCTION RAN1(idum) of numerical recipes p196,
!     except idum is set explicitly below.

!.....RAN1(idum):
!.....Returns a uniform deviate between 0.0 and 1.0.
!     Set idum to any negative value to initialize
!     or reinitialize the sequence

    integer, optional, intent(inout) :: idum_arg
    integer :: idum
    real :: ranf
    integer, parameter :: ia=16807,im=2147483647,iq=127773,ir=2836
    integer, parameter :: ntab=32,ndiv=1+(im-1)/ntab
    real, parameter :: am=1./im,eps=1.2e-7,rnmx=1.-eps
!    parameter (ia=16807,im=2147483647,am=1./im,iq=127773,ir=2836, &
!         ntab=32,ndiv=1+(im-1)/ntab,eps=1.2e-7,rnmx=1.-eps)
    integer :: j,k,iv(ntab),iy
    save iv,iy,idum
    data iv /ntab*0/, iy /0/

!.....To conform to cray, do not pass idum
!CMR, 6/4/2008:  do not reset idum on every call
!CMR             and allow for setting random number seed
    data idum /1/   !CMR: instead of  idum=1
    if(present(idum_arg)) idum=idum_arg
    
    if (idum <= 0.or.iy == 0) then
       idum=max(-idum,1)
       do j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum < 0) idum=idum+IM
          if (j <= NTAB) iv(j)=idum
       enddo
       iy=iv(1)
    endif
    k=idum/IQ
    idum=IA*(idum-k*IQ)-IR*k
    if (idum < 0) idum=idum+IM
    j=1+iy/NDIV
    iy=iv(j)
    iv(j)=idum
!     ran1=min(AM*iy,RNMX) (change name to ranf)
    ranf=min(AM*iy,RNMX)
    if(present(idum_arg)) idum_arg=idum
    
  end function ranf
!  (C) Copr. 1986-92 Numerical Recipes Software.
# else
  function ranf()
    real :: ranf, random
    call random_number(random)
    ranf = random
  end function ranf
# endif
# endif

end module ran

! This file is used when the NAG library is not available.

!program random_test
!  implicit none
!  real g05cae, xdum
!  integer i, idum
!
!  idum=1
!  call g05cbe(idum)
!  print *,    'first number should be 0.38537207475475027
!  do i=1,10
!     print *, 'nag emulation        =', g05cae(xdum)
!  end do
!end program random_test

! this is intended as a private module for use only within this file.
! callable routines follow this module:
module g05cae_mod
  implicit none
contains

! Todo:
! The following random number generator works correctly, and exactly
! reproduces the NAG function.  But it could be sped up by having it 
! calculate 64 random numbers at a time (and saving them in an array
! to be returned one call at a time).  This is the same speed trick
! NAG uses.  Also, it could be made more portable by splitting up the 
! 59-bit random numbers into 4 numbers of 15 bits each, so that each could
! be multiplied in 32 bit registers without overflow, and then doing
! the "carry" operation involved in multiplication by hand.  At present
! this algorithm relies on the fortran compiler supporting 64-bit integers
! and allowing integer multiplication to overflow without stopping and
! leaving the lowest order 64 bit intact (which is true on many computers,
! but perhaps not all).  The random number generators in "Numerical Recipes
! in Fortran 90" actually assume that 32 bit integer arithmetic does the
! "wrap-around" overflow of two's complement arithmetic "that is used on 
! virtually all modern machines" (NR p. 1144).  If we assume this also,
! then for machines that don't support 64 bit integers, we could split the
! 59-bit random numbers into 2 numbers of 30 bits each (instead of 4 numbers
! of 15 bits each as stated above).
!
  function g05cae_work(iseed_io)
    ! real function g05cae_work(iseed_io) generates a random number,
    ! uniformly distributed between 0 and 1.  iseed_io is an optional
    ! argument, which if set, resets the internal seed used.
    !
    ! Uses identical algorithm to NAG's g05cae. This is a linear (or 
    ! multiplicative in this limit) congruential algorithm due to Knuth:
    !
    !          z(i+1) = a * z(i) mod m
    !          r = z/m
    !
    ! with a=13**13,  m = 2**59
    !
    ! In order for this implementation to work, this computer must support
    ! 64-bit integers, and allow 64-bit integer multiplication to overflow 
    ! in such a way to leave the lowest order 64-bits (as is natural).
    !
    ! Written by Greg Hammett, March 6, 1999.

    implicit none
    integer, parameter :: r8 = selected_real_kind(12)
    real(r8) g05cae_work
    integer, intent(in), optional :: iseed_io
    ! need 8-byte/64-bit integers for this to work, 2**63=9.e18 (approx):
    integer, parameter :: i8 = selected_int_kind(18)

    ! local vars: 
    integer(i8), save :: iseed
    integer(i8) ibig, ibig_neg
    integer(i8), save :: a, m
    real(r8), save :: rscale
    logical, save :: first =.true.

    if(first) then
       a=13_i8**13
       m=2_i8**59
       rscale=1.0_r8/m
       iseed = 123456789_i8 * (2_i8**32+1)  ! default initial seed
       first=.false.

       ! test whether 64-bit integers have expected wrap-around properties:
       ibig=huge(1_i8)
       ibig_neg=-ibig-1
       if(ibig+1 .ne. ibig_neg)then
          print *, 'error in 64-bit random number generator in nag_sub.f90'
          stop
       endif
    endif
    if(present(iseed_io)) then
       iseed=2*iseed_io+1
       if(iseed < 1 .or. iseed >= m) then
          print *, "fatal error in g05cae_work, iseed out of range"
          stop
       endif
    endif

    ! The following would fail because iseed*a can overflow to be
    ! negative:
    ! iseed=mod(iseed*a,m)
    ! so instead do the following :

    iseed=iseed*a
    if(iseed < 0) iseed=iseed+2_i8**63  
    ! adding 2**63 turns off the 64th bit (which is a sign bit)
    iseed=mod(iseed,m)

    g05cae_work=iseed*rscale
    return
  end function g05cae_work
end module g05cae_mod

subroutine g05cbe(iseed_io)
  ! initialize random number seed to a repeatable sequence
  use g05cae_mod
  implicit none
  integer iseed_io
  integer, parameter :: r8 = selected_real_kind(12)
  real(r8) xdum
  xdum = g05cae_work(iseed_io)
  return
end subroutine g05cbe

function g05cae(xdum)
  ! real function g05cae() generates random number between 0.0 and 1.0
  use g05cae_mod
  implicit none
  integer, parameter :: r8 = selected_real_kind(12)
  real g05cae
  real, intent(in), optional :: xdum

  g05cae=g05cae_work()
  return
end function g05cae

! modified Bessel function I_0
real function s18cee(rkp2,ifail)
  implicit none
  real rkp2,xdum
  integer ifail
  call gamma_n(rkp2,s18cee,xdum)
  return
end function s18cee

! modified Bessel function I_1
real function s18cfe(rkp2,ifail)
  implicit none
  real rkp2, xdum
  integer ifail
  call gamma_n(rkp2,xdum,s18cfe)
  return
end function s18cfe

      subroutine gamma_n(x_LL,g0_LL,g1_LL)
      
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Evaluates modified Bessel functions
!
!cmnt  gamma_n(x,g0,g1) returns g0=exp(-abs(x))*I(0)(x)
!cmnt  gamma_n(x,g0,g1) returns g1=exp(-abs(x))*I(1)(x)
! (from Kerbel's port of the code).
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      
      integer k
      real tolerance
!      parameter(tolerance=1.0e-16)
      
      real xk,xki,xkp1i,xkp2i,error1,error0
      
      real x_LL
      real g0_LL
      real g1_LL
      
      real tk0_LL
      real tk1_LL
      real xb2_LL
      real xb2sq_LL
      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cmnt  compute Gamma_0 and Gamma_1
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      tolerance=epsilon(1.0)/8.0

      if (x_LL.eq.0.) then
        g0_LL=1.
        g1_LL=0.
        return
      endif

      xb2_LL=.5*x_LL  
      xb2sq_LL=xb2_LL**2
      tk0_LL=exp(-x_LL)
      tk1_LL=tk0_LL*xb2_LL
      g0_LL=tk0_LL
      g1_LL=tk1_LL
      
      xk=1.
      
 10   continue
      
      xki=1./xk
      xkp1i=1./(xk+1.)
      xkp2i=1./(xk+2.)
      tk0_LL=tk0_LL*xb2sq_LL*xki*xki
      tk1_LL=tk1_LL*xb2sq_LL*xki*xkp1i
      g0_LL=g0_LL+tk0_LL
      g1_LL=g1_LL+tk1_LL
      xk=xk+1.
      error0=abs(tk0_LL/g0_LL)
      error1=abs(tk1_LL/g1_LL)
      if (error0.gt.tolerance .or. error1.gt.tolerance) goto 10
            
      return
      end

! elliptic integrals:
real function s21bbe(arg1, arg2, arg3, ifail)
  implicit none
  real arg1, arg2, arg3, rfl
  integer ifail
  s21bbe = rfl(arg1,arg2,arg3,ifail)
  return
end function s21bbe

! elliptic integrals:
real function s21bce(arg1,arg2,arg3,ifail)
  implicit none
  real arg1, arg2, arg3, rdl
  integer ifail
  s21bce= rdl(arg1,arg2,arg3,ifail)
  return
end function s21bce

! complex-to-complex FFT
subroutine c06fce(rx,ry,n,work,ifail)
  use fft_work
  implicit none
  real rx(n), ry(n), work(n)
  integer n, ifail

  ! local vars:
  real scale
  complex zin(n),zout(n)
  real, allocatable :: table(:)
  real, allocatable :: cwork(:)
  integer, save :: n_first = -1

  scale=1.0/sqrt(float(n))
  zin = cmplx(rx,ry)

  if (n .ne. n_first) then
     ! initialize fft work arrays:
     if(.not. allocated(table)) allocate(table(ccfft_table0+ccfft_table1*n))
     if(.not. allocated(cwork)) allocate(cwork(ccfft_work0+ccfft_work1*n))
     call ccfft(0,n,scale,zin,zout,table,cwork,0)
     n_first=n
  endif

  call ccfft(-1,n,scale,zin,zout,table,cwork,0)
  rx=real(zout)
  ry=aimag(zout)

return
end

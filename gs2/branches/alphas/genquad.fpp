module genquad

! Allows a general weight function (as specified by an external function) for integration.
!   Returns the Gaussian quadrature weights and abscissae corresponding to the 
!   set of orthogonal polynomials associated with this weight.
! See: 
! Gnauchi, in: "Orthogonal Polynomials: Theory and Practice", Nevai (Ed.). Kluwer, 1990.
! Wheeler, Rocky Mountain Journal of Mathematics, Vol 4, p. 287 (1974)
!
! Written by: George Wilkie (gwilkie@umd.edu)
!

  implicit none

  public :: get_quadrature_rule

  private

  real,dimension(:),allocatable:: a_poly, b_poly
  logical:: inf_flag1, inf_flag2
  real:: epsfac = 1.0 
  integer:: poly_type, poly_order
 

contains
!> Accepts a weighting function (with only one argument), defined as external, 
!! calculates the modified moments, and returns the abscissae and weights for a new 
!! scheme based on polynomials orthogonal to the given weight function. The method is 
!! that of Wheeler. Also described in Numerical Recipes, 3rd ed. section 4.6.2 and 4.6.3.
  subroutine get_quadrature_rule ( &
    wgt_func, &        !< Input: Suitable integration weights as a single-argument function.
                       !!   Defined as "external" in calling routine!
    N_out,&            !< Input: The number of integration points desired
    v0,&               !< Input: The left boundary of the integration domain
    vf,&               !< Input: The right boundary of the integration domain
    out_abscissae, &   !< Output: The abscissae (quadrature points) for integration
    out_wgts, &        !< Output: The integration weights
    inf_flag1_in, &    !< Optional: Set to true to make lower bound -inf
    inf_flag2_in)      !< Optional: Set to true to make upper bound +inf

    use quadpack, only: qagi, qags

    implicit none
    real, external:: wgt_func
    integer,intent(in):: N_out
    real,intent(in):: v0, vf
    real,dimension(1:), intent(inout):: out_abscissae, out_wgts
    logical,optional,intent(in):: inf_flag1_in, inf_flag2_in
    real,dimension(0:N_out):: alpha,beta
    integer:: i, nerr, ier, neval, n
    real:: errest, moment1, moment2

#ifdef LAPACK

    if( present(inf_flag1_in)) then
       inf_flag1 = inf_flag1_in
    else
       inf_flag1 = .false.
    end if
    if( present(inf_flag2_in)) then
       inf_flag2 = inf_flag2_in
    else
       inf_flag2 = .false.
    end if

    ! Initialize arrays
    call init_general_quad(N_out)

    call set_poly_type()


    ! Calculate the recursion coefficients for the polynomials orthogonal with respect to 
    ! new_wgt_fcn 
    call calculate_new_coeffs(N_out)

!    b_poly(0) = 0.0

    if (poly_type .EQ. 0) then
       call qags( wgt_x,v0,vf,epsfac*epsilon(0.0),epsfac*epsilon(0.0),moment1,errest,neval,ier)
       call qags( wgt_func,v0,vf,epsfac*epsilon(0.0),epsfac*epsilon(0.0),moment2,errest,neval,ier)
       b_poly(0) = moment2
       a_poly(0) =  moment1/moment2
       poly_order = 1
       call qags( poly_sqx,v0,vf,epsfac*epsilon(0.0),epsfac*epsilon(0.0),moment1,errest,neval,ier)
       call qags( poly_sq,v0,vf,epsfac*epsilon(0.0),epsfac*epsilon(0.0),moment2,errest,neval,ier)
       a_poly(1) = moment1/moment2
       poly_order = 0
       call qags( poly_sq,v0,vf,epsfac*epsilon(0.0),epsfac*epsilon(0.0),moment1,errest,neval,ier)
       b_poly(1) = moment2/moment1
       do n = 2,N_out-1
          poly_order = n
          call qags( poly_sqx,v0,vf,epsfac*epsilon(0.0),epsfac*epsilon(0.0),moment1,errest,neval,ier)
          call qags( poly_sq,v0,vf,epsfac*epsilon(0.0),epsfac*epsilon(0.0),moment2,errest,neval,ier)
          a_poly(n) = moment1/moment2
          poly_order = n-1
          call qags( poly_sq,v0,vf,epsfac*epsilon(0.0),epsfac*epsilon(0.0),moment1,errest,neval,ier)
          b_poly(n) = moment2/moment1
       end do

    else if (poly_type .EQ. 1) then
       call qagi( wgt_x,v0,1,epsfac*epsilon(0.0),epsfac*epsilon(0.0),moment1,errest,neval,ier)
       call qagi( wgt_func,v0,1,epsfac*epsilon(0.0),epsfac*epsilon(0.0),moment2,errest,neval,ier)
       b_poly(0) = moment2
       a_poly(0) =  moment1/moment2
       poly_order = 1
       call qagi( poly_sqx,v0,1,epsfac*epsilon(0.0),epsfac*epsilon(0.0),moment1,errest,neval,ier)
       call qagi( poly_sq,v0,1,epsfac*epsilon(0.0),epsfac*epsilon(0.0),moment2,errest,neval,ier)

       a_poly(1) = moment1/moment2
       poly_order = 0
       call qagi( poly_sq,v0,1,epsfac*epsilon(0.0),epsfac*epsilon(0.0),moment1,errest,neval,ier)
       b_poly(1) = moment2/moment1
       do n = 2,N_out-1
          poly_order = n
          call qagi( poly_sqx,v0,1,epsfac*epsilon(0.0),epsfac*epsilon(0.0),moment1,errest,neval,ier)
          call qagi( poly_sq,v0,1,epsfac*epsilon(0.0),epsfac*epsilon(0.0),moment2,errest,neval,ier)
          a_poly(n) = moment1/moment2
          poly_order = n-1
          call qagi( poly_sq,v0,1,epsfac*epsilon(0.0),epsfac*epsilon(0.0),moment1,errest,neval,ier)
          b_poly(n) = moment2/moment1
       end do

    else if (poly_type .EQ. 2) then
       call qagi( wgt_x,v0,2,epsfac*epsilon(0.0),epsfac*epsilon(0.0),moment1,errest,neval,ier)
       call qagi( wgt_func,v0,2,epsfac*epsilon(0.0),epsfac*epsilon(0.0),moment2,errest,neval,ier)
       b_poly(0) = moment2
       a_poly(0) =  moment1/moment2
       poly_order = 1
       call qagi( poly_sqx,v0,2,epsfac*epsilon(0.0),epsfac*epsilon(0.0),moment1,errest,neval,ier)
       call qagi( poly_sq,v0,2,epsfac*epsilon(0.0),epsfac*epsilon(0.0),moment2,errest,neval,ier)
       a_poly(1) = moment1/moment2
       poly_order = 0
       call qagi( poly_sq,v0,2,epsfac*epsilon(0.0),epsfac*epsilon(0.0),moment1,errest,neval,ier)
       b_poly(1) = moment2/moment1
       do n = 2,N_out-1
          poly_order = n
          call qagi( poly_sqx,v0,2,epsfac*epsilon(0.0),epsfac*epsilon(0.0),moment1,errest,neval,ier)
          call qagi( poly_sq,v0,2,epsfac*epsilon(0.0),epsfac*epsilon(0.0),moment2,errest,neval,ier)
          a_poly(n) = moment1/moment2
          poly_order = n-1
          call qagi( poly_sq,v0,2,epsfac*epsilon(0.0),epsfac*epsilon(0.0),moment1,errest,neval,ier)
          b_poly(n) = moment2/moment1
       end do
    end if

    call solve_jacobi_matrix(N_out,out_abscissae,out_wgts)

    call finish_general_quad()

contains

   real function wgt_x(x)
     implicit none
     real,intent(in):: x
     wgt_x = wgt_func(x)*x
     return
   end function wgt_x

   real function poly_sqx(x)
     implicit none
     real,intent(in):: x
     real:: ym1,y,temp
     integer:: i,j

     ym1 = 0.0
     y = 1.0
     do i = 1,poly_order
       temp = y
       y = (x-a_poly(i-1))*y - b_poly(i-1)*ym1
       ym1 = temp
     end do 
     poly_sqx = wgt_func(x)*y*y*x
     return
   end function poly_sqx


   real function poly_sq(x)
     implicit none
     real,intent(in):: x
     real:: ym1,y,temp
     integer:: i,j

     ym1 = 0.0
     y = 1.0
     do i = 1,poly_order
       temp = y
       y = (x-a_poly(i-1))*y - b_poly(i-1)*ym1
       ym1 = temp
     end do 
     poly_sq = wgt_func(x)*y*y
     return
   end function poly_sq

#else
   write(*,*) "You must link to LAPACK for genquad to work."
   stop 1
#endif

  end subroutine get_quadrature_rule

#ifdef LAPACK
   subroutine set_poly_type()
     implicit none

     if ( (.NOT. inf_flag1) .AND. (.NOT. inf_flag2)) then
        ! Finite domain
        poly_type = 0
     else if ( (.NOT. inf_flag1) .AND. inf_flag2 ) then
        ! Semi-infinite domain
        poly_type = 1
     else if ( inf_flag1 .AND. (.NOT. inf_flag2) ) then
        ! Semi-infinite domain
        poly_type = -1
     else
       ! Infinite domain
       poly_type = 2
     end if
 
   end subroutine set_poly_type

!> Evaluates a polynomial given recursion coefficients
   real function evaluate_poly(a,b,n,x)
     implicit none
     real,dimension(0:),intent(in):: a,b
     real,intent(in):: x
     integer,intent(in):: n
     real:: y, ym1, temp
     integer:: i,j
     
     ym1 = 0.0
     y = 1.0
     do i = 1,n 
        temp = y
        y = (x-a(i-1))*y - b(i-1)*ym1
        ym1 = temp
     end do 

     evaluate_poly = y
     return

   end function evaluate_poly

!> Uses Wheeler's algorithm to obtain the recurrence relation for the new polynomials
!! given the recurrence relation and modified moments of a given set of polynomials (here
!! we use Legendre)
   subroutine calculate_new_coeffs(N)
     implicit none
     integer, intent(in):: N
     integer:: k
     real:: moment1, moment2

     b_poly(0) = 0.0 

     if (poly_type .EQ. 1) then
   
     end if

     contains


    end subroutine calculate_new_coeffs

   subroutine init_general_quad (N)
     implicit none
     integer,intent(in):: N
     
     allocate(a_poly(0:N-1)); a_poly = 0.0
     allocate(b_poly(0:N-1)); b_poly = 0.0

   end subroutine init_general_quad

!> Finds the eigenvalues and eigenvectors for a symmetric tridiagonal matrix and interprets
!! them as the absissae and weights for Gaussian quadrature. See Numerical Recipes, 3rd Ed. 
!! section 4.6.2. Uses the LAPACK routine, given as an external source file.
  subroutine solve_jacobi_matrix(N,abscissae,wgts)
    implicit none
    integer,intent(in):: N
    real,dimension(1:),intent(inout):: abscissae,wgts

    real*8, dimension(0:N-1):: a_use, b_use
    real*8 ,dimension(1:N,1:N):: eigenvectors   
    real*8,dimension(1:4*N):: workspace

    integer:: info, i

    b_use = sqrt(b_poly(0:N-1))
    a_use = a_poly(0:N-1)

    call DSTEQR('I',N,a_use(0:N-1),b_use(1:N-1),eigenvectors(1:N,1:N),N,workspace(1:4*N),info)

    if (info .NE. 0) then
       write(*,*) "ERROR: LAPACK returned error code ", info, " in genquad routine."
       write(*,*) "QUADPACK probably couldn't handle large numbers for resolution = ", N
       write(*,*) "Try compiling with quadruple (16-bit) precision."
       write(*,*) "diagonal = ", a_poly(0:N-1)
       write(*,*) "off-diagonal = sqrt ", b_poly(0:N-1)
       write(*,*) "output abscissae = ", a_use
       write(*,*) "output wights = ", (eigenvectors(1,N-i+1)**2,i=1,N)
       stop
    end if

    do i = 1,N 
       abscissae(i) = a_use(N-i)
       wgts(i) = b_poly(0)* eigenvectors(1,N-i+1)**2
    end do

  end subroutine

   subroutine finish_general_quad
     implicit none

     if(allocated(a_poly)) deallocate(a_poly)
     if(allocated(b_poly)) deallocate(b_poly)
   
   end subroutine finish_general_quad

#endif
end module genquad
  

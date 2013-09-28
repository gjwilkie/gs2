module general_quad

! Allows a general weight function (as specified on a table) for integration.
!   Returns the Gaussian quadrature weights and abscissae corresponding to the 
!   set of orthogonal polynomials associated with this weight.
! See: 
! Gnauchi, in: "Orthogonal Polynomials: Theory and Practice", Nevai (Ed.). Kluwer, 1990.
! Wheeler, Rocky Mountain Journal of Mathematics, Vol 4, p. 287 (1974)
!
! Written by: George Wilkie (gwilkie@umd.edu)
!

  implicit none

  public :: get_general_weights_from_grid
  public :: get_general_weights_from_moments

  private

  real,dimension(:),allocatable:: mod_moments
  integer:: N_overresolved
 

contains
!> Accepts an over-resolved integration scheme complete with abscissae and weights, 
!! calculates the modified moments, and returns the abscissae and weights for a new 
!! scheme based on polynomials orthogonal to the given weight function. The method is 
!! that of Wheeler. Also described in Numerical Recipes, 3rd ed. section 4.6.2 and 4.6.3.
  subroutine get_general_weights_from_grid ( &
    N_in, &            !< Input: Size of the input weight function table
    v_table_in,&       !< Input: Indepepdent variable, v, in the following two tables
                       !!   Assumed to be monotonically increasing.                     
    wgt_table_in, &    !< Input: Suitable integration weights corresponding to v_table
    new_wgt_fcn, &     !< Input: The new orthogonality weight function. Assumed to be 
                       !!   multiplied by whatever function will be integrated
    v0,&               !< Input: The left boundary of the integration domain
    vf,&               !< Input: The right boundary of the integration domain
    N_out,&            !< Input: The number of integration points desired
    out_abscissae, &   !< Output: The abscissae (quadrature points) for integration
    out_wgts)          !< Output: The integration weights

    implicit none
    integer,intent(in):: N_in, N_out
    real,dimension(1:),intent(in):: v_table_in, wgt_table_in, new_wgt_fcn
    real,dimension(1:), intent(inout):: out_abscissae, out_wgts
    real,dimension(0:2*N_out-1):: a_legendre,b_legendre
    real,dimension(0:N_out):: alpha,beta
    real,intent(in):: v0, vf
    integer:: i

    N_overresolved = N_in

#ifdef LAPACK
    ! Initialize arrays
    call init_general_quad(N_out)

    ! Calculate recursion coefficients for Legendre polynomials
    call populate_legendre_coeffs(2*N_out-1,a_legendre,b_legendre,v0,vf)

    call calculate_modified_moments(2*N_out-1,v_table_in,wgt_table_in, &
           new_wgt_fcn,a_legendre,b_legendre)

    ! Calculate the recursion coefficients for the polynomials orthogonal with respect to 
    ! new_wgt_fcn 
    call calculate_new_coeffs(N_out,a_legendre,b_legendre,alpha,beta)

    ! Find the zeroes and weights of the new polynomials
    call solve_jacobi_matrix(N_out,alpha(0:N_out-1),beta(0:N_out-1),out_abscissae,out_wgts)

    call finish_general_quad()
#else

    write (*,*) "You need LAPACK for this to work!"
    stop 1
#endif

  end subroutine get_general_weights_from_grid

!> Accepts the "modified moments":
!! Integral between (v0,vf): LegendreP[j,2*v/(vf-v0) + 1 - 2*vf/(vf-v0)] * F0(v)
!! where j ranges from 0 to 2*N_out-1
!! Returns the abscissae and weights for a new scheme based on polynomials
!! orthogonal to the given weight function. The method is that of Wheeler, as modified for
!! conditioning by Gnauchi. Also described in Numerical Recipes, 3rd ed. section 4.6.2 and 4.6.3.
  subroutine get_general_weights_from_moments ( &
    moments_in, &      !< Input: Modified moments of F0
    v0,&               !< Input: The left boundary of the integration domain
    vf,&               !< Input: The right boundary of the integration domain
    N_out,&            !< Input: The number of integration points desired
    out_abscissae, &   !< Output: The abscissae (quadrature points) for integration
    out_wgts)          !< Output: The integration weights

    implicit none
    integer,intent(in):: N_out
    real,dimension(0:),intent(in):: moments_in
    real,dimension(1:), intent(inout):: out_abscissae, out_wgts
    real,dimension(0:2*N_out-1):: a_legendre,b_legendre
    real,dimension(0:N_out):: alpha,beta
    real,intent(in):: v0, vf
    integer:: i

    mod_moments(0:2*N_out-1) = moments_in(0:2*N_out-1)

#ifdef LAPACK
    ! Initialize arrays
    call init_general_quad(N_out)

    ! Calculate recursion coefficients for Legendre polynomials
    call populate_legendre_coeffs(2*N_out-1,a_legendre,b_legendre,v0,vf)

    ! Calculate the recursion coefficients for the polynomials orthogonal with respect to 
    ! new_wgt_fcn 
    call calculate_new_coeffs(N_out,a_legendre,b_legendre,alpha,beta)

    ! Find the zeroes and weights of the new polynomials
    call solve_jacobi_matrix(N_out,alpha(0:N_out-1),beta(0:N_out-1),out_abscissae,out_wgts)

    call finish_general_quad()
#else

    write (*,*) "You need LAPACK for this to work!"
    stop 1
#endif

   end subroutine get_general_weights_from_moments

#ifdef LAPACK
!> Evaluates a polynomial given recursion coefficients
   real function evaluate_poly_recursion(a,b,n,x)
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

     evaluate_poly_recursion = y
     return

   end function evaluate_poly_recursion

!> Calculates the "modified moments" from Gnauchi. If w(v) is the desired integration measure,
!! and P_n are the Legendre polynomials (for instance), the modified moments are Integral( w(v)*P_n(v))
!! over the desired integration domain. Here, a and b are the monic recursion coefficients that define
!! the Legendre polynomials in the appropriate domain
   subroutine calculate_modified_moments(N,vtab,wgts,f,a,b)
     implicit none
     real,dimension(1:),intent(in):: vtab,wgts,f
     real,dimension(0:),intent(in):: a,b
     integer,intent(in):: N
     integer:: i, j
     real:: sum, arg, conv
     
     do i = 0,N
        mod_moments(i) = 0.0
        do j = 1,N_overresolved
           mod_moments(i) = mod_moments(i) + evaluate_poly_recursion(a,b,i,vtab(j))*f(j)*wgts(j)
        end do
     end do

   end subroutine calculate_modified_moments
!
    subroutine populate_legendre_coeffs(N,a,b,v0,vf)
      implicit none
      integer,intent(in):: N
      real,dimension(0:),intent(inout):: a,b
      real,intent(in)::v0,vf
      integer:: k

      a = 0.0 
        
      b(0) = 0.0
      do k = 1,N
         b(k) = 1.0/(4.0 - 1.0/real(k*k))
      end do

      ! The following is a required change to recursion coefficients due to the 
      ! shifting of the domain of integration
      a = (a + 2.0*vf/(vf-v0) - 1.0)*(vf-v0)*0.5
      b = b*((vf-v0)*0.5)**2
     
    end subroutine populate_legendre_coeffs

!> Uses Wheeler's algorithm to obtain the recurrence relation for the new polynomials
!! given the recurrence relation and modified moments of a given set of polynomials (here
!! we use Legendre)
   subroutine calculate_new_coeffs(N,a,b,alpha,beta)
     implicit none
     integer, intent(in):: N
     real,dimension(0:),intent(in):: a,b
     real,dimension(0:),intent(inout):: alpha,beta
  
     real,dimension(-1:N,0:2*N):: sigma
     integer:: k,l
 
     sigma(:,:) =  0.0
     sigma(0,0:2*N-1) = mod_moments(0:2*N-1)
     alpha(:) = 0.0
     alpha(0) = a(0) + mod_moments(1)/mod_moments(0)
     beta(:) = 0.0 

!     do k=1,(N+1)
     do k=1,N
        do l = k,2*N-k
           sigma(k,l) = sigma(k-1,l+1) - (alpha(k-1) - a(l))*sigma(k-1,l) & 
                       - beta(k-1)*sigma(k-2,l) + b(l)*sigma(k-1,l-1)
        end do
        alpha(k) = a(k) - sigma(k-1,k)/sigma(k-1,k-1) + sigma(k,k+1)/sigma(k,k)
        beta(k) = sigma(k,k)/sigma(k-1,k-1)
     end do

   end subroutine calculate_new_coeffs

   subroutine init_general_quad (N)
     implicit none
     integer,intent(in):: N
     
     allocate(mod_moments(0:2*(N+1)-1)); mod_moments = 0.0

   end subroutine init_general_quad

!> Finds the eigenvalues and eigenvectors for a symmetric tridiagonal matrix and interprets
!! them as the absissae and weights for Gaussian quadrature. See Numerical Recipes, 3rd Ed. 
!! section 4.6.2. Uses the 1 routine, given as an external source file.
  subroutine solve_jacobi_matrix(N,a,b,abscissae,wgts)
    implicit none
    integer,intent(in):: N
    real,dimension(0:),intent(inout):: a,b
    real,dimension(1:),intent(inout):: abscissae,wgts

    real, dimension(0:N-1):: a_use, b_use
    real ,dimension(1:N,1:N):: eigenvectors   
    real,dimension(1:4*N):: workspace

    integer:: info, i

    b_use = sqrt(b)
    a_use = a

    call DPTEQR('I',N,a_use(0:N-1),b_use(1:N-1),eigenvectors(1:N,1:N),N,workspace(1:4*N),info)

    if (info .NE. 0) then
       write(*,*) "LAPACK returned error code ", info
       stop
    end if

    do i = 1,N 
       abscissae(i) = a_use(N-i)
       wgts(i) = mod_moments(0)* eigenvectors(1,N-i+1)*eigenvectors(1,N-i+1)
    end do

  end subroutine

   subroutine finish_general_quad
     implicit none

     if(allocated(mod_moments)) deallocate(mod_moments)
   end subroutine finish_general_quad

#endif
end module general_quad
  

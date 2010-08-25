module vdimstub

contains
subroutine vdim (nr, md)

  integer nr, md

  nr = 1
  md = 1
  
end subroutine vdim

subroutine vdat(eqpsi_0, eqpsi_a, B_T, I_0, p_0, &
            r_min, r_max, Rmag, z_max, Zmag, avgrmid, &
            xp, pbar, fp, q, ba, &
            pres, cur, dvdrho, &
            area, psi, vB_psi, c)
  
  integer, parameter :: m1v = 1, mt = 1 

  real eqpsi_0, eqpsi_a
  real :: B_T, I_0, p_0, r_min, r_max, Rmag, z_max, Zmag, avgrmid

  real, dimension(m1v) :: xp, pbar, fp, q, ba, pres, cur, psi, dvdrho, area
  real vB_psi(m1v, mt)
  real c(m1v, 6)
  
  integer i, j

  eqpsi_0 = 0.
  eqpsi_a = 0.
  B_T = 0.
  I_0 = 0.
  p_0 = 0.
  r_min = 0.
  r_max = 0.
  Rmag = 0.
  z_max = 0.
  Zmag = 0.
  avgrmid = 0. 

  do i = 1,m1v
     xp(i) = 0.
     pbar(i) = 0.
     fp(i) = 0.
     q(i) = 0.
     ba(i) = 0.
     pres(i) = 0.
     cur(i) = 0.
     dvdrho(i) = 0.
     area(i) = 0.
     psi(i) = 0.
     c(i,1) = 0.
     c(i,2) = 0.
     c(i,3) = 0.
     c(i,4) = 0.
     c(i,5) = 0.
     c(i,6) = 0.
     do j= 1,mt
        vB_psi(i,j) = 0.
     enddo
  enddo

end subroutine vdat

end module vdimstub

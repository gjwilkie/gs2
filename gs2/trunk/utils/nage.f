      subroutine aazf06( srname, info )
c     mark 12 release. nag copyright 1986.
c     .. scalar arguments ..
      integer            info
      character*13       srname
c     ..
c
c  purpose
c  =======
c
c  aazf06  is an error handler for the level 2 blas routines.
c
c  it is called by the level 2 blas routines if an input parameter is
c  invalid.
c
c  parameters
c  ==========
c
c  srname - character*13.
c           on entry, srname specifies the name of the routine which
c           called aazf06.
c
c  info   - integer.
c           on entry, info specifies the position of the invalid
c           parameter in the parameter-list of the calling routine.
c
c
c  auxiliary routine for level 2 blas.
c
c  written on 20-july-1986.
c
c     .. local scalars ..
      integer            ifail
      character*80       rec( 1 )
c     .. external functions ..
      integer            p01abe
      external           p01abe
c     ..
c     .. executable statements ..
      write( rec( 1 ), 9999 )srname, info
      ifail = 0
      ifail = p01abe( ifail, -1, srname( 1: 6 ), 1, rec )
c
      return
c
 9999 format( ' ** on entry to ', a13, ' parameter number ', i2, ' had',
     $      ' an illegal value' )
c
c     end of aazf06.
c
      end
      subroutine abzp01
c     mark 11.5(f77) release. nag copyright 1986.
c
c     terminates execution when a hard failure occurs.
c
c     ******************** implementation note ********************
c     the following stop statement may be replaced by a call to an
c     implementation-dependent routine to display a message and/or
c     to abort the program.
c     *************************************************************
c     .. executable statements ..
      stop
      end
      subroutine agyd02(y,x,h,n,n1,aux,wspace,wspac1,param)
c
c     uses the technique of nag library
c     procedure d02aae with the
c     specification of aux changed
c     all implicitly declared reals may
c     be declared double precision
c     nag copyright 1975
c     mark 4.5 revised
c     mark 11.5(f77) revised. (sept 1985.)
c     .. scalar arguments ..
      real              h, x
      integer           n, n1
c     .. array arguments ..
      real              param(n1), wspac1(n), wspace(n,9), y(n)
c     .. subroutine arguments ..
      external          aux
c     .. local scalars ..
      real              c1, c2, dum, u, v, w
      integer           i
c     .. intrinsic functions ..
      intrinsic         abs
c     .. executable statements ..
      call aux(wspac1,y,x,param)
      c1 = 1.0/3.0
      c2 = 1.0/6.0
      do 20 i = 1, n
         wspace(i,3) = wspac1(i)
   20 continue
      u = c1*h
      do 40 i = 1, n
         wspace(i,4) = y(i)
         y(i) = y(i) + u*wspace(i,3)
   40 continue
      call aux(wspac1,y,u+x,param)
      do 60 i = 1, n
         wspace(i,5) = wspac1(i)
   60 continue
      v = h*c2
      do 80 i = 1, n
         y(i) = wspace(i,4) + v*(wspace(i,3)+wspace(i,5))
   80 continue
      call aux(wspac1,y,u+x,param)
      do 100 i = 1, n
         wspace(i,5) = wspac1(i)
  100 continue
      u = h*0.125
      v = h*0.375
      do 120 i = 1, n
         y(i) = wspace(i,4) + wspace(i,3)*u + wspace(i,5)*v
  120 continue
      u = 0.5*h
      v = 1.5*h
      w = 2.0*h
      call aux(wspac1,y,u+x,param)
      do 140 i = 1, n
         y(i) = wspace(i,4) + wspace(i,3)*u - wspace(i,5)*v + wspac1(i)
     *          *w
         wspace(i,5) = wspac1(i)
  140 continue
      x = x + h
      call aux(wspac1,y,x,param)
      u = h*c2
      v = 2.0*h*c1
      do 160 i = 1, n
         w = wspace(i,4) + u*(wspace(i,3)+wspac1(i)) + wspace(i,5)*v
         dum = w - y(i)
         wspace(i,2) = 0.2*abs(dum)
         y(i) = w
  160 continue
      return
      end
      subroutine agzd02(x,y,g,t,n,n1,m,h0,h,aux,wspace,wspac1,param)
c
c     uses the logic of nag library routine d02abe
c     eps and dum may be used double length.  eps is the
c     smallest real such that 1+eps>eps.  smax
c     is the largest integer. dum,err and hs maybe declared
c     double precision
c     nag copyright 1975
c     mark 4.5 revised
c     mark 11.5(f77) revised. (sept 1985.)
c     .. scalar arguments ..
      real              h, h0, x
      integer           m, n, n1, t
c     .. array arguments ..
      real              g(n), param(n1), wspac1(n), wspace(n,9), y(n)
c     .. subroutine arguments ..
      external          aux
c     .. local scalars ..
      real              dum, eps, err, hs, p, q
      integer           d, i, j, s, smax
c     .. external functions ..
      real              x02aae
      integer           x02bbe
      external          x02aae, x02bbe
c     .. external subroutines ..
      external          agyd02
c     .. intrinsic functions ..
      intrinsic         abs, aint, max, real, int
c     .. executable statements ..
      eps = x02aae(eps)
      smax = x02bbe(eps)
      m = 0
      dum = max(abs(x),abs(x+h0))
      if (abs(h0).le.eps*dum) return
      dum = h/h0
      if (abs(dum).gt.1.0 .or. h.eq.0.0) go to 20
      dum = abs(h0/h+0.9)
      if (dum.gt.real(smax)) go to 160
      h = h0/aint(dum)
      go to 40
   20 h = h0
   40 p = 1.0
      if (t.eq.3) p = 0.0
      q = 1.0
      if (t.eq.2) q = 0.0
      dum = h0/h + 0.1
      s = int(dum)
      hs = 1.0e-4*abs(h)
   60 do 80 i = 1, n
         wspace(i,9) = y(i)
   80 continue
      call agyd02(y,x,h,n,n1,aux,wspace,wspac1,param)
      d = 0
      do 100 i = 1, n
         err = g(i)*(p+q*abs(y(i)))
         if (wspace(i,2).gt.err) go to 120
         if (40.0*wspace(i,2).gt.err) d = 1
  100 continue
      s = s - 1
      if (s.eq.0) return
      dum = real(s)/2.0 + 0.1
      j = int(dum)*2
      if (d.ne.0 .or. j.ne.s) go to 60
      h = 2.0*h
      s = s/2
      go to 60
  120 x = x - h
      do 140 i = 1, n
         y(i) = wspace(i,9)
  140 continue
      if (s.gt.smax/2) go to 160
      s = 2*s
      h = 0.5*h
      if (abs(h).gt.hs) go to 60
  160 m = 1
      return
      end
      subroutine d02age(h,error,parerr,param,c,n,n1,m1,aux,bcaux,raaux,
     *                  prsol,mat,copy,wspace,wspac1,wspac2,ifail)
c     nag copyright 1975
c     mark 4.5 revised
c     mark 11.5(f77) revised. (sept 1985.)
c
c     solves a general boundary value
c     problem for n differential equations
c     in n1 parameters using a shooting
c     and matching technique.  eps is the
c     largest real variable such that 1+eps=1
c     all implicitly declared reals may be used double-length
c     the array copy is redundant
c     .. parameters ..
      character*6       srname
      parameter         (srname='d02age')
c     .. scalar arguments ..
      real              h
      integer           ifail, m1, n, n1
c     .. array arguments ..
      real              c(m1,n), copy(n1,n1), error(n), mat(n1,n1),
     *                  param(n1), parerr(n1), wspac1(n), wspac2(n),
     *                  wspace(n,9)
c     .. subroutine arguments ..
      external          aux, bcaux, prsol, raaux
c     .. local scalars ..
      real              c1, d, dist, dum, eps, h0, oldres, pert, r,
     *                  resid, x, x1
      integer           count, count1, ct, em, i, itest, j, k, m, one
c     .. local arrays ..
      character*1       p01rec(1)
c     .. external functions ..
      real              x02aae
      integer           p01abe
      external          x02aae, p01abe
c     .. external subroutines ..
      external          agzd02, f03afe, f04aje
c     .. intrinsic functions ..
      intrinsic         abs, real
c     .. executable statements ..
      eps = x02aae(eps)
      m = m1 - 1
      if (n1.le.n) go to 20
      em = 1
      go to 940
   20 count = 0
      h0 = h
      one = 1
      em = -1
c
c     forms the residuals at the
c     matching point
   40 call raaux(x,x1,r,param)
      if ((x-r)*(x1-r).le.0.0) go to 60
      em = 3
      go to 940
   60 if (h0*(x1-x).lt.0.0) h0 = -h0
      call bcaux(wspac1,wspac2,param)
      h = h0
      do 80 i = 1, n
         wspace(i,8) = wspac2(i)
   80 continue
      i = 1
      call agzd02(x,wspac1,error,one,n,n1,i,r-x,h,aux,wspace,wspac2,
     *            param)
      if (i.eq.0) go to 100
      em = 4
      go to 940
  100 do 120 i = 1, n
         dum = wspace(i,8)
         wspace(i,8) = -wspac1(i)
         wspac1(i) = dum
  120 continue
      h = -h0
      i = 1
      call agzd02(x1,wspac1,error,one,n,n1,i,r-x1,h,aux,wspace,wspac2,
     *            param)
      if (i.eq.0) go to 140
      em = 4
      go to 940
  140 resid = 0.0
      ct = 0
      do 160 i = 1, n1
         d = wspac1(i)
         dum = wspace(i,8)
         wspace(i,8) = d + dum
         d = 1.0 + abs(dum) + abs(d)
         dum = wspace(i,8)
         if (abs(dum).lt.error(i)*d) ct = ct + 1
         wspac1(i) = dum
         resid = resid + dum*dum
  160 continue
      call prsol(param,resid,n1,wspac1)
      if (em.ne.-1) go to 620
  180 count = count + 1
      if (count.ne.12) go to 200
      em = 7
      go to 940
c
c     forms the jacobian by numerical
c     differentiation
  200 do 360 k = 1, n1
         pert = 10.0*parerr(k)*(1.0+abs(param(k)))
         param(k) = pert + param(k)
         call raaux(x,x1,r,param)
         if ((x-r)*(x1-r).le.0.0) go to 220
         em = 3
         go to 940
  220    if (h0*(x1-x).lt.0.0) h0 = -h0
         h = h0
         call bcaux(wspac1,wspac2,param)
         do 240 i = 1, n
            wspace(i,7) = wspac2(i)
  240    continue
         i = 1
         call agzd02(x,wspac1,error,one,n,n1,i,r-x,h,aux,wspace,wspac2,
     *               param)
         if (i.eq.0) go to 260
         em = 2
         go to 940
  260    do 280 i = 1, n1
            mat(i,k) = wspac1(i)
  280    continue
         h = -h0
         i = 1
         do 300 i = 1, n
            wspac1(i) = wspace(i,7)
  300    continue
         call agzd02(x1,wspac1,error,one,n,n1,i,r-x1,h,aux,wspace,
     *               wspac2,param)
         if (i.eq.0) go to 320
         em = 2
         go to 940
  320    do 340 i = 1, n1
            mat(i,k) = (mat(i,k)-wspac1(i)+wspace(i,8))/pert
            if (abs(mat(i,k)).lt.5.0*eps*abs(wspace(i,8))/pert) mat(i,k)
     *           = 0.0
  340    continue
         param(k) = param(k) - pert
  360 continue
      itest = 1
      em = -3
c
c     performs column scaling on the jacobian
c     and forms a triangular decomposition
      do 440 i = 1, n1
         d = 0.0
         do 380 j = 1, n1
            if (abs(mat(j,i)).gt.d) d = abs(mat(j,i))
  380    continue
         if (d.ne.0.0) go to 400
         em = 5
         go to 940
  400    do 420 j = 1, n1
            mat(j,i) = mat(j,i)/d
  420    continue
         wspace(i,7) = d
  440 continue
      i = 1
      call f03afe(n1,eps,mat,n1,d,j,wspac1,i)
      if (i.eq.0) go to 460
      em = 5
      go to 940
  460 do 480 i = 1, n1
         wspace(i,6) = wspac1(i)
  480 continue
c
c     uses a generalised newton raphson
c     technique to solve the nonlinear
c     equations at the matching point
  500 oldres = resid
      count1 = 0
      do 520 i = 1, n1
         wspac1(i) = wspace(i,6)
         wspace(i,1) = wspace(i,8)
  520 continue
      call f04aje(n1,one,mat,n1,wspac1,wspace,n)
      do 540 i = 1, n1
         wspace(i,1) = wspace(i,1)/wspace(i,7)
         param(i) = param(i) + wspace(i,1)
  540 continue
      if (ct.lt.n1) go to 580
      do 560 i = 1, n1
         if (abs(wspace(i,1)).gt.parerr(i)*(1.0+abs(param(i))))
     *       go to 580
  560 continue
      em = -5
      go to 740
  580 do 600 i = 1, n1
         wspace(i,1) = -wspace(i,1)
  600 continue
      go to 40
  620 if (count1.ne.0) go to 640
      if (resid.ge.oldres/10.0) go to 640
      em = -2
      itest = 0
      go to 500
  640 if (resid.lt.oldres) go to 180
      if (count1.ne.3) go to 700
      if (itest.eq.0) go to 660
      em = 6
      go to 940
  660 do 680 i = 1, n1
         param(i) = param(i) + wspace(i,1)
  680 continue
      em = -1
      go to 40
  700 count1 = count1 + 1
      em = -4
      do 720 i = 1, n1
         wspace(i,1) = wspace(i,1)/2.0
         param(i) = param(i) + wspace(i,1)
  720 continue
      go to 40
c
c     calculates the final solution
  740 if (m.le.0) go to 940
      call raaux(x,x1,r,param)
      if ((x-r)*(x1-r).le.0.0) go to 760
      em = 3
      go to 940
  760 if (h0*(x1-x).lt.0.0) h0 = -h0
      h = h0
      call bcaux(wspac1,wspac2,param)
      do 780 i = 1, n
         wspace(i,7) = wspac2(i)
  780 continue
      dist = (x1-x)/real(m)
      j = 1
      c1 = x
      k = 1
  800 do 820 i = 1, n
         c(j,i) = wspac1(i)
  820 continue
      if ((r-c1-0.25*dist)*dist.le.0.0) go to 860
      i = 1
      call agzd02(c1,wspac1,error,one,n,n1,i,dist,h,aux,wspace,wspac2,
     *            param)
      if (i.eq.0) go to 840
      em = 4
      go to 940
  840 j = j + k
      go to 800
  860 if (k.eq.-1) go to 900
      dist = -dist
      c1 = x1
      h = -h0
      j = m1
      k = -1
      do 880 i = 1, n
         wspac1(i) = wspace(i,7)
  880 continue
      go to 800
  900 call bcaux(wspac1,wspac2,param)
      do 920 i = 1, n
         c(1,i) = wspac1(i)
         c(m1,i) = wspac2(i)
  920 continue
  940 if (em.le.0) ifail = 0
      if (em.gt.0) ifail = p01abe(ifail,em,srname,0,p01rec)
      return
c     end of d02age
      end
      subroutine f03afe(n,eps,a,ia,d1,id,p,ifail)
c     mark 2 release. nag copyright 1972
c     mark 3 revised.
c     mark 4.5 revised
c     mark 11 revised. vectorisation (jan 1984).
c     mark 11.5(f77) revised. (sept 1985.)
c     mark 12 revised (level 2 blas) (mar 1986)
c     mark 12 revised. extended blas (june 1986)
c
c     unsymdet
c     the unsymmetric matrix, a, is stored in the n*n array a(i,j),
c     i=1,n, j=1,n. the decomposition a=lu, where l is a
c     lower triangular matrix and u is a unit upper triangular
c     matrix, is performed and overwritten on a, omitting the unit
c     diagonal of u. a record of any interchanges made to the rows
c     of a is kept in p(i), i=1,n, such that the i-th row and
c     the p(i)-th row were interchanged at the i-th step. the
c     determinant, d1 * 2.0**id, of a is also computed. the
c     subroutine
c     will fail if a, modified by the rounding errors, is singular
c     or almost singular. sets ifail = 0 if successful else ifail =
c     1.
c     1st december 1971
c
c     .. parameters ..
      character*6       srname
      parameter         (srname='f03afe')
c     .. scalar arguments ..
      real              d1, eps
      integer           ia, id, ifail, n
c     .. array arguments ..
      real              a(ia,n), p(n)
c     .. local scalars ..
      real              x, y
      integer           i, j, k, l
c     .. local arrays ..
      character*1       p01rec(1)
c     .. external functions ..
      integer           p01abe
      external          p01abe
c     .. external subroutines ..
      external          sgemv, strsv
c     .. intrinsic functions ..
      intrinsic         abs, sqrt
c     .. executable statements ..
      do 20 i = 1, n
         p(i) = 0.0
   20 continue
      do 60 j = 1, n
         do 40 i = 1, n
            p(i) = p(i) + a(i,j)**2
   40    continue
   60 continue
      do 80 i = 1, n
         if (p(i).le.0.0) go to 240
         p(i) = 1.0/sqrt(p(i))
   80 continue
      d1 = 1.0
      id = 0
      do 220 k = 1, n
         l = k
         x = 0.0
         do 100 i = k, n
            y = abs(a(i,k)*p(i))
            if (y.le.x) go to 100
            x = y
            l = i
  100    continue
         if (l.eq.k) go to 140
         d1 = -d1
         do 120 j = 1, n
            y = a(k,j)
            a(k,j) = a(l,j)
            a(l,j) = y
  120    continue
         p(l) = p(k)
  140    p(k) = l
         d1 = d1*a(k,k)
         if (x.lt.8.0*eps) go to 240
  160    if (abs(d1).lt.1.0) go to 180
         d1 = d1*0.0625
         id = id + 4
         go to 160
  180    if (abs(d1).ge.0.0625) go to 200
         d1 = d1*16.0
         id = id - 4
         go to 180
  200    if (k.lt.n) then
            call strsv('l','n','n',k,a,ia,a(1,k+1),1)
            call sgemv('n',n-k,k,-1.0,a(k+1,1),ia,a(1,k+1),1,1.0,a(k+1,
     *                 k+1),1)
         end if
  220 continue
      ifail = 0
      return
  240 ifail = p01abe(ifail,1,srname,0,p01rec)
      return
      end
      subroutine f04aje(n,ir,a,ia,p,b,ib)
c     mark 2 release. nag copyright 1972
c     mark 4 revised.
c     mark 4.5 revised
c     mark 11 revised. vectorisation (jan 1984).
c     mark 11.5(f77) revised. (sept 1985.)
c     mark 12 revised. extended blas (june 1986)
c
c     unsymsol
c     solves ax=b, where a is an unsymmetric matrix and b is an
c     n*ir
c     matrix of ir right-hand sides. the subroutine f04aje must by
c     preceded by f03afe in which l and u are produced in a(i,j),
c     from a, and the record of the interchanges is produced in
c     p(i). ax=b is solved in three steps, interchange the
c     elements of b, ly=b and ux=y. the matrices y and then x are
c     overwritten on b.
c     1st august 1971
c
c     .. scalar arguments ..
      integer           ia, ib, ir, n
c     .. array arguments ..
      real              a(ia,n), b(ib,ir), p(n)
c     .. local scalars ..
      real              x
      integer           i, i1, k
c     .. external subroutines ..
      external          strsv
c     .. executable statements ..
c     interchanging elements of b
      do 40 i = 1, n
         i1 = p(i) + 0.5
         if (i1.eq.i) go to 40
         do 20 k = 1, ir
            x = b(i,k)
            b(i,k) = b(i1,k)
            b(i1,k) = x
   20    continue
   40 continue
      do 60 k = 1, ir
c        solution of ly= b
         call strsv('l','n','n',n,a,ia,b(1,k),1)
c        solution of ux= y
         call strsv('u','n','u',n,a,ia,b(1,k),1)
   60 continue
      return
      end
      subroutine f06pae( trans, m, n, alpha, a, lda, x, incx, beta, y,
     $                   incy )
c     mark 12 release. nag copyright 1986.
c     .. entry points ..
      entry              sgemv( trans, m, n, alpha, a, lda, x, incx,
     $                   beta, y, incy )
c     .. scalar arguments ..
      real               alpha, beta
      integer            incx, incy, lda, m, n
      character*1        trans
c     .. array arguments ..
      real               a( lda, * ), x( * ), y( * )
c     ..
c
c  purpose
c  =======
c
c  sgemv  performs one of the matrix-vector operations
c
c     y := alpha*a*x + beta*y,   or   y := alpha*a'*x + beta*y,
c
c  where alpha and beta are scalars, x and y are vectors and a is an
c  m by n matrix.
c
c  parameters
c  ==========
c
c  trans  - character*1.
c           on entry, trans specifies the operation to be performed as
c           follows:
c
c              trans = 'n' or 'n'   y := alpha*a*x + beta*y.
c
c              trans = 't' or 't'   y := alpha*a'*x + beta*y.
c
c              trans = 'c' or 'c'   y := alpha*a'*x + beta*y.
c
c           unchanged on exit.
c
c  m      - integer.
c           on entry, m specifies the number of rows of the matrix a.
c           m must be at least zero.
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the number of columns of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  alpha  - real            .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  a      - real             array of dimension ( lda, n ).
c           before entry, the leading m by n part of the array a must
c           contain the matrix of coefficients.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           max( 1, m ).
c           unchanged on exit.
c
c  x      - real             array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ) when trans = 'n' or 'n'
c           and at least
c           ( 1 + ( m - 1 )*abs( incx ) ) otherwise.
c           before entry, the incremented array x must contain the
c           vector x.
c           unchanged on exit.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c  beta   - real            .
c           on entry, beta specifies the scalar beta. when beta is
c           supplied as zero then y need not be set on input.
c           unchanged on exit.
c
c  y      - real             array of dimension at least
c           ( 1 + ( m - 1 )*abs( incy ) ) when trans = 'n' or 'n'
c           and at least
c           ( 1 + ( n - 1 )*abs( incy ) ) otherwise.
c           before entry with beta non-zero, the incremented array y
c           must contain the vector y. on exit, y is overwritten by the
c           updated vector y.
c
c  incy   - integer.
c           on entry, incy specifies the increment for the elements of
c           y. incy must not be zero.
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 20-july-1986.
c     jack dongarra, argonne national laboratory.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      real               one, zero
      parameter          ( one = 1.0e+0, zero = 0.0e+0 )
c     .. local scalars ..
      real               temp
      integer            i, info, ix, iy, j, jx, jy, kx, ky, lenx, leny
c     .. external subroutines ..
      external           aazf06
c     .. intrinsic functions ..
      intrinsic          max
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if( ( trans.ne.'n' ).and.( trans.ne.'n' ).and.
     $    ( trans.ne.'t' ).and.( trans.ne.'t' ).and.
     $    ( trans.ne.'c' ).and.( trans.ne.'c' ) )then
         info = 1
      else if( m.lt.0 )then
         info = 2
      else if( n.lt.0 )then
         info = 3
      else if( lda.lt.max( 1, m ) )then
         info = 6
      else if( incx.eq.0 )then
         info = 8
      else if( incy.eq.0 )then
         info = 11
      end if
      if( info.ne.0 )then
         call aazf06( 'f06pae/sgemv ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( m.eq.0 ).or.( n.eq.0 ).or.
     $    ( ( alpha.eq.zero ).and.( beta.eq.one ) ) )return
c
c     set lenx and leny, the lengths of the vectors x and y.
c
      if( ( trans.eq.'n' ).or.( trans.eq.'n' ) )then
         lenx = n
         leny = m
      else
         lenx = m
         leny = n
      end if
c
c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through a.
c
c     first form  y := beta*y  and set up the start points in x and y if
c     the increments are not both unity.
c
      if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
         if( beta.ne.one )then
            if( beta.eq.zero )then
               do 10, i = 1, leny
                  y( i ) = zero
   10          continue
            else
               do 20, i = 1, leny
                  y( i ) = beta*y( i )
   20          continue
            end if
         end if
      else
         if( incx.gt.0 )then
            kx = 1
         else
            kx = 1 - ( lenx - 1 )*incx
         end if
         if( incy.gt.0 )then
            ky = 1
         else
            ky = 1 - ( leny - 1 )*incy
         end if
         if( beta.ne.one )then
            iy = ky
            if( beta.eq.zero )then
               do 30, i = 1, leny
                  y( iy ) = zero
                  iy = iy + incy
   30          continue
            else
               do 40, i = 1, leny
                  y( iy ) = beta*y( iy )
                  iy = iy + incy
   40          continue
            end if
         end if
      end if
      if( alpha.eq.zero )
     $   return
      if( ( trans.eq.'n' ).or.( trans.eq.'n' ) )then
c
c        form  y := alpha*a*x + y.
c
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 60, j = 1, n
               if( x( j ).ne.zero )then
                  temp = alpha*x( j )
                  do 50, i = 1, m
                     y( i ) = y( i ) + temp*a( i, j )
   50             continue
               end if
   60       continue
         else
            jx = kx
            do 80, j = 1, n
               if( x( jx ).ne.zero )then
                  temp = alpha*x( jx )
                  iy = ky
                  do 70, i = 1, m
                     y( iy ) = y( iy ) + temp*a( i, j )
                     iy = iy + incy
   70             continue
               end if
               jx = jx + incx
   80       continue
         end if
      else
c
c        form  y := alpha*a'*x + y.
c
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 100, j = 1, n
               temp = zero
               do 90, i = 1, m
                  temp = temp + a( i, j )*x( i )
   90          continue
               y( j ) = y( j ) + alpha*temp
  100       continue
         else
            jy = ky
            do 120, j = 1, n
               temp = zero
               ix = kx
               do 110, i = 1, m
                  temp = temp + a( i, j )*x( ix )
                  ix = ix + incx
  110          continue
               y( jy ) = y( jy ) + alpha*temp
               jy = jy + incy
  120       continue
         end if
      end if
c
      return
c
c     end of f06pae (sgemv ).
c
      end
      subroutine f06pje( uplo, trans, diag, n, a, lda, x, incx )
c     mark 12 release. nag copyright 1986.
c     .. entry points ..
      entry              strsv( uplo, trans, diag, n, a, lda, x, incx )
c     .. scalar arguments ..
      integer            incx, lda, n
      character*1        diag, trans, uplo
c     .. array arguments ..
      real               a( lda, * ), x( * )
c     ..
c
c  purpose
c  =======
c
c  strsv  solves one of the systems of equations
c
c     a*x = b,   or   a'*x = b,
c
c  where b and x are n element vectors and a is an n by n unit, or
c  non-unit, upper or lower triangular matrix.
c
c  no test for singularity or near-singularity is included in this
c  routine. such tests must be performed before calling this routine.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the matrix is an upper or
c           lower triangular matrix as follows:
c
c              uplo = 'u' or 'u'   a is an upper triangular matrix.
c
c              uplo = 'l' or 'l'   a is a lower triangular matrix.
c
c           unchanged on exit.
c
c  trans  - character*1.
c           on entry, trans specifies the equations to be solved as
c           follows:
c
c              trans = 'n' or 'n'   a*x = b.
c
c              trans = 't' or 't'   a'*x = b.
c
c              trans = 'c' or 'c'   a'*x = b.
c
c           unchanged on exit.
c
c  diag   - character*1.
c           on entry, diag specifies whether or not a is unit
c           triangular as follows:
c
c              diag = 'u' or 'u'   a is assumed to be unit triangular.
c
c              diag = 'n' or 'n'   a is not assumed to be unit
c                                  triangular.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  a      - real             array of dimension ( lda, n ).
c           before entry with  uplo = 'u' or 'u', the leading n by n
c           upper triangular part of the array a must contain the upper
c           triangular matrix and the strictly lower triangular part of
c           a is not referenced.
c           before entry with uplo = 'l' or 'l', the leading n by n
c           lower triangular part of the array a must contain the lower
c           triangular matrix and the strictly upper triangular part of
c           a is not referenced.
c           note that when  diag = 'u' or 'u', the diagonal elements of
c           a are not referenced either, but are assumed to be unity.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           max( 1, n ).
c           unchanged on exit.
c
c  x      - real             array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the n
c           element right-hand side vector b. on exit, x is overwritten
c           with the solution vector x.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 20-july-1986.
c     jack dongarra, argonne national laboratory.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      real               zero
      parameter          ( zero = 0.0e+0 )
c     .. local scalars ..
      integer            i, info, ix, j, jx, kx
      logical            nounit
c     .. external subroutines ..
      external           aazf06
c     .. intrinsic functions ..
      intrinsic          max
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if( ( uplo.ne.'u' ).and.( uplo.ne.'u' ).and.( uplo.ne.'l' ).and.
     $    ( uplo.ne.'l' ) )then
         info = 1
      else if( ( trans.ne.'n' ).and.( trans.ne.'n' ).and.
     $         ( trans.ne.'t' ).and.( trans.ne.'t' ).and.
     $         ( trans.ne.'c' ).and.( trans.ne.'c' ) )then
         info = 2
      else if( ( diag.ne.'u' ).and.( diag.ne.'u' ).and.
     $         ( diag.ne.'n' ).and.( diag.ne.'n' ) )then
         info = 3
      else if( n.lt.0 )then
         info = 4
      else if( lda.lt.max( 1, n ) )then
         info = 6
      else if( incx.eq.0 )then
         info = 8
      end if
      if( info.ne.0 )then
         call aazf06( 'f06pje/strsv ', info )
         return
      end if
c
c     quick return if possible.
c
      if( n.eq.0 )
     $   return
c
      nounit = ( diag.eq.'n' ).or.( diag.eq.'n' )
c
c     set up the start point in x if the increment is not unity. this
c     will be  ( n - 1 )*incx  too small for descending loops.
c
      if( incx.le.0 )then
         kx = 1 - ( n - 1 )*incx
      else if( incx.ne.1 )then
         kx = 1
      end if
c
c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through a.
c
      if( ( trans.eq.'n' ).or.( trans.eq.'n' ) )then
c
c        form  x := inv( a )*x.
c
         if( ( uplo.eq.'u' ).or.( uplo.eq.'u' ) )then
            if( incx.eq.1 )then
               do 20, j = n, 1, -1
                  if( x( j ).ne.zero )then
                     if( nounit )
     $                  x( j ) = x( j )/a( j, j )
                     do 10, i = j - 1, 1, -1
                        x( i ) = x( i ) - x( j )*a( i, j )
   10                continue
                  end if
   20          continue
            else
               jx = kx + ( n - 1 )*incx
               do 40, j = n, 1, -1
                  if( x( jx ).ne.zero )then
                     if( nounit )
     $                  x( jx ) = x( jx )/a( j, j )
                     ix = jx
                     do 30, i = j - 1, 1, -1
                        ix = ix - incx
                        x( ix ) = x( ix ) - x( jx )*a( i, j )
   30                continue
                  end if
                  jx = jx - incx
   40          continue
            end if
         else
            if( incx.eq.1 )then
               do 60, j = 1, n
                  if( x( j ).ne.zero )then
                     if( nounit )
     $                  x( j ) = x( j )/a( j, j )
                     do 50, i = j + 1, n
                        x( i ) = x( i ) - x( j )*a( i, j )
   50                continue
                  end if
   60          continue
            else
               jx = kx
               do 80, j = 1, n
                  if( x( jx ).ne.zero )then
                     if( nounit )
     $                  x( jx ) = x( jx )/a( j, j )
                     ix = jx
                     do 70, i = j + 1, n
                        ix = ix + incx
                        x( ix ) = x( ix ) - x( jx )*a( i, j )
   70                continue
                  end if
                  jx = jx + incx
   80          continue
            end if
         end if
      else
c
c        form  x := inv( a' )*x.
c
         if( ( uplo.eq.'u' ).or.( uplo.eq.'u' ) )then
            if( incx.eq.1 )then
               do 100, j = 1, n
                  do 90, i = 1, j - 1
                     x( j ) = x( j ) - a( i, j )*x( i )
   90             continue
                  if( nounit )
     $               x( j ) = x( j )/a( j, j )
  100          continue
            else
               jx = kx
               do 120, j = 1, n
                  ix = kx
                  do 110, i = 1, j - 1
                     x( jx ) = x( jx ) - a( i, j )*x( ix )
                     ix = ix + incx
  110             continue
                  if( nounit )
     $               x( jx ) = x( jx )/a( j, j )
                  jx = jx + incx
  120          continue
            end if
         else
            if( incx.eq.1 )then
               do 140, j = n, 1, -1
                  do 130, i = n, j + 1, -1
                     x( j ) = x( j ) - a( i, j )*x( i )
  130             continue
                  if( nounit )
     $               x( j ) = x( j )/a( j, j )
  140          continue
            else
               kx = kx + ( n - 1 )*incx
               jx = kx
               do 160, j = n, 1, -1
                  ix = kx
                  do 150, i = n, j + 1, -1
                     x( jx ) = x( jx ) - a( i, j )*x( ix )
                     ix = ix - incx
  150             continue
                  if( nounit )
     $               x( jx ) = x( jx )/a( j, j )
                  jx = jx - incx
  160          continue
            end if
         end if
      end if
c
      return
c
c     end of f06pje (strsv ).
c
      end
      integer function p01abe(ifail,ierror,srname,nrec,rec)
c     mark 11.5(f77) release. nag copyright 1986.
c
c     p01abe is the error-handling routine for the nag library.
c
c     p01abe either returns the value of ierror through the routine
c     name (soft failure), or terminates execution of the program
c     (hard failure). diagnostic messages may be output.
c
c     if ierror = 0 (successful exit from the calling routine),
c     the value 0 is returned through the routine name, and no
c     message is output
c
c     if ierror is non-zero (abnormal exit from the calling routine),
c     the action taken depends on the value of ifail.
c
c     ifail = 1: soft failure, silent exit (i.e. no messages are output)
c     ifail =-1: soft failure, noisy exit (i.e. messages are output)
c     ifail = 0: hard failure, noisy exit
c
c     for compatibility with certain routines included before mark 12
c     p01abe also allows an alternative specification of ifail in which
c     it is regarded as a decimal integer with least significant digits
c     cba. then
c
c     a = 0: hard failure  a = 1: soft failure
c     b = 0: silent exit   b = 1: noisy exit
c
c     except that hard failure now always implies a noisy exit.
c
c     s.hammarling, m.p.hooper and j.j.du croz, nag central office.
c
c     .. scalar arguments ..
      integer                 ierror, ifail, nrec
      character*(*)           srname
c     .. array arguments ..
      character*(*)           rec(*)
c     .. local scalars ..
      integer                 i, nerr
      character*72            mess
c     .. external subroutines ..
      external                abzp01, x04aae, x04bae
c     .. intrinsic functions ..
      intrinsic               abs, mod
c     .. executable statements ..
      if (ierror.ne.0) then
c        abnormal exit from calling routine
         if (ifail.eq.-1 .or. ifail.eq.0 .or.
     *       (ifail.gt.0 .and. mod(ifail/10,10).ne.0)) then
c           noisy exit
            call x04aae(0,nerr)
            do 20 i = 1, nrec
               call x04bae(nerr,rec(i))
   20       continue
            write (mess,fmt=99999) srname, ierror
            call x04bae(nerr,mess)
            if (abs(mod(ifail,10)).ne.1) then
c              hard failure
               call x04bae(nerr,
     *                     ' ** nag hard failure - execution terminated'
     *                     )
               call abzp01
            else
c              soft failure
               call x04bae(nerr,
     *                     ' ** nag soft failure - control returned')
            end if
         end if
      end if
      p01abe = ierror
      return
c
99999 format (' ** abnormal exit from nag library routine ',a,': ifail',
     *       ' =',i6)
      end
      real function x02aae(x)
c     mark 12 re-issue. nag copyright 1986.
c
c     returns  (1/2)*b**(1-p)  if rounds is .true.
c     returns  b**(1-p)  otherwise
c     i.e. returns the same value as x02aje
c
c     .. scalar arguments ..
      real                 x
c     .. external functions ..
      real                 x02aje
      external             x02aje
c     .. executable statements ..
      x02aae = x02aje()
      return
      end
      real function x02abe(x)
c     mark 12 re-issue. nag copyright 1986.
c
c     returns  b**(emin-1)  (the smallest positive model number)
c     i.e. returns the same value as x02ake
c
c     .. scalar arguments ..
      real                 x
c     .. external functions ..
      real                 x02ake
      external             x02ake
c     .. executable statements ..
      x02abe = x02ake()
      return
      end
      real function x02ace(x)
c     mark 12 re-issue. nag copyright 1986.
c
c     returns  (1 - b**(-p)) * b**emax  (the largest positive model
c     number)
c     i.e. returns the same value as x02ale
c
c     .. scalar arguments ..
      real                 x
c     .. external functions ..
      real                 x02ale
      external             x02ale
c     .. executable statements ..
      x02ace = x02ale()
      return
      end
      real function x02ade(x)
c     mark 12 re-issue. nag copyright 1986.
c
c     returns  x02ake/x02aje
c
c     .. scalar arguments ..
      real                 x
c     .. external functions ..
      real                 x02aje, x02ake
      external             x02aje, x02ake
c     .. executable statements ..
      x02ade = x02ake()/x02aje()
      return
      end
      real function x02aee(x)
c     mark 7 release. nag copyright 1978.
c     mark 11.5(f77) revised. (sept 1985.)
c     * minimum (largest negative) argument for exp *
c     returns the largest negative real number eneg such that
c     exp(eneg) can be successfully computed without underflow by
c     the compiler-supplied exp routine.
c     the value cannot be less than (i.e. more negative than)
c     alog(x02abe(x)), but may be greater because of specific
c     restrictions or limitations in the exp routine.
c     .. scalar arguments ..
      real                 x
c.(gbc)
cxx      data eneg /'ac50c3ae'x/
	data eneg/-87.33655/
c     (about  -87.3365)
c..
c     .. executable statements ..
      x02aee = eneg
      return
      end
      real function x02afe(x)
c     mark 9 release. nag copyright 1981.
c     mark 11.5(f77) revised. (sept 1985.)
c
c     * maximum argument for exp *
c     returns the largest positive real number epos such that
c     exp(epos) can be successfully computed without overflow
c     by the compiler supplied exp routine.
c
c     .. scalar arguments ..
      real                 x
c.(gbc)
cxx      data epos /'ac5043ae'x/
	data epos/87.33655/
c     (about 87.3365)
c..
c     .. executable statements ..
      x02afe = epos
      return
      end
      real function x02age(x)
c     mark 8 release. nag copyright 1980.
c     mark 11.5(f77) revised. (sept 1985.)
c
c     returns the smallest positive floating-point number  r
c     exactly representable on the computer such that  -r, 1.0/r,
c     and -1.0/r can all be computed without overflow or underflow.
c     on many machines the correct value can be derived from those
c     of x02aae, x02abe and x02ace as follows
c
c     if (x02abe(x)*x02ace(x).ge.1.0) x02age = x02abe(x)
c     if (x02abe(x)*x02ace(x).lt.1.0)
c    *                            x02age = (1.0+x02aae(x))/x02ace(x)
c
c     the correct value should be defined as a constant,
c     possibly in some binary, octal or hexadecimal representation,
c     and inserted into the assignment statement below.
c
c     x is a dummy argument
c
c     .. scalar arguments ..
      real                 x
c.(gbc)
cxx      data r /'10180'x/
	data r/1.1754945e-38/
c     (about 0.11754945e-37)
c..
c     .. executable statements ..
      x02age = r
      return
      end
      real function x02ahe(x)
c     mark 9 release. nag copyright 1981.
c     mark 11.5(f77) revised. (sept 1985.)
c
c     * maximum argument for sin and cos *
c     returns the largest positive real number maxsc such that
c     sin(maxsc) and cos(maxsc) can be successfully computed
c     by the compiler supplied sin and cos routines.
c
c     .. scalar arguments ..
      real                 x
c.(gbc)
      real maxsc
cxx      data maxsc /'4c00'x/
	data maxsc/0.8388608e7/
c     (about 0.83886e7)
c..
c     .. executable statements ..
      x02ahe = maxsc
      return
      end
      real function x02aje()
c     mark 12 release. nag copyright 1986.
c
c     returns  (1/2)*b**(1-p)  if rounds is .true.
c     returns  b**(1-p)  otherwise
c
c.(gbc)
cxx      data x /'3500'x/
c     (about 0.119209e-6)
	 data x /1.1920929e-07/
c     .. executable statements ..
      x02aje = x
c..
      return
      end
      real function x02ake()
c     mark 12 release. nag copyright 1986.
c
c     returns  b**(emin-1)  (the smallest positive model number)
c
c.(gbc)
cxxx      data x /'80'x/
cxxx	data x/1.1754945e-38/
cxxx	data x/2.9387359e-39/
c     (about 0.29387e-38)
      external x02ame
c     .. executable statements ..
      x02ake = x02ame()
c..
      return
      end
      real function x02ale()
c     mark 12 release. nag copyright 1986.
c
c     returns  (1 - b**(-p)) * b**emax  (the largest positive model
c     number)
c
c.(gbc)
cxx      data x /'ffff7f7f'x/
	data x/8.5070587e37/
c     (about 0.850706e38)
c     .. executable statements ..
      x02ale = x
c..
      return
      end
      real function x02ame()
c     mark 12 release. nag copyright 1986.
c
c     returns the 'safe range' parameter
c     i.e. the smallest positive model number z such that
c     for any x which satisfies x.ge.z and x.le.1/z
c     the following can be computed without overflow, underflow or other
c     error
c
c        -x
c        1.0/x
c        sqrt(x)
c        log(x)
c        exp(log(x))
c        y**(log(x)/log(y)) for any y
c
c.(gbc)
cxx      data x /'10180'x/
	data x/1.1754945e-38/
c     (about 0.11754945e-37)
c     .. executable statements ..
      x02ame = x
c..
      return
      end
      integer function x02bae(x)
c     mark 12 re-issue. nag copyright 1986.
c
c     returns the base of the arithmetic used on the computer
c     i.e. returns the same value as x02bhe
c
c     .. scalar arguments ..
      real                    x
c     .. external functions ..
      integer                 x02bhe
      external                x02bhe
c     .. executable statements ..
      x02bae = x02bhe()
      return
      end

      integer function x02bbe(x)
c     nag copyright 1975
c     mark 4.5 release
c     mark 11.5(f77) revised. (sept 1985.)
c     * maxint *
c     returns the largest integer representable on the computer
c     the x parameter is not used
c     .. scalar arguments ..
      real x
c.(gbc)
c      data maxint /'7fffffff'/
c     (2,147,483,647)
cc     .. executable statements ..
      x02bbe = 2147483647
cc..
      return
      end

      integer function x02bce(x)
c     nag copyright 1975
c     mark 4.5 release
c     mark 11.5(f77) revised. (sept 1985.)
c     * maxpw2 *
c     returns the largest integer to which power 2.0 may be raised
c     the x parameter is not used
c     .. scalar arguments ..
      real                    x
c     .. executable statements ..
c.(gbc)
      x02bce = 126
c..
      return
      end
      integer function x02bde(x)
c     nag copyright 1975
c     mark 4.5 release
c     mark 11.5(f77) revised. (sept 1985.)
c     * minpw2 *
c     returns the largest negative integer to which power 2.0 may
c     be raised
c     the x parameter is not used
c     .. scalar arguments ..
      real                    x
c     .. executable statements ..
c.(gbc)
      x02bde = -126
c..
      return
      end
      integer function x02bee(x)
c     nag copyright 1975
c     mark 4.5 release
c     mark 11.5(f77) revised. (sept 1985.)
c     * maxdec *
c     returns the maximum number of decimal digits that can be
c     accurately represented in the computer over the whole range
c     of floating-point numbers.
c     the x parameter is not used.
c     .. scalar arguments ..
      real                    x
c     .. executable statements ..
c.(gbc)
      x02bee = 7
c..
      return
      end
      integer function x02bhe()
c     mark 12 release. nag copyright 1986.
c
c     returns the model parameter, b.
c
c     .. executable statements ..
c.(gbc)
      x02bhe = 2
c..
      return
      end
      integer function x02bje()
c     mark 12 release. nag copyright 1986.
c
c     returns the model parameter, p.
c
c     .. executable statements ..
c.(gbc)
      x02bje = 24
c..
      return
      end
      integer function x02bke()
c     mark 12 release. nag copyright 1986.
c
c     returns the model parameter, emin.
c
c     .. executable statements ..
c.(gbc)
      x02bke = -127
c..
      return
      end
      integer function x02ble()
c     mark 12 release. nag copyright 1986.
c
c     returns the model parameter, emax.
c
c     .. executable statements ..
c.(gbc)
      x02ble = 126
c..
      return
      end
      integer function x02cae(x)
c     mark 7 release. nag copyright 1978.
c     mark 11.5(f77) revised. (sept 1985.)
c     * active set size *
c     in a non-paged environment, returns zero.
c     in a paged environment, returns a safe underestimate of the
c     amount of actual (as opposed to virtual) storage that is
c     expected to be available to a program running under typical
c     conditions (the expected active set size or working set size).
c     the amount of storage is given as the number of real
c     variables that can be stored in it. it is important that
c     the value should be an underestimate rather than an
c     overestimate. the value is used by linear algebra routines to
c     determine how many columns of a matrix can be handled together
c     as a block in the expectation that the whole block can fit
c     into actual storage.
c     the precise value is not critical. it does not affect the
c     numerical results, but merely the amount of paging performed
c     by some routines when solving large problems in a paged
c     environment. if in doubt, set the value to zero.
c     ***** the value may vary between installations *****
c     .. scalar arguments ..
      real                    x
c     .. executable statements ..
      x02cae = 0
      return
      end
      logical function x02dae(x)
c     mark 8 release. nag copyright 1980.
c     mark 11.5(f77) revised. (sept 1985.)
c
c     returns .false. if the system sets underflowing quantities
c     to zero, without any error indication or undesirable warning
c     or system overhead.
c     returns .true. otherwise, in which case certain library
c     routines will take special precautions to avoid underflow
c     (usually at some cost in efficiency).
c
c     x is a dummy argument
c
c     .. scalar arguments ..
      real                    x
c     .. executable statements ..
c.(gbc)
      x02dae = .false.
c..
      return
      end
      logical function x02dje()
c     mark 12 release. nag copyright 1986.
c
c     returns the model parameter, rounds.
c
c     .. executable statements ..
c.(gbc)
      x02dje = .false.
c..
      return
      end
      subroutine x04aae(i,nerr)
c     mark 7 release. nag copyright 1978
c     mark 7c revised ier-190 (may 1979)
c     mark 11.5(f77) revised. (sept 1985.)
c     if i = 0, sets nerr to current error message unit number
c     (stored in nerr1).
c     if i = 1, changes current error message unit number to
c     value specified by nerr.
c
c     .. scalar arguments ..
      integer           i, nerr
c     .. local scalars ..
      integer           nerr1
c     .. save statement ..
      save              nerr1
c     .. data statements ..
c.(gbc)
      data              nerr1/6/
c..
c     .. executable statements ..
      if (i.eq.0) nerr = nerr1
      if (i.eq.1) nerr1 = nerr
      return
      end
      subroutine x04bae(nout,rec)
c     mark 11.5(f77) release. nag copyright 1986.
c
c     x04bae writes the contents of rec to the unit defined by nout.
c
c     trailing blanks are not output, except that if rec is entirely
c     blank, a single blank character is output.
c     if nout.lt.0, i.e. if nout is not a valid fortran unit identifier,
c     then no output occurs.
c
c     .. scalar arguments ..
      integer           nout
      character*(*)     rec
c     .. local scalars ..
      integer           i
c     .. intrinsic functions ..
      intrinsic         len
c     .. executable statements ..
      if (nout.ge.0) then
c        remove trailing blanks
         do 20 i = len(rec), 2, -1
            if (rec(i:i).ne.' ') go to 40
   20    continue
c        write record to external file
   40    write (nout,fmt=99999) rec(1:i)
      end if
      return
c
99999 format (a)
      end


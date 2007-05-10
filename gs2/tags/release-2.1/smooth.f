      subroutine dsmooth (n, xin, yin, yout)
      implicit none
cin
      integer n
      real xin(n),yin(n)
cout
      real yout(n)
c
      integer max
      parameter (max=1000)
      real yp(max), ypp(max), xypp(max)
      real ypb, yppb, ypppb, xb
      integer maxiter
      parameter (maxiter=20)
      integer iter,i,j
      real diff
      parameter (diff=0.4)
c
      if (n .gt. max) then
         write (*,*) 'error in dsmooth: recompile with bigger max'
         stop
      endif
      do i=1,n
         yout(i)=yin(i)
      enddo
      do iter=1,maxiter
         do i=1,n-1
            yp(i) = (yout(i+1)-yout(i))/(xin(i+1)-xin(i))
         enddo
         do i=2,n-1
            xypp(i) = yout(i+1)-2.0*yout(i)+yout(i-1)
            ypp(i) = 2.0*(yp(i)-yp(i-1))/(xin(i+1)-xin(i-1))
         enddo
         ypb = 0.0
         yppb = 0.0
         xb = xin(n-2)-xin(2)
         do i=2,n-2
            ypb = ypb+abs(yp(i)/xb)
            yppb = yppb+abs(ypp(i)/xb)
         enddo
         if (ypb/yppb .gt. 2.0*xb/float(n)) goto 1
         do i=2,n-1
            yout(i) = yout(i) + diff*xypp(i)
         enddo
      enddo
 1    continue
      return
      end

      subroutine d2smooth (n, xin, yin, yout)
      implicit none
cin
      integer n
      real xin(n),yin(n)
cout
      real yout(n)
c
      integer max
      parameter (max=1000)
      real yp(max), ypp(max), yppp(max), xypp(max)
      real ypb, yppb, ypppb, xb
      integer maxiter
      parameter (maxiter=20)
      integer iter,i,j
      real diff
      parameter (diff=0.4)
c
      if (n .gt. max) then
         write (*,*) 'error in d2smooth: recompile with bigger max'
         stop
      endif
      do i=1,n
         yout(i)=yin(i)
      enddo
      do iter=1,maxiter
         do i=1,n-1
            yp(i) = (yout(i+1)-yout(i))/(xin(i+1)-xin(i))
         enddo
         do i=2,n-1
            xypp(i) = yout(i+1)-2.0*yout(i)+yout(i-1)
            ypp(i) = 2.0*(yp(i)-yp(i-1))/(xin(i+1)-xin(i-1))
         enddo
         do i=2,n-2
            yppp(i) = (ypp(i+1)-ypp(i))/(xin(i+1)-xin(i))
         enddo
         ypb = 0.0
         yppb = 0.0
         ypppb = 0.0
         xb = xin(n-2)-xin(2)
         do i=2,n-2
            ypb = ypb+abs(yp(i)/xb)
            yppb = yppb+abs(ypp(i)/xb)
            ypppb = ypppb+abs(yppp(i)/xb)
         enddo
         if (ypb/yppb .gt. 4.0*xb/float(n)
     >        .and. yppb/ypppb .gt. 4.0*xb/float(n))
     >             goto 1
         do i=2,n-1
            yout(i) = yout(i) + diff*xypp(i)
         enddo
      enddo
 1    continue
      return
      end

      subroutine smooth (n, xin, yin, var, yout, ifail)
      implicit none
cin
      integer n
      real xin(n),yin(n),var
cout
      real yout(n)
      integer ifail
c
      integer i
      integer max
      parameter (max=5000)
      double precision x(max),f(max),y(max),df(max),c(3*max),se(1)
      double precision wk(max+2,7)
      double precision dvar
c
      if (n .gt. max) then
         write (*,*) 'error in smooth: recompile with bigger max'
         stop
      endif
      ifail = 0
      do i=1,n
         x(i) = xin(i)
         f(i) = yin(i)
         df(i) = 1d0
      enddo
      dvar=var
      call cubgcv (x,f,df,n,y,c,n,dvar,0,se,wk,ifail)
      var=dvar
      do i=1,n
         yout(i) = y(i)
      enddo
      return
      end

C     ALGORITHM 642 COLLECTED ALGORITHMS FROM ACM.
C     ALGORITHM APPEARED IN ACM-TRANS. MATH. SOFTWARE, VOL.12, NO. 2,
C     JUN., 1986, P. 150.
C   SUBROUTINE NAME     - CUBGCV
C
C--------------------------------------------------------------------------
C
C   COMPUTER            - VAX/DOUBLE
C
C   AUTHOR              - M.F.HUTCHINSON
C                         CSIRO DIVISION OF MATHEMATICS AND STATISTICS
C                         P.O. BOX 1965
C                         CANBERRA, ACT 2601
C                         AUSTRALIA
C
C   LATEST REVISION     - 15 AUGUST 1985
C
C   PURPOSE             - CUBIC SPLINE DATA SMOOTHER
C
C   USAGE               - CALL CUBGCV (X,F,DF,N,Y,C,IC,VAR,JOB,SE,WK,IER)
C
C   ARGUMENTS    X      - VECTOR OF LENGTH N CONTAINING THE
C                           ABSCISSAE OF THE N DATA POINTS
C                           (X(I),F(I)) I=1..N. (INPUT) X
C                           MUST BE ORDERED SO THAT
C                           X(I) .LT. X(I+1).
C                F      - VECTOR OF LENGTH N CONTAINING THE
C                           ORDINATES (OR FUNCTION VALUES)
C                           OF THE N DATA POINTS (INPUT).
C                DF     - VECTOR OF LENGTH N. (INPUT/OUTPUT)
C                           DF(I) IS THE RELATIVE STANDARD DEVIATION
C                           OF THE ERROR ASSOCIATED WITH DATA POINT I.
C                           EACH DF(I) MUST BE POSITIVE.  THE VALUES IN
C                           DF ARE SCALED BY THE SUBROUTINE SO THAT
C                           THEIR MEAN SQUARE VALUE IS 1, AND UNSCALED
C                           AGAIN ON NORMAL EXIT.
C                           THE MEAN SQUARE VALUE OF THE DF(I) IS RETURNED
C                           IN WK(7) ON NORMAL EXIT.
C                           IF THE ABSOLUTE STANDARD DEVIATIONS ARE KNOWN,
C                           THESE SHOULD BE PROVIDED IN DF AND THE ERROR
C                           VARIANCE PARAMETER VAR (SEE BELOW) SHOULD THEN
C                           BE SET TO 1.
C                           IF THE RELATIVE STANDARD DEVIATIONS ARE UNKNOWN,
C                           SET EACH DF(I)=1.
C                N      - NUMBER OF DATA POINTS (INPUT).
C                           N MUST BE .GE. 3.
C                Y,C    - SPLINE COEFFICIENTS. (OUTPUT) Y
C                           IS A VECTOR OF LENGTH N. C IS
C                           AN N-1 BY 3 MATRIX. THE VALUE
C                           OF THE SPLINE APPROXIMATION AT T IS
C                           S(T)=((C(I,3)*D+C(I,2))*D+C(I,1))*D+Y(I)
C                           WHERE X(I).LE.T.LT.X(I+1) AND
C                           D = T-X(I).
C                IC     - ROW DIMENSION OF MATRIX C EXACTLY
C                           AS SPECIFIED IN THE DIMENSION
C                           STATEMENT IN THE CALLING PROGRAM. (INPUT)
C                VAR    - ERROR VARIANCE. (INPUT/OUTPUT)
C                           IF VAR IS NEGATIVE (I.E. UNKNOWN) THEN
C                           THE SMOOTHING PARAMETER IS DETERMINED
C                           BY MINIMIZING THE GENERALIZED CROSS VALIDATION
C                           AND AN ESTIMATE OF THE ERROR VARIANCE IS
C                           RETURNED IN VAR.
C                           IF VAR IS NON-NEGATIVE (I.E. KNOWN) THEN THE
C                           SMOOTHING PARAMETER IS DETERMINED TO MINIMIZE
C                           AN ESTIMATE, WHICH DEPENDS ON VAR, OF THE TRUE
C                           MEAN SQUARE ERROR, AND VAR IS UNCHANGED.
C                           IN PARTICULAR, IF VAR IS ZERO, THEN AN
C                           INTERPOLATING NATURAL CUBIC SPLINE IS CALCULATED.
C                           VAR SHOULD BE SET TO 1 IF ABSOLUTE STANDARD
C                           DEVIATIONS HAVE BEEN PROVIDED IN DF (SEE ABOVE).
C                JOB    - JOB SELECTION PARAMETER. (INPUT)
C                         JOB = 0 SHOULD BE SELECTED IF POINT STANDARD ERROR
C                           ESTIMATES ARE NOT REQUIRED IN SE.
C                         JOB = 1 SHOULD BE SELECTED IF POINT STANDARD ERROR
C                           ESTIMATES ARE REQUIRED IN SE.
C                SE     - VECTOR OF LENGTH N CONTAINING BAYESIAN STANDARD
C                           ERROR ESTIMATES OF THE FITTED SPLINE VALUES IN Y.
C                           SE IS NOT REFERENCED IF JOB=0. (OUTPUT)
C                WK     - WORK VECTOR OF LENGTH 7*(N + 2). ON NORMAL EXIT THE
C                           FIRST 7 VALUES OF WK ARE ASSIGNED AS FOLLOWS:-
C
C                           WK(1) = SMOOTHING PARAMETER (= RHO/(RHO + 1))
C                           WK(2) = ESTIMATE OF THE NUMBER OF DEGREES OF
C                                   FREEDOM OF THE RESIDUAL SUM OF SQUARES
C                           WK(3) = GENERALIZED CROSS VALIDATION
C                           WK(4) = MEAN SQUARE RESIDUAL
C                           WK(5) = ESTIMATE OF THE TRUE MEAN SQUARE ERROR
C                                   AT THE DATA POINTS
C                           WK(6) = ESTIMATE OF THE ERROR VARIANCE
C                           WK(7) = MEAN SQUARE VALUE OF THE DF(I)
C
C                           IF WK(1)=0 (RHO=0) AN INTERPOLATING NATURAL CUBIC
C                           SPLINE HAS BEEN CALCULATED.
C                           IF WK(1)=1 (RHO=INFINITE) A LEAST SQUARES
C                           REGRESSION LINE HAS BEEN CALCULATED.
C                           WK(2) IS AN ESTIMATE OF THE NUMBER OF DEGREES OF
C                           FREEDOM OF THE RESIDUAL WHICH REDUCES TO THE
C                           USUAL VALUE OF N-2 WHEN A LEAST SQUARES REGRESSION
C                           LINE IS CALCULATED.
C                           WK(3),WK(4),WK(5) ARE CALCULATED WITH THE DF(I)
C                           SCALED TO HAVE MEAN SQUARE VALUE 1.  THE
C                           UNSCALED VALUES OF WK(3),WK(4),WK(5) MAY BE
C                           CALCULATED BY DIVIDING BY WK(7).
C                           WK(6) COINCIDES WITH THE OUTPUT VALUE OF VAR IF
C                           VAR IS NEGATIVE ON INPUT.  IT IS CALCULATED WITH
C                           THE UNSCALED VALUES OF THE DF(I) TO FACILITATE
C                           COMPARISONS WITH A PRIORI VARIANCE ESTIMATES.
C
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129, IC IS LESS THAN N-1.
C                           IER = 130, N IS LESS THAN 3.
C                           IER = 131, INPUT ABSCISSAE ARE NOT
C                             ORDERED SO THAT X(I).LT.X(I+1).
C                           IER = 132, DF(I) IS NOT POSITIVE FOR SOME I.
C                           IER = 133, JOB IS NOT 0 OR 1.
C
C   PRECISION/HARDWARE  - DOUBLE
C
C   REQUIRED ROUTINES   - SPINT1,SPFIT1,SPCOF1,SPERR1
C
C   REMARKS      THE NUMBER OF ARITHMETIC OPERATIONS REQUIRED BY THE
C                SUBROUTINE IS PROPORTIONAL TO N.  THE SUBROUTINE
C                USES AN ALGORITHM DEVELOPED BY M.F. HUTCHINSON AND
C                F.R. DE HOOG, 'SMOOTHING NOISY DATA WITH SPLINE
C                FUNCTIONS', NUMER. MATH. (IN PRESS)
C
C-----------------------------------------------------------------------
C
      SUBROUTINE CUBGCV(X,F,DF,N,Y,C,IC,VAR,JOB,SE,WK,IER)
C
C---SPECIFICATIONS FOR ARGUMENTS---
      INTEGER N,IC,JOB,IER
      DOUBLE PRECISION X(N),F(N),DF(N),Y(N),C(IC,3),SE(N),VAR,
     .                 WK(0:N+1,7)
C
C---SPECIFICATIONS FOR LOCAL VARIABLES---
      DOUBLE PRECISION DELTA,ERR,GF1,GF2,GF3,GF4,R1,R2,R3,R4,TAU,RATIO,
     .                 AVH,AVDF,AVAR,ZERO,ONE,STAT(6),P,Q
C
      DATA RATIO/2.0D0/
      DATA TAU/1.618033989D0/
      DATA ZERO,ONE/0.0D0,1.0D0/
C
C---INITIALIZE---
      IER = 133
      IF (JOB.LT.0 .OR. JOB.GT.1) GO TO 140
      CALL SPINT1(X,AVH,F,DF,AVDF,N,Y,C,IC,WK,WK(0,4),IER)
      IF (IER.NE.0) GO TO 140
      AVAR = VAR
      IF (VAR.GT.ZERO) AVAR = VAR*AVDF*AVDF
C
C---CHECK FOR ZERO VARIANCE---
      IF (VAR.NE.ZERO) GO TO 10
      R1 = ZERO
      GO TO 90
C
C---FIND LOCAL MINIMUM OF GCV OR THE EXPECTED MEAN SQUARE ERROR---
   10 R1 = ONE
      R2 = RATIO*R1
      CALL SPFIT1(X,AVH,DF,N,R2,P,Q,GF2,AVAR,STAT,Y,C,IC,WK,WK(0,4),
     .            WK(0,6),WK(0,7))
   20 CALL SPFIT1(X,AVH,DF,N,R1,P,Q,GF1,AVAR,STAT,Y,C,IC,WK,WK(0,4),
     .            WK(0,6),WK(0,7))
      IF (GF1.GT.GF2) GO TO 30
C
C---EXIT IF P ZERO---
      IF (P.LE.ZERO) GO TO 100
      R2 = R1
      GF2 = GF1
      R1 = R1/RATIO
      GO TO 20

   30 R3 = RATIO*R2
   40 CALL SPFIT1(X,AVH,DF,N,R3,P,Q,GF3,AVAR,STAT,Y,C,IC,WK,WK(0,4),
     .            WK(0,6),WK(0,7))
      IF (GF3.GT.GF2) GO TO 50
C
C---EXIT IF Q ZERO---
      IF (Q.LE.ZERO) GO TO 100
      R2 = R3
      GF2 = GF3
      R3 = RATIO*R3
      GO TO 40

   50 R2 = R3
      GF2 = GF3
      DELTA = (R2-R1)/TAU
      R4 = R1 + DELTA
      R3 = R2 - DELTA
      CALL SPFIT1(X,AVH,DF,N,R3,P,Q,GF3,AVAR,STAT,Y,C,IC,WK,WK(0,4),
     .            WK(0,6),WK(0,7))
      CALL SPFIT1(X,AVH,DF,N,R4,P,Q,GF4,AVAR,STAT,Y,C,IC,WK,WK(0,4),
     .            WK(0,6),WK(0,7))
C
C---GOLDEN SECTION SEARCH FOR LOCAL MINIMUM---
   60 IF (GF3.GT.GF4) GO TO 70
      R2 = R4
      GF2 = GF4
      R4 = R3
      GF4 = GF3
      DELTA = DELTA/TAU
      R3 = R2 - DELTA
      CALL SPFIT1(X,AVH,DF,N,R3,P,Q,GF3,AVAR,STAT,Y,C,IC,WK,WK(0,4),
     .            WK(0,6),WK(0,7))
      GO TO 80

   70 R1 = R3
      GF1 = GF3
      R3 = R4
      GF3 = GF4
      DELTA = DELTA/TAU
      R4 = R1 + DELTA
      CALL SPFIT1(X,AVH,DF,N,R4,P,Q,GF4,AVAR,STAT,Y,C,IC,WK,WK(0,4),
     .            WK(0,6),WK(0,7))
   80 ERR = (R2-R1)/ (R1+R2)
      IF (ERR*ERR+ONE.GT.ONE .AND. ERR.GT.1.0D-6) GO TO 60
      R1 = (R1+R2)*0.5D0
C
C---CALCULATE SPLINE COEFFICIENTS---
   90 CALL SPFIT1(X,AVH,DF,N,R1,P,Q,GF1,AVAR,STAT,Y,C,IC,WK,WK(0,4),
     .            WK(0,6),WK(0,7))
  100 CALL SPCOF1(X,AVH,F,DF,N,P,Q,Y,C,IC,WK(0,6),WK(0,7))
C
C---OPTIONALLY CALCULATE STANDARD ERROR ESTIMATES---
      IF (VAR.GE.ZERO) GO TO 110
      AVAR = STAT(6)
      VAR = AVAR/ (AVDF*AVDF)
  110 IF (JOB.EQ.1) CALL SPERR1(X,AVH,DF,N,WK,P,AVAR,SE)
C
C---UNSCALE DF---
      DO 120 I = 1,N
         DF(I) = DF(I)*AVDF
  120 CONTINUE
C
C--PUT STATISTICS IN WK---
      DO 130 I = 0,5
         WK(I,1) = STAT(I+1)
  130 CONTINUE
      WK(5,1) = STAT(6)/ (AVDF*AVDF)
      WK(6,1) = AVDF*AVDF
      GO TO 150
C
C---CHECK FOR ERROR CONDITION---
  140 CONTINUE
C     IF (IER.NE.0) CONTINUE
  150 RETURN
      END
      SUBROUTINE SPINT1(X,AVH,Y,DY,AVDY,N,A,C,IC,R,T,IER)
C
C INITIALIZES THE ARRAYS C, R AND T FOR ONE DIMENSIONAL CUBIC
C SMOOTHING SPLINE FITTING BY SUBROUTINE SPFIT1.  THE VALUES
C DF(I) ARE SCALED SO THAT THE SUM OF THEIR SQUARES IS N
C AND THE AVERAGE OF THE DIFFERENCES X(I+1) - X(I) IS CALCULATED
C IN AVH IN ORDER TO AVOID UNDERFLOW AND OVERFLOW PROBLEMS IN
C SPFIT1.
C
C SUBROUTINE SETS IER IF ELEMENTS OF X ARE NON-INCREASING,
C IF N IS LESS THAN 3, IF IC IS LESS THAN N-1 OR IF DY(I) IS
C NOT POSITIVE FOR SOME I.
C
C---SPECIFICATIONS FOR ARGUMENTS---
      INTEGER N,IC,IER
      DOUBLE PRECISION X(N),Y(N),DY(N),A(N),C(IC,3),R(0:N+1,3),
     .                 T(0:N+1,2),AVH,AVDY
C
C---SPECIFICATIONS FOR LOCAL VARIABLES---
      INTEGER I
      DOUBLE PRECISION E,F,G,H,ZERO
      DATA ZERO/0.0D0/
C
C---INITIALIZATION AND INPUT CHECKING---
      IER = 0
      IF (N.LT.3) GO TO 60
      IF (IC.LT.N-1) GO TO 70
C
C---GET AVERAGE X SPACING IN AVH---
      G = ZERO
      DO 10 I = 1,N - 1
         H = X(I+1) - X(I)
         IF (H.LE.ZERO) GO TO 80
         G = G + H
   10 CONTINUE
      AVH = G/ (N-1)
C
C---SCALE RELATIVE WEIGHTS---
      G = ZERO
      DO 20 I = 1,N
         IF (DY(I).LE.ZERO) GO TO 90
         G = G + DY(I)*DY(I)
   20 CONTINUE
      AVDY = DSQRT(G/N)
C
      DO 30 I = 1,N
         DY(I) = DY(I)/AVDY
   30 CONTINUE
C
C---INITIALIZE H,F---
      H = (X(2)-X(1))/AVH
      F = (Y(2)-Y(1))/H
C
C---CALCULATE A,T,R---
      DO 40 I = 2,N - 1
         G = H
         H = (X(I+1)-X(I))/AVH
         E = F
         F = (Y(I+1)-Y(I))/H
         A(I) = F - E
         T(I,1) = 2.0D0* (G+H)/3.0D0
         T(I,2) = H/3.0D0
         R(I,3) = DY(I-1)/G
         R(I,1) = DY(I+1)/H
         R(I,2) = -DY(I)/G - DY(I)/H
   40 CONTINUE
C
C---CALCULATE C = R'*R---
      R(N,2) = ZERO
      R(N,3) = ZERO
      R(N+1,3) = ZERO
      DO 50 I = 2,N - 1
         C(I,1) = R(I,1)*R(I,1) + R(I,2)*R(I,2) + R(I,3)*R(I,3)
         C(I,2) = R(I,1)*R(I+1,2) + R(I,2)*R(I+1,3)
         C(I,3) = R(I,1)*R(I+2,3)
   50 CONTINUE
      RETURN
C
C---ERROR CONDITIONS---
   60 IER = 130
      RETURN

   70 IER = 129
      RETURN

   80 IER = 131
      RETURN

   90 IER = 132
      RETURN
      END
      SUBROUTINE SPFIT1(X,AVH,DY,N,RHO,P,Q,FUN,VAR,STAT,A,C,IC,R,T,U,V)
C
C FITS A CUBIC SMOOTHING SPLINE TO DATA WITH RELATIVE
C WEIGHTING DY FOR A GIVEN VALUE OF THE SMOOTHING PARAMETER
C RHO USING AN ALGORITHM BASED ON THAT OF C.H. REINSCH (1967),
C NUMER. MATH. 10, 177-183.
C
C THE TRACE OF THE INFLUENCE MATRIX IS CALCULATED USING AN
C ALGORITHM DEVELOPED BY M.F.HUTCHINSON AND F.R.DE HOOG (NUMER.
C MATH., IN PRESS), ENABLING THE GENERALIZED CROSS VALIDATION
C AND RELATED STATISTICS TO BE CALCULATED IN ORDER N OPERATIONS.
C
C THE ARRAYS A, C, R AND T ARE ASSUMED TO HAVE BEEN INITIALIZED
C BY THE SUBROUTINE SPINT1.  OVERFLOW AND UNDERFLOW PROBLEMS ARE
C AVOIDED BY USING P=RHO/(1 + RHO) AND Q=1/(1 + RHO) INSTEAD OF
C RHO AND BY SCALING THE DIFFERENCES X(I+1) - X(I) BY AVH.
C
C THE VALUES IN DF ARE ASSUMED TO HAVE BEEN SCALED SO THAT THE
C SUM OF THEIR SQUARED VALUES IS N.  THE VALUE IN VAR, WHEN IT IS
C NON-NEGATIVE, IS ASSUMED TO HAVE BEEN SCALED TO COMPENSATE FOR
C THE SCALING OF THE VALUES IN DF.
C
C THE VALUE RETURNED IN FUN IS AN ESTIMATE OF THE TRUE MEAN SQUARE
C WHEN VAR IS NON-NEGATIVE, AND IS THE GENERALIZED CROSS VALIDATION
C WHEN VAR IS NEGATIVE.
C
C---SPECIFICATIONS FOR ARGUMENTS---
      INTEGER IC,N
      DOUBLE PRECISION X(N),DY(N),RHO,STAT(6),A(N),C(IC,3),R(0:N+1,3),
     .                 T(0:N+1,2),U(0:N+1),V(0:N+1),FUN,VAR,AVH,P,Q
C
C---LOCAL VARIABLES---
      INTEGER I
      DOUBLE PRECISION E,F,G,H,ZERO,ONE,TWO,RHO1
      DATA ZERO,ONE,TWO/0.0D0,1.0D0,2.0D0/
C
C---USE P AND Q INSTEAD OF RHO TO PREVENT OVERFLOW OR UNDERFLOW---
      RHO1 = ONE + RHO
      P = RHO/RHO1
      Q = ONE/RHO1
      IF (RHO1.EQ.ONE) P = ZERO
      IF (RHO1.EQ.RHO) Q = ZERO
C
C---RATIONAL CHOLESKY DECOMPOSITION OF P*C + Q*T---
      F = ZERO
      G = ZERO
      H = ZERO
      DO 10 I = 0,1
         R(I,1) = ZERO
   10 CONTINUE
      DO 20 I = 2,N - 1
         R(I-2,3) = G*R(I-2,1)
         R(I-1,2) = F*R(I-1,1)
         R(I,1) = ONE/ (P*C(I,1)+Q*T(I,1)-F*R(I-1,2)-G*R(I-2,3))
         F = P*C(I,2) + Q*T(I,2) - H*R(I-1,2)
         G = H
         H = P*C(I,3)
   20 CONTINUE
C
C---SOLVE FOR U---
      U(0) = ZERO
      U(1) = ZERO
      DO 30 I = 2,N - 1
         U(I) = A(I) - R(I-1,2)*U(I-1) - R(I-2,3)*U(I-2)
   30 CONTINUE
      U(N) = ZERO
      U(N+1) = ZERO
      DO 40 I = N - 1,2,-1
         U(I) = R(I,1)*U(I) - R(I,2)*U(I+1) - R(I,3)*U(I+2)
   40 CONTINUE
C
C---CALCULATE RESIDUAL VECTOR V---
      E = ZERO
      H = ZERO
      DO 50 I = 1,N - 1
         G = H
         H = (U(I+1)-U(I))/ ((X(I+1)-X(I))/AVH)
         V(I) = DY(I)* (H-G)
         E = E + V(I)*V(I)
   50 CONTINUE
      V(N) = DY(N)* (-H)
      E = E + V(N)*V(N)
C
C---CALCULATE UPPER THREE BANDS OF INVERSE MATRIX---
      R(N,1) = ZERO
      R(N,2) = ZERO
      R(N+1,1) = ZERO
      DO 60 I = N - 1,2,-1
         G = R(I,2)
         H = R(I,3)
         R(I,2) = -G*R(I+1,1) - H*R(I+1,2)
         R(I,3) = -G*R(I+1,2) - H*R(I+2,1)
         R(I,1) = R(I,1) - G*R(I,2) - H*R(I,3)
   60 CONTINUE
C
C---CALCULATE TRACE---
      F = ZERO
      G = ZERO
      H = ZERO
      DO 70 I = 2,N - 1
         F = F + R(I,1)*C(I,1)
         G = G + R(I,2)*C(I,2)
         H = H + R(I,3)*C(I,3)
   70 CONTINUE
      F = F + TWO* (G+H)
C
C---CALCULATE STATISTICS---
      STAT(1) = P
      STAT(2) = F*P
      STAT(3) = N*E/ (F*F)
      STAT(4) = E*P*P/N
      STAT(6) = E*P/F
      IF (VAR.GE.ZERO) GO TO 80
      STAT(5) = STAT(6) - STAT(4)
      FUN = STAT(3)
      GO TO 90

   80 STAT(5) = DMAX1(STAT(4)-TWO*VAR*STAT(2)/N+VAR,ZERO)
      FUN = STAT(5)
   90 RETURN
      END
      SUBROUTINE SPERR1(X,AVH,DY,N,R,P,VAR,SE)
C
C CALCULATES BAYESIAN ESTIMATES OF THE STANDARD ERRORS OF THE FITTED
C VALUES OF A CUBIC SMOOTHING SPLINE BY CALCULATING THE DIAGONAL ELEMENTS
C OF THE INFLUENCE MATRIX.
C
C---SPECIFICATIONS FOR ARGUMENTS---
      INTEGER N
      DOUBLE PRECISION X(N),DY(N),R(0:N+1,3),SE(N),AVH,P,VAR
C
C---SPECIFICATIONS FOR LOCAL VARIABLES---
      INTEGER I
      DOUBLE PRECISION F,G,H,F1,G1,H1,ZERO,ONE
      DATA ZERO,ONE/0.0D0,1.0D0/
C
C---INITIALIZE---
      H = AVH/ (X(2)-X(1))
      SE(1) = ONE - P*DY(1)*DY(1)*H*H*R(2,1)
      R(1,1) = ZERO
      R(1,2) = ZERO
      R(1,3) = ZERO
C
C---CALCULATE DIAGONAL ELEMENTS---
      DO 10 I = 2,N - 1
         F = H
         H = AVH/ (X(I+1)-X(I))
         G = -F - H
         F1 = F*R(I-1,1) + G*R(I-1,2) + H*R(I-1,3)
         G1 = F*R(I-1,2) + G*R(I,1) + H*R(I,2)
         H1 = F*R(I-1,3) + G*R(I,2) + H*R(I+1,1)
         SE(I) = ONE - P*DY(I)*DY(I)* (F*F1+G*G1+H*H1)
   10 CONTINUE
      SE(N) = ONE - P*DY(N)*DY(N)*H*H*R(N-1,1)
C
C---CALCULATE STANDARD ERROR ESTIMATES---
      DO 20 I = 1,N
         SE(I) = DSQRT(DMAX1(SE(I)*VAR,ZERO))*DY(I)
   20 CONTINUE
      RETURN
      END
      SUBROUTINE SPCOF1(X,AVH,Y,DY,N,P,Q,A,C,IC,U,V)
C
C CALCULATES COEFFICIENTS OF A CUBIC SMOOTHING SPLINE FROM
C PARAMETERS CALCULATED BY SUBROUTINE SPFIT1.
C
C---SPECIFICATIONS FOR ARGUMENTS---
      INTEGER IC,N
      DOUBLE PRECISION X(N),Y(N),DY(N),P,Q,A(N),C(IC,3),U(0:N+1),
     .                 V(0:N+1),AVH
C
C---SPECIFICATIONS FOR LOCAL VARIABLES---
      INTEGER I
      DOUBLE PRECISION H,QH
C
C---CALCULATE A---
      QH = Q/ (AVH*AVH)
      DO 10 I = 1,N
         A(I) = Y(I) - P*DY(I)*V(I)
         U(I) = QH*U(I)
   10 CONTINUE
C
C---CALCULATE C---
      DO 20 I = 1,N - 1
         H = X(I+1) - X(I)
         C(I,3) = (U(I+1)-U(I))/ (3.0D0*H)
         C(I,1) = (A(I+1)-A(I))/H - (H*C(I,3)+U(I))*H
         C(I,2) = U(I)
   20 CONTINUE
      RETURN
      END

C---C CUBGCV TEST DRIVER
C---C ------------------
C---C
C---C AUTHOR          - M.F.HUTCHINSON
C---C                   CSIRO DIVISION OF WATER AND LAND RESOURCES
C---C                   GPO BOX 1666
C---C                   CANBERRA ACT 2601
C---C                   AUSTRALIA
C---C
C---C LATEST REVISION - 7 AUGUST 1986
C---C
C---C COMPUTER        - VAX/DOUBLE
C---C
C---C USAGE           - MAIN PROGRAM
C---C
C---C REQUIRED ROUTINES - CUBGCV,SPINT1,SPFIT1,SPCOF1,SPERR1,GGRAND
C---C
C---C REMARKS   USES SUBROUTINE CUBGCV TO FIT A CUBIC SMOOTHING SPLINE
C---C           TO 50 DATA POINTS WHICH ARE GENERATED BY ADDING A RANDOM
C---C           VARIABLE WITH UNIFORM DENSITY IN THE INTERVAL [-0.3,0.3]
C---C           TO 50 POINTS SAMPLED FROM THE CURVE  Y=SIN(3*PI*X/2).
C---C           RANDOM DEVIATES IN THE INTERVAL [0,1] ARE GENERATED BY THE
C---C           DOUBLE PRECISION FUNCTION GGRAND (SIMILAR TO IMSL FUNCTION
C---C           GGUBFS).  THE ABSCISSAE ARE UNEQUALLY SPACED IN [0,1].
C---C
C---C           POINT STANDARD ERROR ESTIMATES ARE RETURNED IN SE BY
C---C           SETTING JOB=1.  THE ERROR VARIANCE ESTIMATE IS RETURNED
C---C           IN VAR.  IT COMPARES FAVOURABLY WITH THE TRUE VALUE OF 0.03.
C---C           SUMMARY STATISTICS FROM THE ARRAY WK ARE WRITTEN TO
C---C           UNIT 6.  DATA VALUES AND FITTED VALUES WITH ESTIMATED
C---C           STANDARD ERRORS ARE ALSO WRITTEN TO UNIT 6.
C---C
C---      PARAMETER (N=50, IC=49)
C---C
C---      INTEGER            JOB,IER
C---      DOUBLE PRECISION   X(N),F(N),Y(N),DF(N),C(IC,3),WK(7*(N+2)),
C---     *                   VAR,SE(N)
C---      DOUBLE PRECISION   GGRAND,DSEED
C---C
C---C---INITIALIZE---
C---      DSEED=1.2345D4
C---      JOB=1
C---      VAR=-1.0D0
C---C
C---C---CALCULATE DATA POINTS---
C---      DO 10 I=1,N
C---      X(I)=(I - 0.5)/N + (2.0*GGRAND(DSEED) - 1.0)/(3.0*N)
C---      F(I)=DSIN(4.71238*X(I)) + (2.0*GGRAND(DSEED) - 1.0)*0.3
C---      DF(I)=1.0D0
C---  10  CONTINUE
C---C
C---C---FIT CUBIC SPLINE---
C---      CALL CUBGCV(X,F,DF,N,Y,C,IC,VAR,JOB,SE,WK,IER)
C---C
C---C---WRITE OUT RESULTS---
C---      WRITE(6,20)
C---  20  FORMAT(' CUBGCV TEST DRIVER RESULTS:')
C---      WRITE(6,30) IER,VAR,WK(3),WK(4),WK(2)
C---  30  FORMAT(/' IER =',I4/' VAR =',F7.4/
C---     *        ' GENERALIZED CROSS VALIDATION =',F7.4/
C---     *        ' MEAN SQUARE RESIDUAL         =',F7.4/
C---     *        ' RESIDUAL DEGREES OF FREEDOM  =',F7.2)
C---      WRITE(6,40)
C---  40  FORMAT(/' INPUT DATA',17X,'OUTPUT RESULTS'//
C---     *         '   I    X(I)    F(I)',6X,'    Y(I)   SE(I)',
C---     *          '      C(I,1)      C(I,2)      C(I,3)')
C---      DO 60 I=1,N-1
C---      WRITE(6,50) I,X(I),F(I),Y(I),SE(I),(C(I,J),J=1,3)
C---  50  FORMAT(I4,2F8.4,6X,2F8.4,3E12.4)
C---  60  CONTINUE
C---      WRITE(6,50) N,X(N),F(N),Y(N),SE(N)
C---      STOP
C---      END
C---      DOUBLE PRECISION FUNCTION GGRAND(DSEED)
C---C
C---C DOUBLE PRECISION UNIFORM RANDOM NUMBER GENERATOR
C---C
C---C CONSTANTS: A = 7**5
C---C            B = 2**31 - 1
C---C            C = 2**31
C---C
C---C REFERENCE: IMSL MANUAL, CHAPTER G - GENERATION AND TESTING OF
C---C                                     RANDOM NUMBERS
C---C
C---C---SPECIFICATIONS FOR ARGUMENTS---
C---      DOUBLE PRECISION DSEED
C---C
C---C---SPECIFICATIONS FOR LOCAL VARIABLES---
C---      DOUBLE PRECISION A,B,C,S
C---C
C---      DATA A,B,C/16807.0D0, 2147483647.0D0, 2147483648.0D0/
C---C
C---      S=DSEED
C---      S=DMOD(A*S, B)
C---      GGRAND=S/C
C---      DSEED=S
C---      RETURN
C---      END
C---
C---CUBGCV TEST DRIVER RESULTS:
C---
C---IER =   0
C---VAR = 0.0279
C---GENERALIZED CROSS VALIDATION = 0.0318
C---MEAN SQUARE RESIDUAL         = 0.0246
C---RESIDUAL DEGREES OF FREEDOM  =  43.97
C---
C---INPUT DATA                 OUTPUT RESULTS
C---
C---  I    X(I)    F(I)          Y(I)   SE(I)      C(I,1)      C(I,2)      C(I,3)
C---  1  0.0046  0.2222        0.0342  0.1004  0.3630E+01  0.0000E+00  0.2542E+02
C---  2  0.0360 -0.1098        0.1488  0.0750  0.3705E+01  0.2391E+01 -0.9537E+01
C---  3  0.0435 -0.0658        0.1767  0.0707  0.3740E+01  0.2175E+01 -0.4233E+02
C---  4  0.0735  0.3906        0.2900  0.0594  0.3756E+01 -0.1642E+01 -0.2872E+02
C---  5  0.0955  0.6054        0.3714  0.0558  0.3642E+01 -0.3535E+01  0.2911E+01
C---  6  0.1078  0.3034        0.4155  0.0549  0.3557E+01 -0.3428E+01 -0.1225E+02
C---  7  0.1269  0.7386        0.4822  0.0544  0.3412E+01 -0.4131E+01  0.2242E+02
C---  8  0.1565  0.4616        0.5800  0.0543  0.3227E+01 -0.2143E+01  0.6415E+01
C---  9  0.1679  0.4315        0.6165  0.0543  0.3180E+01 -0.1923E+01 -0.1860E+02
C--- 10  0.1869  0.5716        0.6762  0.0544  0.3087E+01 -0.2985E+01 -0.3274E+02
C--- 11  0.2149  0.6736        0.7595  0.0542  0.2843E+01 -0.5733E+01 -0.4435E+02
C--- 12  0.2356  0.7388        0.8155  0.0539  0.2549E+01 -0.8486E+01 -0.5472E+02
C--- 13  0.2557  1.1953        0.8630  0.0537  0.2139E+01 -0.1180E+02 -0.9784E+01
C--- 14  0.2674  1.0299        0.8864  0.0536  0.1860E+01 -0.1214E+02  0.9619E+01
C--- 15  0.2902  0.7981        0.9225  0.0534  0.1322E+01 -0.1149E+02 -0.7202E+01
C--- 16  0.3155  0.8973        0.9485  0.0532  0.7269E+00 -0.1203E+02 -0.1412E+02
C--- 17  0.3364  1.2695        0.9583  0.0530  0.2040E+00 -0.1292E+02  0.2796E+02
C--- 18  0.3557  0.7253        0.9577  0.0527 -0.2638E+00 -0.1130E+02 -0.3453E+01
C--- 19  0.3756  1.2127        0.9479  0.0526 -0.7176E+00 -0.1151E+02  0.3235E+02
C--- 20  0.3881  0.7304        0.9373  0.0525 -0.9889E+00 -0.1030E+02  0.4381E+01
C--- 21  0.4126  0.9810        0.9069  0.0525 -0.1486E+01 -0.9977E+01  0.1440E+02
C--- 22  0.4266  0.7117        0.8842  0.0525 -0.1756E+01 -0.9373E+01 -0.8925E+01
C--- 23  0.4566  0.7203        0.8227  0.0524 -0.2344E+01 -0.1018E+02 -0.2278E+02
C--- 24  0.4704  0.9242        0.7884  0.0524 -0.2637E+01 -0.1112E+02 -0.4419E+01
C--- 25  0.4914  0.7345        0.7281  0.0523 -0.3110E+01 -0.1140E+02 -0.3562E+01
C--- 26  0.5084  0.7378        0.6720  0.0524 -0.3500E+01 -0.1158E+02  0.5336E+01
C--- 27  0.5277  0.7441        0.6002  0.0525 -0.3941E+01 -0.1127E+02  0.2479E+02
C--- 28  0.5450  0.5612        0.5286  0.0527 -0.4310E+01 -0.9980E+01  0.2920E+02
C--- 29  0.5641  0.5049        0.4429  0.0529 -0.4659E+01 -0.8309E+01  0.3758E+02
C--- 30  0.5857  0.4725        0.3390  0.0531 -0.4964E+01 -0.5878E+01  0.5563E+02
C--- 31  0.6159  0.1380        0.1850  0.0531 -0.5167E+01 -0.8307E+00  0.4928E+02
C--- 32  0.6317  0.1412        0.1036  0.0531 -0.5157E+01  0.1499E+01  0.5437E+02
C--- 33  0.6446 -0.1110        0.0371  0.0531 -0.5091E+01  0.3614E+01  0.3434E+02
C--- 34  0.6707 -0.2605       -0.0927  0.0532 -0.4832E+01  0.6302E+01  0.1164E+02
C--- 35  0.6853 -0.1284       -0.1619  0.0533 -0.4640E+01  0.6812E+01  0.1617E+02
C--- 36  0.7064 -0.3452       -0.2564  0.0536 -0.4332E+01  0.7834E+01  0.4164E+01
C--- 37  0.7310 -0.5527       -0.3582  0.0538 -0.3939E+01  0.8141E+01 -0.2214E+02
C--- 38  0.7531 -0.3459       -0.4415  0.0540 -0.3611E+01  0.6674E+01 -0.9205E+01
C--- 39  0.7686 -0.5902       -0.4961  0.0541 -0.3410E+01  0.6245E+01 -0.2193E+02
C--- 40  0.7952 -0.7644       -0.5828  0.0541 -0.3125E+01  0.4494E+01 -0.4649E+02
C--- 41  0.8087 -0.5392       -0.6242  0.0541 -0.3029E+01  0.2614E+01 -0.3499E+02
C--- 42  0.8352 -0.4247       -0.7031  0.0539 -0.2964E+01 -0.1603E+00  0.2646E+01
C--- 43  0.8501 -0.6327       -0.7476  0.0538 -0.2967E+01 -0.4132E-01  0.1817E+02
C--- 44  0.8726 -0.9983       -0.8139  0.0538 -0.2942E+01  0.1180E+01 -0.6774E+01
C--- 45  0.8874 -0.9082       -0.8574  0.0542 -0.2911E+01  0.8778E+00 -0.1364E+02
C--- 46  0.9139 -0.8930       -0.9340  0.0566 -0.2893E+01 -0.2044E+00 -0.8094E+01
C--- 47  0.9271 -1.0233       -0.9723  0.0593 -0.2903E+01 -0.5258E+00 -0.1498E+02
C--- 48  0.9473 -0.8839       -1.0313  0.0665 -0.2942E+01 -0.1433E+01  0.4945E+01
C--- 49  0.9652 -1.0172       -1.0843  0.0766 -0.2989E+01 -0.1168E+01  0.1401E+02
C--- 50  0.9930 -1.2715       -1.1679  0.0998
C---Documentation:
C---C   COMPUTER            - VAX/DOUBLE
C---C
C---C   AUTHOR              - M.F.HUTCHINSON
C---C                         CSIRO DIVISION OF MATHEMATICS AND STATISTICS
C---C                         P.O. BOX 1965
C---C                         CANBERRA, ACT 2601
C---C                         AUSTRALIA
C---C
C---C   LATEST REVISION     - 15 AUGUST 1985
C---C
C---C   PURPOSE             - CUBIC SPLINE DATA SMOOTHER
C---C
C---C   USAGE               - CALL CUBGCV (X,F,DF,N,Y,C,IC,VAR,JOB,SE,WK,IER)
C---C
C---C   ARGUMENTS    X      - VECTOR OF LENGTH N CONTAINING THE
C---C                           ABSCISSAE OF THE N DATA POINTS
C---C                           (X(I),F(I)) I=1..N. (INPUT) X
C---C                           MUST BE ORDERED SO THAT
C---C                           X(I) .LT. X(I+1).
C---C                F      - VECTOR OF LENGTH N CONTAINING THE
C---C                           ORDINATES (OR FUNCTION VALUES)
C---C                           OF THE N DATA POINTS (INPUT).
C---C                DF     - VECTOR OF LENGTH N. (INPUT/OUTPUT)
C---C                           DF(I) IS THE RELATIVE STANDARD DEVIATION
C---C                           OF THE ERROR ASSOCIATED WITH DATA POINT I.
C---C                           EACH DF(I) MUST BE POSITIVE.  THE VALUES IN
C---C                           DF ARE SCALED BY THE SUBROUTINE SO THAT
C---C                           THEIR MEAN SQUARE VALUE IS 1, AND UNSCALED
C---C                           AGAIN ON NORMAL EXIT.
C---C                           THE MEAN SQUARE VALUE OF THE DF(I) IS RETURNED
C---C                           IN WK(7) ON NORMAL EXIT.
C---C                           IF THE ABSOLUTE STANDARD DEVIATIONS ARE KNOWN,
C---C                           THESE SHOULD BE PROVIDED IN DF AND THE ERROR
C---C                           VARIANCE PARAMETER VAR (SEE BELOW) SHOULD THEN
C---C                           BE SET TO 1.
C---C                           IF THE RELATIVE STANDARD DEVIATIONS ARE UNKNOWN,
C---C                           SET EACH DF(I)=1.
C---C                N      - NUMBER OF DATA POINTS (INPUT).
C---C                           N MUST BE .GE. 3.
C---C                Y,C    - SPLINE COEFFICIENTS. (OUTPUT) Y
C---C                           IS A VECTOR OF LENGTH N. C IS
C---C                           AN N-1 BY 3 MATRIX. THE VALUE
C---C                           OF THE SPLINE APPROXIMATION AT T IS
C---C                           S(T)=((C(I,3)*D+C(I,2))*D+C(I,1))*D+Y(I)
C---C                           WHERE X(I).LE.T.LT.X(I+1) AND
C---C                           D = T-X(I).
C---C                IC     - ROW DIMENSION OF MATRIX C EXACTLY
C---C                           AS SPECIFIED IN THE DIMENSION
C---C                           STATEMENT IN THE CALLING PROGRAM. (INPUT)
C---C                VAR    - ERROR VARIANCE. (INPUT/OUTPUT)
C---C                           IF VAR IS NEGATIVE (I.E. UNKNOWN) THEN
C---C                           THE SMOOTHING PARAMETER IS DETERMINED
C---C                           BY MINIMIZING THE GENERALIZED CROSS VALIDATION
C---C                           AND AN ESTIMATE OF THE ERROR VARIANCE IS
C---C                           RETURNED IN VAR.
C---C                           IF VAR IS NON-NEGATIVE (I.E. KNOWN) THEN THE
C---C                           SMOOTHING PARAMETER IS DETERMINED TO MINIMIZE
C---C                           AN ESTIMATE, WHICH DEPENDS ON VAR, OF THE TRUE
C---C                           MEAN SQUARE ERROR, AND VAR IS UNCHANGED.
C---C                           IN PARTICULAR, IF VAR IS ZERO, THEN AN
C---C                           INTERPOLATING NATURAL CUBIC SPLINE IS CALCULATED.
C---C                           VAR SHOULD BE SET TO 1 IF ABSOLUTE STANDARD
C---C                           DEVIATIONS HAVE BEEN PROVIDED IN DF (SEE ABOVE).
C---C                JOB    - JOB SELECTION PARAMETER. (INPUT)
C---C                         JOB = 0 SHOULD BE SELECTED IF POINT STANDARD ERROR
C---C                           ESTIMATES ARE NOT REQUIRED IN SE.
C---C                         JOB = 1 SHOULD BE SELECTED IF POINT STANDARD ERROR
C---C                           ESTIMATES ARE REQUIRED IN SE.
C---C                SE     - VECTOR OF LENGTH N CONTAINING BAYESIAN STANDARD
C---C                           ERROR ESTIMATES OF THE FITTED SPLINE VALUES IN Y.
C---C                           SE IS NOT REFERENCED IF JOB=0. (OUTPUT)
C---C                WK     - WORK VECTOR OF LENGTH 7*(N + 2). ON NORMAL EXIT THE
C---C                           FIRST 7 VALUES OF WK ARE ASSIGNED AS FOLLOWS:-
C---C
C---C                           WK(1) = SMOOTHING PARAMETER (= RHO/(RHO + 1))
C---C                           WK(2) = ESTIMATE OF THE NUMBER OF DEGREES OF
C---C                                   FREEDOM OF THE RESIDUAL SUM OF SQUARES
C---C                           WK(3) = GENERALIZED CROSS VALIDATION
C---C                           WK(4) = MEAN SQUARE RESIDUAL
C---C                           WK(5) = ESTIMATE OF THE TRUE MEAN SQUARE ERROR
C---C                                   AT THE DATA POINTS
C---C                           WK(6) = ESTIMATE OF THE ERROR VARIANCE
C---C                           WK(7) = MEAN SQUARE VALUE OF THE DF(I)
C---C
C---C                           IF WK(1)=0 (RHO=0) AN INTERPOLATING NATURAL CUBIC
C---C                           SPLINE HAS BEEN CALCULATED.
C---C                           IF WK(1)=1 (RHO=INFINITE) A LEAST SQUARES
C---C                           REGRESSION LINE HAS BEEN CALCULATED.
C---C                           WK(2) IS AN ESTIMATE OF THE NUMBER OF DEGREES OF
C---C                           FREEDOM OF THE RESIDUAL WHICH REDUCES TO THE
C---C                           USUAL VALUE OF N-2 WHEN A LEAST SQUARES REGRESSION
C---C                           LINE IS CALCULATED.
C---C                           WK(3),WK(4),WK(5) ARE CALCULATED WITH THE DF(I)
C---C                           SCALED TO HAVE MEAN SQUARE VALUE 1.  THE
C---C                           UNSCALED VALUES OF WK(3),WK(4),WK(5) MAY BE
C---C                           CALCULATED BY DIVIDING BY WK(7).
C---C                           WK(6) COINCIDES WITH THE OUTPUT VALUE OF VAR IF
C---C                           VAR IS NEGATIVE ON INPUT.  IT IS CALCULATED WITH
C---C                           THE UNSCALED VALUES OF THE DF(I) TO FACILITATE
C---C                           COMPARISONS WITH A PRIORI VARIANCE ESTIMATES.
C---C
C---C                IER    - ERROR PARAMETER. (OUTPUT)
C---C                         TERMINAL ERROR
C---C                           IER = 129, IC IS LESS THAN N-1.
C---C                           IER = 130, N IS LESS THAN 3.
C---C                           IER = 131, INPUT ABSCISSAE ARE NOT
C---C                             ORDERED SO THAT X(I).LT.X(I+1).
C---C                           IER = 132, DF(I) IS NOT POSITIVE FOR SOME I.
C---C                           IER = 133, JOB IS NOT 0 OR 1.
C---C
C---C   PRECISION/HARDWARE  - DOUBLE
C---C
C---C   REQUIRED ROUTINES   - SPINT1,SPFIT1,SPCOF1,SPERR1
C---C
C---C   REMARKS      THE NUMBER OF ARITHMETIC OPERATIONS REQUIRED BY THE
C---C                SUBROUTINE IS PROPORTIONAL TO N.  THE SUBROUTINE
C---C                USES AN ALGORITHM DEVELOPED BY M.F. HUTCHINSON AND
C---C                F.R. DE HOOG, 'SMOOTHING NOISY DATA WITH SPLINE
C---C                FUNCTIONS', NUMER. MATH. (IN PRESS)
C---C
C---C-----------------------------------------------------------------------
C---C
C---
C---
C---ALGORITHM
C---
C---CUBGCV calculates a natural cubic spline curve which smoothes a given set
C---of data points, using statistical considerations to determine the amount
C---of smoothing required, as described in reference 2.  If the error variance
C---is known, it should be supplied to the routine in VAR.  The degree of
C---smoothing is then determined by minimizing an unbiased estimate of the true
C---mean square error.  On the other hand, if the error variance is not known, VAR
C---should be set to -1.0.  The routine then determines the degree of smoothing
C---by minimizing the generalized cross validation.  This is asymptotically the
C---same as minimizing the true mean square error (see reference 1).  In this
C---case, an estimate of the error variance is returned in VAR which may be
C---compared with any a priori approximate estimates.  In either case, an
C---estimate of the true mean square error is returned in WK(5).  This estimate,
C---however, depends on the error variance estimate, and should only be accepted
C---if the error variance estimate is reckoned to be correct.
C---
C---If job is set to 1, bayesian estimates of the standard error of each
C---smoothed data value are returned in the array SE.  These also depend on
C---the error variance estimate and should only be accepted if the error
C---variance estimate is reckoned to be correct.  See reference 4.
C---
C---The number of arithmetic operations and the amount of storage required by
C---the routine are both proportional to N, so that very large data sets may be
C---analysed.  The data points do not have to be equally spaced or uniformly
C---weighted.  The residual and the spline coefficients are calculated in the
C---manner described in reference 3, while the trace and various statistics,
C---including the generalized cross validation, are calculated in the manner
C---described in reference 2.
C---
C---When VAR is known, any value of N greater than 2 is acceptable.  It is
C---advisable, however, for N to be greater than about 20 when VAR is unknown.
C---If the degree of smoothing done by CUBGCV when VAR is unknown is not
C---satisfactory, the user should try specifying the degree of smoothing by
C---setting VAR to a reasonable value.
C---
C---References:
C---
C---1.  Craven, Peter and Wahba, Grace, "Smoothing noisy data with spline
C---    functions", Numer. Math. 31, 377-403 (1979).
C---2.  Hutchinson, M.F. and de Hoog, F.R., "Smoothing noisy data with spline
C---    functions", Numer. Math. (in press).
C---3.  Reinsch, C.H., "Smoothing by spline functions", Numer. Math. 10,
C---    177-183 (1967).
C---4.  Wahba, Grace, "Bayesian 'confidence intervals' for the cross-validated
C---    smoothing spline", J.R.Statist. Soc. B 45, 133-150 (1983).
C---
C---
C---Example
C---
C---A sequence of 50 data points are generated by adding a random variable with
C---uniform density in the interval [-0.3,0.3] to the curve y=sin(3*pi*x/2).
C---The abscissae are unequally spaced in [0,1].  Point standard error estimates
C---are returned in SE by setting JOB to 1.  The error variance estimate is
C---returned in VAR.  It compares favourably with the true value of 0.03.
C---The IMSL function GGUBFS is used to generate sample values of a uniform
C---variable on [0,1].
C---
C---
C---INPUT:
C---
C---      INTEGER          N,IC,JOB,IER
C---      DOUBLE PRECISION X(50),F(50),Y(50),DF(50),C(49,3),WK(400),
C---     *                 VAR,SE(50)
C---      DOUBLE PRECISION GGUBFS,DSEED,DN
C---      DATA DSEED/1.2345D4/
C---C
C---      N=50
C---      IC=49
C---      JOB=1
C---      VAR=-1.0D0
C---      DN=N
C---C
C---      DO 10 I=1,N
C---      X(I)=(I - 0.5)/DN + (2.0*GGUBFS(DSEED) - 1.0)/(3.0*DN)
C---      F(I)=DSIN(4.71238*X(I)) + (2.0*GGUBFS(DSEED) - 1.0)*0.3
C---      DF(I)=1.0D0
C---  10  CONTINUE
C---      CALL CUBGCV(X,F,DF,N,Y,C,IC,VAR,JOB,SE,WK,IER)
C---       .
C---       .
C---       .
C---      END
C---
C---OUTPUT:
C---
C---IER = 0
C---VAR = 0.0279
C---GENERALIZED CROSS VALIDATION = 0.0318
C---MEAN SQUARE RESIDUAL = 0.0246
C---RESIDUAL DEGREES OF FREEDOM = 43.97
C---FOR CHECKING PURPOSES THE FOLLOWING OUTPUT IS GIVEN:
C---
C---X(1)  = 0.0046    F(1)  =  0.2222     Y(1)  =  0.0342     SE(1)  = 0.1004
C---X(21) = 0.4126    F(21) =  0.9810     Y(21) =  0.9069     SE(21) = 0.0525
C---X(41) = 0.8087    F(41) = -0.5392     Y(41) = -0.6242     SE(41) = 0.0541
C---

       PROGRAM AVGFLX
C t3e:  f90 -o xavg avgflux.f
c arc:  lf95 -o xavg avgflux.f
       PARAMETER (MXT=10000,MXFIN=8,MXPIN=8)
       REAL TFLX(MXT),TPFLX(MXT),QI(MXT),QE(MXT),SUMI(MXT),SUME(MXT),
     1 DELT(MXT),PI(MXT),PE(MXT),ZTOT(MXT),SUMZT(MXT),ZTOTF(MXT),
     2 SUMZTF(MXT)
       integer ncarg,IUNIN/17/
       INTEGER ITTYWR/6/
C       INTEGER ITTYRD/5/,ITTYWR/6/
       CHARACTER*20 INFLX,run_name,time1,time2,short
       LOGICAL LOUT_SHORT/.FALSE./
C
       ncarg=iargc()
C LIB$SYS_TRNLOG is an obsolete System Service routine. Invoking it as
C a function returns a success code of 1 (1577: nonexistent logical name)
ccc       NTEST=LIB$SYS_TRNLOG('AVG_RUN',NTEST2,RUN_NAME)
cdbg       WRITE(6,'('' NTEST='',I5)') NTEST
c	NCARG=0
c       IF(NTEST.EQ.1) NCARG=1
       IF(NCARG.LT.1) THEN
         WRITE(ITTYWR,
     x   '('' Need command line arguments:'',/,
     2 ''%xavg run_name [time1] [time2] [oneline]'',/,
     3 ''  run_name (e.g, g14) is required'',/,
     4 ''  time1 is the averaging start time [optional]'',/,
     5 ''  time2 is the averaging end time [optional]'',/,
     6 ''  "oneline" reduces output for bulk processing'',/,
     7 '' Output:'',/,''  <tot> averages over all times'',/,
     8 ''  tdef is an automatically determined start time'',/,
     9 ''  <def> is average for tdef < t < end of file'',/,
     x ''  t1 & t2 are user input time1 & time2'',/,
     1 ''  <1-2> is average for t1 < t < t2'',/,
     2 ''  ts & te are start and end times for final averages'')')
         stop
       ENDIF
c t3e    CALL pxfgetarg(1,run_name,ilen,ierr)
       CALL getarg(1,run_name)
       IBLNK=INDEX(run_name,' ')
       INFLX=run_name(1:IBLNK-1)//'.flux'
       OPEN(UNIT=IUNIN,file=INFLX,status='old')
ccc       NTEST=LIB$SYS_TRNLOG('AVG_T1',NTEST2,time1)
cdbg       WRITE(6,'('' NTEST='',I5)') NTEST
ccc       IF(NTEST.EQ.1) NCARG=2
       IF(NCARG.GE.2) then
c t3e         CALL pxfgetarg(2,time1,ilen,ierr)
         CALL getarg(2,time1)
         IF(TIME1.EQ.'oneline') THEN
           LOUT_SHORT=.TRUE.
           T1=0.
         ELSE
           READ(time1,*) T1
         ENDIF
       ELSE
         T1=0.
       ENDIF
ccc       NTEST=LIB$SYS_TRNLOG('AVG_T2',NTEST2,time2)
cdbg       WRITE(6,'('' NTEST='',I5)') NTEST
ccc       IF(NTEST.EQ.1) NCARG=3
       IF(NCARG.GE.3) then
c t3e         CALL pxfgetarg(3,time2,ilen,ierr)
         CALL getarg(3,time2)
         IF(TIME2.EQ.'oneline') THEN
           LOUT_SHORT=.TRUE.
           T2=0.
         ELSE
           READ(time2,*) T2
         ENDIF
       ELSE
         T2=0.
       ENDIF
ccc       NTEST=LIB$SYS_TRNLOG('AVG_ONE',NTEST2,SHORT)
cdbg       WRITE(6,'('' NTEST='',I5)') NTEST
ccc       IF(NTEST.EQ.1) NCARG=4
       IF(NCARG.GE.4) then
c t3e      CALL pxfgetarg(4,short,ilen,ierr)
         CALL getarg(4,short)
         IF(short.EQ.'oneline') LOUT_SHORT=.TRUE.
         IF(short.EQ.'ONELINE') LOUT_SHORT=.TRUE.
       ENDIF
C
       ICOUNT=0
 100   CONTINUE
C12345678(1)2345678(2)2345678(3)2345678(4)2345678(5)2345678(6)2345678(7)
       READ(IUNIN,'(T3,E17.10,T55,E11.4,T66,E11.4)',
     X ERR=900,END=200) ZTIME,ZQI,ZQE
       ICOUNT=ICOUNT+1
       IF(ICOUNT.GT.MXT) THEN
         WRITE(ITTYWR,'('' Arrays too short; need larger MXT'')')
         stop
       ENDIF
       TFLX(ICOUNT)=ZTIME
       QI(ICOUNT)=ZQI
       QE(ICOUNT)=ZQE
       GO TO 100
 200   CONTINUE
       IF(T1.GT.TFLX(1)) THEN
         IT1=0
       ELSE
         IT1=1
       ENDIF
       IF(T2.LT.TFLX(ICOUNT)) THEN
         IT2=0
         IF(T2.LT.TFLX(IT1)) IT2=ICOUNT
         IF(T2.LT.0.) IT2=ICOUNT
       ELSE
         IT2=ICOUNT
       ENDIF
       SUMI(1)=0.
       SUME(1)=0.
       DO J=1,ICOUNT-1
         DELT(J)=TFLX(J+1)-TFLX(J)
         SUMI(J+1)=SUMI(J)+QI(J)*DELT(J)
         SUME(J+1)=SUME(J)+QE(J)*DELT(J)
         IF(IT1.EQ.0 .AND. TFLX(J).LE.T1 .AND. TFLX(J+1).GT.T1) IT1=J
         IF(IT2.EQ.0 .AND. TFLX(J).LE.T2 .AND. TFLX(J+1).GT.T2) IT2=J
       ENDDO
 210   CONTINUE
       AVGTOTI=SUMI(ICOUNT)/(TFLX(ICOUNT)-TFLX(1))
       AVGTOTE=SUME(ICOUNT)/(TFLX(ICOUNT)-TFLX(1))
       IF(IT2.EQ.0) IT2=ICOUNT
       AVGUI=(SUMI(IT2)-SUMI(IT1))/(TFLX(IT2)-TFLX(IT1))
       AVGUE=(SUME(IT2)-SUME(IT1))/(TFLX(IT2)-TFLX(IT1))
C
       QIPEAK=0.
       ITIPEAK=1
       ITIDEF=1
       DO J=1,ICOUNT-1
         AVGDEFI=(SUMI(ICOUNT)-SUMI(J))/(TFLX(ICOUNT)-TFLX(J))
C  Don't look for default time until the initial peak has been passed
C  if this is a run from t=0:
         IF(TFLX(1).LT.10.) THEN
           IF(QIPEAK.LE.0. .AND. QI(J).GT.10. 
     1       .AND. QI(J+1).LT.QI(J)) THEN
             QIPEAK=QI(J)
             ITIPEAK=J
           ENDIF
         ELSE
           QIPEAK=1.
         ENDIF
         IF(QIPEAK.GT.0. .AND. QI(J).LT.AVGDEFI) THEN
           TDEFI=TFLX(J)
           ITIDEF=J
           GO TO 220
         ENDIF
       ENDDO
 220   CONTINUE
       QEPEAK=0.
       ITEDEF=1
       DO J=1,ICOUNT-1
         AVGDEFE=(SUME(ICOUNT)-SUME(J))/(TFLX(ICOUNT)-TFLX(J))
C  Don't look for default time until the initial peak has been passed
C  if this is a run from t=0:
         IF(TFLX(1).LT.10.) THEN
           IF(QEPEAK.LE.0. .AND. QE(J).GT.10. .AND.
     1        QE(J+1).LT.QE(J)) QEPEAK=QE(J)
         ELSE
           QEPEAK=1.
         ENDIF
         IF(QEPEAK.GT.0. .AND. QE(J).LT.AVGDEFE) THEN
           TDEFE=TFLX(J)
           ITEDEF=J
           GO TO 230
         ENDIF
       ENDDO
 230   CONTINUE
       IF(LOUT_SHORT) GO TO 299
       WRITE(ITTYWR,'(5X,''<tot>  tdef  <def>    t1   t2  <1-2>'')')
       WRITE(ITTYWR,'('' Qi'',0PF7.2,2X,I4,0PF7.2,2X,I4,X,I4,0PF7.2)')
     1 AVGTOTI,IFIX(TDEFI+0.5),AVGDEFI,IFIX(TFLX(IT1)+0.5),
     2 IFIX(TFLX(IT2)+0.5),AVGUI
       WRITE(ITTYWR,'('' Qe'',0PF7.2,2X,I4,0PF7.2,2X,I4,X,I4,0PF7.2)')
     1 AVGTOTE,IFIX(TDEFE+0.5),AVGDEFE,IFIX(TFLX(IT1)+0.5),
     2 IFIX(TFLX(IT2)+0.5),AVGUE
C
 299   CONTINUE
       CLOSE(IUNIN)
       PI=0.
       PE=0.
       INFLX=run_name(1:IBLNK-1)//'.plux'
       OPEN(UNIT=IUNIN,file=INFLX,status='old',err=400)
       ICNTP=0
 300   CONTINUE
C12345678(1)2345678(2)2345678(3)2345678(4)2345678(5)2345678(6)2345678(7)
       READ(IUNIN,'(T3,E17.10,T55,E11.4,T66,E11.4)',
     X ERR=900,END=400) ZTIME,ZPI,ZPE
       ICNTP=ICNTP+1
       IF(ICNTP.GT.MXT) THEN
         WRITE(ITTYWR,'('' Arrays too short; need larger MXT'')')
         stop
       ENDIF
       TPFLX(ICNTP)=ZTIME
       PI(ICNTP)=ZPI
       PE(ICNTP)=ZPE
       ZTOT(ICNTP)=ABS(ZPI-ZPE)
       ZTOTF(ICNTP)=ABS(ZPI-ZPE)/AMAX1(ABS(ZPE),ABS(ZPI))
       GO TO 300
 400   CONTINUE
       SUMI(1)=0.
       SUME(1)=0.
       SUMZT(1)=0.
       SUMZTF(1)=0.
       DO J=1,ICNTP-1
         DELT(J)=TPFLX(J+1)-TPFLX(J)
         SUMI(J+1)=SUMI(J)+PI(J)*DELT(J)
         SUME(J+1)=SUME(J)+PE(J)*DELT(J)
         SUMZT(J+1)=SUMZT(J)+ZTOT(J)*DELT(J)
         SUMZTF(J+1)=SUMZTF(J)+ZTOTF(J)*DELT(J)
       ENDDO
       ITP1=IT1
       ITP2=MIN0(ICNTP,IT2)
       AVGTPI=SUMI(ICNTP)/(TPFLX(ICNTP)-TPFLX(1))
       AVGTPE=SUME(ICNTP)/(TPFLX(ICNTP)-TPFLX(1))
       AVGUPI=(SUMI(ITP2)-SUMI(ITP1))/(TPFLX(ITP2)-TPFLX(ITP1))
       AVGUPE=(SUME(ITP2)-SUME(ITP1))/(TPFLX(ITP2)-TPFLX(ITP1))
       AVGTZT=SUMZT(ICNTP)/(TPFLX(ICNTP)-TPFLX(1))
       AVGUZT=(SUMZT(ITP2)-SUMZT(ITP1))/(TPFLX(ITP2)-TPFLX(ITP1))
       AVGTZTF=SUMZTF(ICNTP)/(TPFLX(ICNTP)-TPFLX(1))
       AVGUZTF=(SUMZTF(ITP2)-SUMZTF(ITP1))/(TPFLX(ITP2)-TPFLX(ITP1))
C
       DO J=ITIDEF,ICNTP-1
         AVGDPI=(SUMI(ICNTP)-SUMI(J))/(TPFLX(ICNTP)-TPFLX(J))
         IF(PI(J).LT.AVGDPI) THEN
           TDEFPI=TPFLX(J)
           AVGDZT=(SUMZT(ICNTP)-SUMZT(J))/(TPFLX(ICNTP)-TPFLX(J))
           AVGDZTF=(SUMZTF(ICNTP)-SUMZTF(J))/(TPFLX(ICNTP)-TPFLX(J))
           GO TO 420
         ENDIF
       ENDDO
 420   CONTINUE
       DO J=1,ICNTP-1
         AVGDPE=(SUME(ICNTP)-SUME(J))/(TPFLX(ICNTP)-TPFLX(J))
         IF(PE(J).LT.AVGDPE) THEN
           TDEFPE=TPFLX(J)
           GO TO 430
         ENDIF
       ENDDO
 430   CONTINUE
       IF(LOUT_SHORT) GO TO 499
       WRITE(ITTYWR,'('' Pi'',0PF7.2,2X,I4,0PF7.2,2X,I4,X,I4,0PF7.2)')
     1 AVGTPI,IFIX(TDEFPI+0.5),AVGDPI,IFIX(TPFLX(IT1)+0.5),
     2 IFIX(TPFLX(IT2)+0.5),AVGUPI
       WRITE(ITTYWR,'('' Zt'',0PF7.2,2X,I4,0PF7.2,2X,I4,X,I4,0PF7.2)')
     1 AVGTZT,IFIX(TDEFPI+0.5),AVGDZT,IFIX(TPFLX(IT1)+0.5),
     2 IFIX(TPFLX(IT2)+0.5),AVGUZT
       WRITE(ITTYWR,'('' Zf'',0PF7.2,2X,I4,0PF7.2,2X,I4,X,I4,0PF7.2)')
     1 AVGTZTF,IFIX(TDEFPI+0.5),AVGDZTF,IFIX(TPFLX(IT1)+0.5),
     2 IFIX(TPFLX(IT2)+0.5),AVGUZTF
C
       WRITE(ITTYWR,
     1 '(/,3X,'' ts   te    Qi     Qe     Pi'',27X,''Zt/Pi   Ztf'')')
C
 499   CONTINUE
       IF(T1.LE.0) THEN
         ST1=TDEFI
         ST2=TFLX(ICOUNT)
         SAVGI=AVGDEFI
         SAVGE=AVGDEFE
         SAVGPI=AVGDPI
         SAVGZT=AVGDZT
         SAVGZTF=AVGDZTF
         CALL GETRMS(TFLX,QI,ITIDEF,ICOUNT,TAVGI,TRMSI)
         CALL GETRMS(TFLX,QE,ITEDEF,ICOUNT,TAVGE,TRMSE)
         CALL GETRMS(TPFLX,PI,ITIDEF,ICNTP,TAVGP,TRMSP)
         CALL GETRMS(TPFLX,ZTOT,ITIDEF,ICNTP,TAVGZT,TRMSZT)
         CALL GETRMS(TPFLX,ZTOTF,ITIDEF,ICNTP,TAVGZTF,TRMSZTF)
       ELSE
         ST1=TFLX(IT1)
         ST2=TFLX(IT2)
         SAVGI=AVGUI
         SAVGE=AVGUE
         SAVGPI=AVGUPI
         SAVGZT=AVGUZT
         SAVGZTF=AVGUZTF
         CALL GETRMS(TFLX,QI,IT1,IT2,TAVGI,TRMSI)
         CALL GETRMS(TFLX,QE,IT1,IT2,TAVGE,TRMSE)
         CALL GETRMS(TPFLX,PI,IT1,IT2,TAVGP,TRMSP)
         CALL GETRMS(TPFLX,ZTOT,IT1,IT2,TAVGZT,TRMSZT)
         CALL GETRMS(TPFLX,ZTOTF,IT1,IT2,TAVGZTF,TRMSZTF)
       ENDIF
c       WRITE(ITTYWR,'(2X,I4,X,I4,X,0P3F7.2,2X,A)')
c     1 IFIX(ST1+0.5),IFIX(ST2+0.5),SAVGI,SAVGE,SAVGPI,RUN_NAME
       WRITE(ITTYWR,'(2X,I4,X,I4,X,0P3F7.2,2X,A,X,0P2F7.3)')
     1 IFIX(ST1+0.5),IFIX(ST2+0.5),TAVGI,TAVGE,TAVGP,RUN_NAME,
     2 TAVGZT/TAVGP,TAVGZTF
       IF(.NOT.LOUT_SHORT)WRITE(ITTYWR,
     1 '('' rms var:   '',0P3F7.2,23X,0P2F7.3)')
     2  TRMSI,TRMSE,TRMSP,TRMSZT/TAVGP,TRMSZTF
C
       stop
 900   CONTINUE
       WRITE(ITTYWR,'('' ERROR reading file '',A20)') INFLX
       stop
       END
C/ MODULE GETRMS
      SUBROUTINE GETRMS(PRTIM,PRFUN,KSTART,KEND,PAVG,PRMS)
C
      INTEGER I,KSTART,KEND
      REAL PRTIM(KEND),PRFUN(KEND),PRMS
C
      PAVG=0.
      DO I=KSTART,KEND-1
        PAVG=PAVG+0.5*(PRFUN(I)+PRFUN(I+1))*(PRTIM(I+1)-PRTIM(I))/
     1  (PRTIM(KEND)-PRTIM(KSTART))
      ENDDO
C
      PRMS=0.
      DO I=KSTART,KEND-1
        PRMS=PRMS+(0.5*(PRFUN(I)+PRFUN(I+1))-PAVG)**2*
     1  (PRTIM(I+1)-PRTIM(I))/(PRTIM(KEND)-PRTIM(KSTART))
      ENDDO
      PRMS=SQRT(PRMS) ! root mean square deviation.
      RETURN
      END


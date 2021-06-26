CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCC wrappers implemented by ANGEL GUTARRA-LEON CCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



      SUBROUTINE TRACE_08_W (x,y,z, datetime, PARMOD,
     *  EXNAME, INNAME,
     *  DIR, DSMAX, ERR, RLIM, R0, IOPT,
     *  XF, YF, ZF, XX, YY, ZZ, L)

!f2py intent(in) :: x,y,z,PARMOD, datetime
!f2py real optional, intent(in) :: DIR=1. ,DSMAX=1.0 ,ERR=0.0001 ,RLIM=60.
!f2py real optional, intent(in) :: R0=1.,IOPT=0
!f2py character(30) optional, intent(in) :: EXNAME="t96_01", INNAME="dip_08"
!f2py intent(out) :: XX, YY, ZZ
!f2py intent(hide) :: L,  XF, YF, ZF
      CHARACTER(30) :: EXNAME, INNAME
      REAL*8 :: X, Y, Z, DIR, DSMAX, ERR, RLIM, R0,XF,YF,ZF,
     * XX(500), YY(500), ZZ(500), PARMOD(10)
      INTEGER*4 :: datetime(5), IOPT
      CALL RECALC_08_W (datetime)
C
      CALL TRACE_08 (x,y,z,DIR,DSMAX,ERR,RLIM,R0,IOPT,
     * PARMOD,EXNAME,INNAME,XF,YF,ZF,XX,YY,ZZ,L,500)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE CONS_MAG(PARMOD,bx, by, bz)

!f2py intent(in) :: PARMOD
!f2py intent(out) :: bx, by, bz
        REAL*8 :: PARMOD(10), bx, by, bz

        bx = PARMOD(1)
        by = PARMOD(2)
        bz = PARMOD(3)
        RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE RECALC_08_W (datetime)
        COMMON /GEOPACK1/ AA(10),SPS,CPS,BB(3),PSI,CC(18)
        INTEGER*4 IYEAR,IDAY,IHOUR,MIN,ISEC, datetime(5)
        REAL*8 :: vx, vy,vz
!f2py intent(in) :: datetime
        vx=-400
        vy=0
        vz=0
        IYEAR = datetime(1)
        IDAY = datetime(2)
        IHOUR = datetime(3)
        min = datetime(4)
        ISEC = datetime(5)

C  RECALC_08 prepares elements of rot matrix and puts in common block
        CALL RECALC_08(IYEAR,IDAY,IHOUR,MIN,ISEC,vx,vy,vz)

        RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE SPHtoCAR_08_V (R,THETA,PHI,X,Y,Z,N)

        REAL*8 R(N), THETA(N), PHI(N), X(N), Y(N), Z(N)
        INTEGER N
!f2py intent(in) :: R, THETA, PHI
!f2py intent(out) :: X, Y, Z
!f2py intent(hide) :: N

        DO i=1, N
          CALL SPHCAR_08 (R(i),THETA(i),PHI(i),X(i),Y(i),Z(i),1)
        END DO
        RETURN
        END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        SUBROUTINE CARtoSPH_08_V (R,THETA,RHO,X,Y,Z,N)

          REAL*8 R(N), THETA(N), RHO(N), X(N), Y(N), Z(N)
          INTEGER N
!f2py intent(in) ::  X, Y, Z
!f2py intent(out) :: R, THETA, RHO
!f2py intent(hide) :: N
          DO  i=1, N
            CALL SPHCAR_08 (R(i),THETA(i),RHO(i),X(i),Y(i),Z(i),-1)
          END DO
          RETURN
          END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE transform (Xin,Yin,Zin, Xout,Yout,Zout,
     1                       trans, datetime,n,mod_psi,cpsi )
C
C----------------------------------------------------------------------
C     Uses: transforms vector from one coordinate system to another one
C           using geopack. Geopack was created by Tsyganenko.
C
C     List of available coordinate systems:
C        GEO - Geocentric
C        MAG - geomagnetic
C        SM  - solar magnetic
C        GSM - Geocentric Solar Magnetic
C        GSE - geocentric solar ecliptic
C        GEI - geocentric intertial
C
C
C-----INPUT PARAMETERS:
C     Xin, Yin, Zin : 1xN arrays of floats in the old coordinate system
C     datetime : 1x5 array of floats
C       datetime(0)   -  YEAR NUMBER (FOUR DIGITS)
C       datetme(1)  -  DAY OF YEAR (DAY 1 = JAN 1)
C       datetime(2) -  HOUR OF DAY (00 TO 23)
C       datetime(3)   -  MINUTE OF HOUR (00 TO 59)
C       datetime(4)  -  SECONDS OF MINUTE (00 TO 59)
C     trans : string that specifies the coordinate transformation.
C             Follows the pattern "GEOtoMAG".
C
C-----OUTPUT PARAMETERS:
C     Xout, Yout, Zout : 1xN arrays of floats in the new coordinate system
C
C
c Tsyganenko uses J to determine the direction of the transformation.
c For examples take SMGSW_08, if J is positive it converts from SM to GSW,
c but if J is negative it goes from GSW to SM. positive J always takes the
c first coordinate to be the input and the second to be the output. The roles
c reverse for regular J. The user of this script does not need to specify J
c because a series of "if-else" statements will determine the proper choice of
c J based on the "trans" keyword i.e. "SM2GSW"

c Tsyganenko uses Geocentric Solar Wind GSW coordinate system instead of
c Geocentric Solar Magnetic GSM. GSW becomes identical to GSM when the solar
c wind velocity's only component is in the -x-direction of the GSE coordinate system.
c This subroutine will always set the solar wind velocity so that the GSW
c becomes identical to the GSM coordinate system.

      COMMON /GEOPACK1/ AA(10),SPS,CPS,BB(3),PSI,CC(18)
      INTEGER*4 :: N, ind,datetime(5), mod_psi
      REAL*8 :: SPS, CPS, PSI, cpsi,
     *  Xin(N), Yin(N), Zin(N), xtmp(N), ytmp(N), ztmp(N),
     *  xtmp2(N), ytmp2(N), ztmp2(N),
     *  Xout(N), Yout(N), Zout(N)
      CHARACTER(10) :: trans
!f2py intent(in) ::  Xin, Yin, Zin, trans, datetime
!f2py optional, intent(in):: mod_psi = 1, cpsi=0
!f2py intent(hide) :: N
!f2py intent(out) ::  Xout, Yout, Zout

      N = SIZE(Xin)
C  RECALC_08 prepares elements of rot matrix and puts in common block
      CALL RECALC_08_W (datetime)
      IF (mod_psi > 0) THEN
        PSI = cpsi * 3.1415/180.
        SPS = SIN(PSI)
        CPS = COS(PSI)
      END IF

      DO ind=1,N

c these transformations only rely on one call.
      IF      (trans=='GEItoGEO') THEN
        CALL GEIGEO_08 (Xin(ind),Yin(ind),Zin(ind),
     1               Xout(ind),Yout(ind),Zout(ind), 1)
      ELSE IF (trans=='GEOtoGEI') THEN
        CALL GEIGEO_08 (Xout(ind),Yout(ind),Zout(ind),
     1               Xin(ind),Yin(ind),Zin(ind), -1)

      ELSE IF (trans=='GEOtoGSM') THEN
        CALL GEOGSW_08 (Xin(ind),Yin(ind),Zin(ind),
     1               Xout(ind),Yout(ind),Zout(ind), 1)
      ELSE IF (trans=='GSMtoGEO') THEN
         CALL GEOGSW_08 (Xout(ind),Yout(ind),Zout(ind),
     1               Xin(ind),Yin(ind),Zin(ind), -1)

      ELSE IF (trans=='GEOtoMAG') THEN
        CALL GEOMAG_08 (Xin(ind),Yin(ind),Zin(ind),
     1               Xout(ind),Yout(ind),Zout(ind), 1)
      ELSE IF (trans=='MAGtoGEO') THEN
        CALL GEOMAG_08 (Xout(ind),Yout(ind),Zout(ind),
     1               Xin(ind),Yin(ind),Zin(ind), -1)

      ELSE IF (trans=='GSMtoGSE') THEN
        CALL GSWGSE_08 (Xin(ind),Yin(ind),Zin(ind),
     1               Xout(ind),Yout(ind),Zout(ind), 1)
      ELSE IF (trans=='GSEtoGSM') THEN
        CALL GSWGSE_08 (Xout(ind),Yout(ind),Zout(ind),
     1               Xin(ind),Yin(ind),Zin(ind), -1)

      ELSE  IF (trans=='MAGtoSM') THEN
        CALL MAGSM_08 (Xin(ind),Yin(ind),Zin(ind),
     1               Xout(ind),Yout(ind),Zout(ind), 1)
      ELSE IF (trans=='SMtoMAG') THEN
        CALL MAGSM_08 (Xout(ind),Yout(ind),Zout(ind),
     1               Xin(ind),Yin(ind),Zin(ind), -1)

      ELSE IF (trans=='SMtoGSM') THEN
        CALL SMGSW_08 (Xin(ind),Yin(ind),Zin(ind),
     1                 Xout(ind),Yout(ind),Zout(ind), 1)
      ELSE IF (trans=='GSMtoSM') THEN
        CALL SMGSW_08 (Xout(ind),Yout(ind),Zout(ind),
     1               Xin(ind),Yin(ind),Zin(ind), -1)

c these transformations rely on multiple calls.
      ELSE IF (trans=='MAGtoGSM') THEN
        CALL MAGSM_08 (Xin(ind),Yin(ind),Zin(ind),
     1               xtmp(ind),ytmp(ind),ztmp(ind), 1)
        CALL SMGSW_08 (xtmp(ind),ytmp(ind),ztmp(ind),
     1               Xout(ind),Yout(ind),Zout(ind), 1)
      ELSE IF (trans=='GSMtoMAG') THEN
        CALL SMGSW_08 (xtmp(ind),ytmp(ind),ztmp(ind),
     1               Xin(ind),Yin(ind),Zin(ind), -1)
        CALL MAGSM_08 (Xout(ind),Yout(ind),Zout(ind),
     1               xtmp(ind),ytmp(ind),ztmp(ind), -1)


      ELSE IF (trans=='GEItoMAG') THEN
        CALL GEIGEO_08 (Xin(ind),Yin(ind),Zin(ind),
     1               xtmp(ind),ytmp(ind),ztmp(ind), 1)
        CALL GEOMAG_08 (xtmp(ind),ytmp(ind),ztmp(ind),
     1               Xout(ind),Yout(ind),Zout(ind), 1)
      ELSE IF (trans=='MAGtoGEI') THEN
        CALL GEOMAG_08 (xtmp(ind),ytmp(ind),ztmp(ind),
     1               Xin(ind),Yin(ind),Zin(ind), -1)
        CALL GEIGEO_08 (Xout(ind),Yout(ind),Zout(ind),
     1               xtmp(ind),ytmp(ind),ztmp(ind), -1)


      ELSE IF (trans=='GEOtoSM') THEN
        CALL GEOMAG_08 (Xin(ind),Yin(ind),Zin(ind),
     1               xtmp(ind),ytmp(ind),ztmp(ind), 1)
        CALL MAGSM_08 (xtmp(ind),ytmp(ind),ztmp(ind),
     1               Xout(ind),Yout(ind),Zout(ind), 1)
      ELSE IF (trans=='SMtoGEO') THEN
        CALL MAGSM_08 (xtmp(ind),ytmp(ind),ztmp(ind),
     1               Xin(ind),Yin(ind),Zin(ind), -1)
        CALL GEOMAG_08 (Xout(ind),Yout(ind),Zout(ind),
     1                xtmp(ind),ytmp(ind),ztmp(ind),-1)


      ELSE IF (trans=='SMtoGSE') THEN
        CALL SMGSW_08 (Xin(ind),Yin(ind),Zin(ind),
     1               xtmp(ind),ytmp(ind),ztmp(ind), 1)
        CALL GSWGSE_08 (xtmp(ind),ytmp(ind),ztmp(ind),
     1               Xout(ind),Yout(ind),Zout(ind), 1)
      ELSE IF (trans=='GSEtoSM') THEN
        CALL GSWGSE_08 (xtmp(ind),ytmp(ind),ztmp(ind),
     1               Xin(ind),Yin(ind),Zin(ind), -1)
        CALL SMGSW_08 (Xout(ind),Yout(ind),Zout(ind),
     1               xtmp(ind),ytmp(ind),ztmp(ind), -1)


      ELSE IF (trans=='GEOtoGSE') THEN
        CALL GEOGSW_08 (Xin(ind),Yin(ind),Zin(ind),
     1               xtmp(ind),ytmp(ind),ztmp(ind), 1)
        CALL GSWGSE_08 (xtmp(ind),ytmp(ind),ztmp(ind),
     1               Xout(ind),Yout(ind),Zout(ind), 1)
      ELSE IF (trans=='GSEtoGEO') THEN
        CALL GSWGSE_08 (xtmp(ind),ytmp(ind),ztmp(ind),
     1               Xin(ind),Yin(ind),Zin(ind), -1)
        CALL GEOGSW_08 (Xout(ind),Yout(ind),Zout(ind),
     1               xtmp(ind),ytmp(ind),ztmp(ind), -1)


      ELSE IF (trans=='GEItoGSM') THEN
        CALL GEIGEO_08 (Xin(ind),Yin(ind),Zin(ind),
     1               xtmp(ind),ytmp(ind),ztmp(ind), 1)
        CALL GEOGSW_08 (xtmp(ind),ytmp(ind),ztmp(ind),
     1               Xout(ind),Yout(ind),Zout(ind), 1)
      ELSE IF (trans=='GSMtoGEI') THEN
        CALL GEOGSW_08 (xtmp(ind),ytmp(ind),ztmp(ind),
     1               Xin(ind),Yin(ind),Zin(ind), -1)
        CALL GEIGEO_08 (Xout(ind),Yout(ind),Zout(ind),
     1               xtmp(ind),ytmp(ind),ztmp(ind), -1)


      ELSE IF (trans=='GEItoSM') THEN
        CALL GEIGEO_08 (Xin(ind),Yin(ind),Zin(ind),
     1               xtmp(ind),ytmp(ind),ztmp(ind), 1)
        CALL GEOMAG_08 (xtmp(ind),ytmp(ind),ztmp(ind),
     1               xtmp2(ind),ytmp2(ind),ztmp2(ind), 1)
        CALL MAGSM_08 (xtmp2(ind),ytmp2(ind),ztmp2(ind),
     1               Xout(ind),Yout(ind),Zout(ind), 1)
      ELSE IF (trans=='SMtoGEI') THEN
        CALL MAGSM_08 (xtmp(ind),ytmp(ind),ztmp(ind),
     1               Xin(ind),Yin(ind),Zin(ind), -1)
        CALL GEOMAG_08 (xtmp2(ind),ytmp2(ind),ztmp2(ind),
     1               xtmp(ind),ytmp(ind),ztmp(ind), -1)
        CALL GEIGEO_08 (Xout(ind),Yout(ind),Zout(ind),
     1               xtmp2(ind),ytmp2(ind),ztmp2(ind), -1)


      ELSE IF (trans=='GEItoGSE') THEN
        CALL GEIGEO_08 (Xin(ind),Yin(ind),Zin(ind),
     1              xtmp(ind),ytmp(ind),ztmp(ind), 1)
        CALL GEOGSW_08 (xtmp(ind),ytmp(ind),ztmp(ind),
     1               xtmp2(ind),ytmp2(ind),ztmp2(ind), 1)
        CALL GSWGSE_08 (xtmp2(ind),ytmp2(ind),ztmp2(ind),
     1               Xout(ind),Yout(ind),Zout(ind), 1)
      ELSE IF (trans=='GSEtoGEI') THEN
        CALL GSWGSE_08 (xtmp(ind),ytmp(ind),ztmp(ind),
     1               Xin(ind),Yin(ind),Zin(ind), -1)
        CALL GEOGSW_08 (xtmp2(ind),ytmp2(ind),ztmp2(ind),
     1               xtmp(ind),ytmp(ind),ztmp(ind), -1)
        CALL GEIGEO_08(Xout(ind),Yout(ind),Zout(ind),
     1               xtmp2(ind),ytmp2(ind),ztmp2(ind), -1)

      ELSE IF (trans=='MAGtoGSE') THEN
        CALL MAGSM_08 (Xin(ind),Yin(ind),Zin(ind),
     1               xtmp(ind),ytmp(ind),ztmp(ind), 1)
        CALL SMGSW_08 (xtmp(ind),ytmp(ind),ztmp(ind),
     1               xtmp2(ind),ytmp2(ind),ztmp2(ind), 1)
        CALL GSWGSE_08 (xtmp2(ind),ytmp2(ind),ztmp2(ind),
     1               Xout(ind),Yout(ind),Zout(ind), 1)

      ELSE IF (trans=='GSEtoMAG') THEN
        CALL GSWGSE_08 (xtmp(ind),ytmp(ind),ztmp(ind),
     1               Xin(ind),Yin(ind),Zin(ind), -1)
        CALL SMGSW_08 (xtmp2(ind),ytmp2(ind),ztmp2(ind),
     1               xtmp(ind),ytmp(ind),ztmp(ind), -1)
        CALL MAGSM_08 (Xout(ind),Yout(ind),Zout(ind),
     1               xtmp2(ind),ytmp2(ind),ztmp2(ind), -1)

      ELSE
       PRINT *, "INCORRECT VALUE: ", trans, "FOR trans"

      END IF
C10    CONTINUE
      END DO
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCC wrappers implemented by ANGEL GUTARRA-LEON CCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

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

      SUBROUTINE transform (indatav, trans, datetime, Nout,
     1                       Nv, Nt, mod_psi, cpsi, outdatav )
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
      INTEGER*4 :: Nv, Nt, Nout, i, datetime(Nt,5), mod_psi
      REAL*8 :: SPS, CPS, PSI, cpsi,
     *  indatav(Nv,3), outdatav(Nout,3),
     *  xtmp, ytmp, ztmp, xtmp2, ytmp2, ztmp2
      CHARACTER(10) :: trans
!f2py intent(in) ::  indatav, trans, datetime, Nout
!f2py optional, intent(in):: mod_psi = 1, cpsi=0
!f2py intent(hide) :: Nv, Nt
!f2py intent(out) ::  outdatav

c these transformations only rely on one call.
      IF      (trans=='GEItoGEO') THEN
        IF (Nv.EQ.Nt) THEN
        DO i=1,Nv
        CALL RECALC_08_W (datetime(i,:))
        CALL GEIGEO_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
        END DO
      ELSE IF (Nv.EQ.1) THEN
        DO i=1,Nt
          CALL RECALC_08_W (datetime(i,:))
          CALL GEIGEO_08(indatav(1,1),indatav(1,2),indatav(1,3),
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
        END DO
      ELSE
        CALL RECALC_08_W (datetime(1,:))
        DO i=1,Nv
          CALL GEIGEO_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
        END DO
      END IF


      ELSE IF (trans=='GEOtoGEI') THEN
        IF (Nv.EQ.Nt) THEN
          DO i=1,Nv
          CALL RECALC_08_W (datetime(i,:))
          CALL GEIGEO_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          END DO
        ELSE IF (Nv.EQ.1) THEN
          DO i=1,Nt
          CALL RECALC_08_W (datetime(i,:))
          CALL GEIGEO_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   indatav(1,1),indatav(1,2),indatav(1,3),-1)
          END DO
        ELSE
          CALL RECALC_08_W (datetime(1,:))
          DO i=1,Nv
          CALL GEIGEO_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          END DO
        END IF

C      ELSE IF (trans=='GEOtoGSM') THEN
C        CALL GEOGSW_08 (Xin(ind),Yin(ind),Zin(ind),
C     1               Xout(ind),Yout(ind),Zout(ind), 1)
C      ELSE IF (trans=='GSMtoGEO') THEN
C         CALL GEOGSW_08 (Xout(ind),Yout(ind),Zout(ind),
C     1               Xin(ind),Yin(ind),Zin(ind), -1)
C
C      ELSE IF (trans=='GEOtoMAG') THEN
C        CALL GEOMAG_08 (Xin(ind),Yin(ind),Zin(ind),
C     1               Xout(ind),Yout(ind),Zout(ind), 1)
C      ELSE IF (trans=='MAGtoGEO') THEN
C        CALL GEOMAG_08 (Xout(ind),Yout(ind),Zout(ind),
C     1               Xin(ind),Yin(ind),Zin(ind), -1)
C
C      ELSE IF (trans=='GSMtoGSE') THEN
C        CALL GSWGSE_08 (Xin(ind),Yin(ind),Zin(ind),
C     1               Xout(ind),Yout(ind),Zout(ind), 1)
C      ELSE IF (trans=='GSEtoGSM') THEN
C        CALL GSWGSE_08 (Xout(ind),Yout(ind),Zout(ind),
C     1               Xin(ind),Yin(ind),Zin(ind), -1)
C
C      ELSE  IF (trans=='MAGtoSM') THEN
C        CALL MAGSM_08 (Xin(ind),Yin(ind),Zin(ind),
C     1               Xout(ind),Yout(ind),Zout(ind), 1)
C      ELSE IF (trans=='SMtoMAG') THEN
C        CALL MAGSM_08 (Xout(ind),Yout(ind),Zout(ind),
C     1               Xin(ind),Yin(ind),Zin(ind), -1)
C
C      ELSE IF (trans=='SMtoGSM') THEN
C        CALL SMGSW_08 (Xin(ind),Yin(ind),Zin(ind),
C     1                 Xout(ind),Yout(ind),Zout(ind), 1)
C      ELSE IF (trans=='GSMtoSM') THEN
C        CALL SMGSW_08 (Xout(ind),Yout(ind),Zout(ind),
C     1               Xin(ind),Yin(ind),Zin(ind), -1)
C
C these transformations rely on multiple calls.
C      ELSE IF (trans=='MAGtoGSM') THEN
C        CALL MAGSM_08 (Xin(ind),Yin(ind),Zin(ind),
C     1               xtmp(ind),ytmp(ind),ztmp(ind), 1)
C        CALL SMGSW_08 (xtmp(ind),ytmp(ind),ztmp(ind),
C     1               Xout(ind),Yout(ind),Zout(ind), 1)
C      ELSE IF (trans=='GSMtoMAG') THEN
C        CALL SMGSW_08 (xtmp(ind),ytmp(ind),ztmp(ind),
C     1               Xin(ind),Yin(ind),Zin(ind), -1)
C        CALL MAGSM_08 (Xout(ind),Yout(ind),Zout(ind),
C     1               xtmp(ind),ytmp(ind),ztmp(ind), -1)
C
C
C      ELSE IF (trans=='GEItoMAG') THEN
C        CALL GEIGEO_08 (Xin(ind),Yin(ind),Zin(ind),
C     1               xtmp(ind),ytmp(ind),ztmp(ind), 1)
C        CALL GEOMAG_08 (xtmp(ind),ytmp(ind),ztmp(ind),
C     1               Xout(ind),Yout(ind),Zout(ind), 1)
C      ELSE IF (trans=='MAGtoGEI') THEN
C        CALL GEOMAG_08 (xtmp(ind),ytmp(ind),ztmp(ind),
C     1               Xin(ind),Yin(ind),Zin(ind), -1)
C        CALL GEIGEO_08 (Xout(ind),Yout(ind),Zout(ind),
C     1               xtmp(ind),ytmp(ind),ztmp(ind), -1)
C
C
C      ELSE IF (trans=='GEOtoSM') THEN
C        CALL GEOMAG_08 (Xin(ind),Yin(ind),Zin(ind),
C     1               xtmp(ind),ytmp(ind),ztmp(ind), 1)
C        CALL MAGSM_08 (xtmp(ind),ytmp(ind),ztmp(ind),
C     1               Xout(ind),Yout(ind),Zout(ind), 1)
C      ELSE IF (trans=='SMtoGEO') THEN
C        CALL MAGSM_08 (xtmp(ind),ytmp(ind),ztmp(ind),
C     1               Xin(ind),Yin(ind),Zin(ind), -1)
C        CALL GEOMAG_08 (Xout(ind),Yout(ind),Zout(ind),
C     1                xtmp(ind),ytmp(ind),ztmp(ind),-1)
C
C
C      ELSE IF (trans=='SMtoGSE') THEN
C        CALL SMGSW_08 (Xin(ind),Yin(ind),Zin(ind),
C     1               xtmp(ind),ytmp(ind),ztmp(ind), 1)
C        CALL GSWGSE_08 (xtmp(ind),ytmp(ind),ztmp(ind),
C     1               Xout(ind),Yout(ind),Zout(ind), 1)
C      ELSE IF (trans=='GSEtoSM') THEN
C        CALL GSWGSE_08 (xtmp(ind),ytmp(ind),ztmp(ind),
C     1               Xin(ind),Yin(ind),Zin(ind), -1)
C        CALL SMGSW_08 (Xout(ind),Yout(ind),Zout(ind),
C     1               xtmp(ind),ytmp(ind),ztmp(ind), -1)
C
C
C      ELSE IF (trans=='GEOtoGSE') THEN
C        CALL GEOGSW_08 (Xin(ind),Yin(ind),Zin(ind),
C     1               xtmp(ind),ytmp(ind),ztmp(ind), 1)
C        CALL GSWGSE_08 (xtmp(ind),ytmp(ind),ztmp(ind),
C     1               Xout(ind),Yout(ind),Zout(ind), 1)
C      ELSE IF (trans=='GSEtoGEO') THEN
C        CALL GSWGSE_08 (xtmp(ind),ytmp(ind),ztmp(ind),
C     1               Xin(ind),Yin(ind),Zin(ind), -1)
C        CALL GEOGSW_08 (Xout(ind),Yout(ind),Zout(ind),
C     1               xtmp(ind),ytmp(ind),ztmp(ind), -1)
C
C
C      ELSE IF (trans=='GEItoGSM') THEN
C        CALL GEIGEO_08 (Xin(ind),Yin(ind),Zin(ind),
C     1               xtmp(ind),ytmp(ind),ztmp(ind), 1)
C        CALL GEOGSW_08 (xtmp(ind),ytmp(ind),ztmp(ind),
C     1               Xout(ind),Yout(ind),Zout(ind), 1)
C      ELSE IF (trans=='GSMtoGEI') THEN
C        CALL GEOGSW_08 (xtmp(ind),ytmp(ind),ztmp(ind),
C     1               Xin(ind),Yin(ind),Zin(ind), -1)
C        CALL GEIGEO_08 (Xout(ind),Yout(ind),Zout(ind),
C     1               xtmp(ind),ytmp(ind),ztmp(ind), -1)
C
C
C      ELSE IF (trans=='GEItoSM') THEN
C        CALL GEIGEO_08 (Xin(ind),Yin(ind),Zin(ind),
C     1               xtmp(ind),ytmp(ind),ztmp(ind), 1)
C        CALL GEOMAG_08 (xtmp(ind),ytmp(ind),ztmp(ind),
C     1               xtmp2(ind),ytmp2(ind),ztmp2(ind), 1)
C        CALL MAGSM_08 (xtmp2(ind),ytmp2(ind),ztmp2(ind),
C     1               Xout(ind),Yout(ind),Zout(ind), 1)
C      ELSE IF (trans=='SMtoGEI') THEN
C        CALL MAGSM_08 (xtmp(ind),ytmp(ind),ztmp(ind),
C     1               Xin(ind),Yin(ind),Zin(ind), -1)
C        CALL GEOMAG_08 (xtmp2(ind),ytmp2(ind),ztmp2(ind),
C     1               xtmp(ind),ytmp(ind),ztmp(ind), -1)
C        CALL GEIGEO_08 (Xout(ind),Yout(ind),Zout(ind),
C     1               xtmp2(ind),ytmp2(ind),ztmp2(ind), -1)
C
C
C      ELSE IF (trans=='GEItoGSE') THEN
C        CALL GEIGEO_08 (Xin(ind),Yin(ind),Zin(ind),
C     1              xtmp(ind),ytmp(ind),ztmp(ind), 1)
C        CALL GEOGSW_08 (xtmp(ind),ytmp(ind),ztmp(ind),
C     1               xtmp2(ind),ytmp2(ind),ztmp2(ind), 1)
C        CALL GSWGSE_08 (xtmp2(ind),ytmp2(ind),ztmp2(ind),
C     1               Xout(ind),Yout(ind),Zout(ind), 1)
C      ELSE IF (trans=='GSEtoGEI') THEN
C        CALL GSWGSE_08 (xtmp(ind),ytmp(ind),ztmp(ind),
C     1               Xin(ind),Yin(ind),Zin(ind), -1)
C        CALL GEOGSW_08 (xtmp2(ind),ytmp2(ind),ztmp2(ind),
C     1               xtmp(ind),ytmp(ind),ztmp(ind), -1)
C        CALL GEIGEO_08(Xout(ind),Yout(ind),Zout(ind),
C     1               xtmp2(ind),ytmp2(ind),ztmp2(ind), -1)
C
C      ELSE IF (trans=='MAGtoGSE') THEN
C        CALL MAGSM_08 (Xin(ind),Yin(ind),Zin(ind),
C     1               xtmp(ind),ytmp(ind),ztmp(ind), 1)
C        CALL SMGSW_08 (xtmp(ind),ytmp(ind),ztmp(ind),
C     1               xtmp2(ind),ytmp2(ind),ztmp2(ind), 1)
C        CALL GSWGSE_08 (xtmp2(ind),ytmp2(ind),ztmp2(ind),
C     1               Xout(ind),Yout(ind),Zout(ind), 1)
C
C      ELSE IF (trans=='GSEtoMAG') THEN
C        CALL GSWGSE_08 (xtmp(ind),ytmp(ind),ztmp(ind),
C     1               Xin(ind),Yin(ind),Zin(ind), -1)
C        CALL SMGSW_08 (xtmp2(ind),ytmp2(ind),ztmp2(ind),
C     1               xtmp(ind),ytmp(ind),ztmp(ind), -1)
C        CALL MAGSM_08 (Xout(ind),Yout(ind),Zout(ind),
C     1               xtmp2(ind),ytmp2(ind),ztmp2(ind), -1)
C
C      ELSE
C       PRINT *, "INCORRECT VALUE: ", trans, "FOR trans"
C
      END IF

      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

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
      IF (trans=='GEItoGEO') THEN
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

      ELSE IF (trans=='GEOtoGSM') THEN
        IF (Nv.EQ.Nt) THEN
          DO i=1,Nv
          CALL RECALC_08_W (datetime(i,:))
          CALL GEOGSW_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        ELSE IF (Nv.EQ.1) THEN
          DO i=1,Nt
          CALL RECALC_08_W (datetime(i,:))
          CALL GEOGSW_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        ELSE
          CALL RECALC_08_W (datetime(1,:))
          DO i=1,Nv
          CALL GEOGSW_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        END IF


      ELSE IF (trans=='GSMtoGEO') THEN
        IF (Nv.EQ.Nt) THEN
          DO i=1,Nv
          CALL RECALC_08_W (datetime(i,:))
          CALL GEOGSW_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          END DO
        ELSE IF (Nv.EQ.1) THEN
          DO i=1,Nt
          CALL RECALC_08_W (datetime(i,:))
          CALL GEOGSW_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   indatav(1,1),indatav(1,2),indatav(1,3),-1)
          END DO
        ELSE
          CALL RECALC_08_W (datetime(1,:))
          DO i=1,Nv
          CALL GEOGSW_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          END DO
        END IF

      ELSE IF (trans=='SMtoGSM') THEN
        IF (Nv.EQ.Nt) THEN
          DO i=1,Nv
          CALL RECALC_08_W (datetime(i,:))
          CALL SMGSW_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        ELSE IF (Nv.EQ.1) THEN
          DO i=1,Nt
          CALL RECALC_08_W (datetime(i,:))
          CALL SMGSW_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        ELSE
          CALL RECALC_08_W (datetime(1,:))
          DO i=1,Nv
          CALL SMGSW_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        END IF


      ELSE IF (trans=='GSMtoSM') THEN
        IF (Nv.EQ.Nt) THEN
          DO i=1,Nv
          CALL RECALC_08_W (datetime(i,:))
          CALL SMGSW_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          END DO
        ELSE IF (Nv.EQ.1) THEN
          DO i=1,Nt
          CALL RECALC_08_W (datetime(i,:))
          CALL SMGSW_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   indatav(1,1),indatav(1,2),indatav(1,3),-1)
          END DO
        ELSE
          CALL RECALC_08_W (datetime(1,:))
          DO i=1,Nv
          CALL SMGSW_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          END DO
        END IF




      ELSE IF (trans=='GEOtoMAG') THEN
        IF (Nv.EQ.Nt) THEN
          DO i=1,Nv
          CALL RECALC_08_W (datetime(i,:))
          CALL GEOMAG_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        ELSE IF (Nv.EQ.1) THEN
          DO i=1,Nt
          CALL RECALC_08_W (datetime(i,:))
          CALL GEOMAG_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        ELSE
          CALL RECALC_08_W (datetime(1,:))
          DO i=1,Nv
          CALL GEOMAG_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        END IF


      ELSE IF (trans=='MAGtoGEO') THEN
        IF (Nv.EQ.Nt) THEN
          DO i=1,Nv
          CALL RECALC_08_W (datetime(i,:))
          CALL GEOMAG_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          END DO
        ELSE IF (Nv.EQ.1) THEN
          DO i=1,Nt
          CALL RECALC_08_W (datetime(i,:))
          CALL GEOMAG_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   indatav(1,1),indatav(1,2),indatav(1,3),-1)
          END DO
        ELSE
          CALL RECALC_08_W (datetime(1,:))
          DO i=1,Nv
          CALL GEOMAG_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          END DO
        END IF


      ELSE IF (trans=='GSMtoGSE') THEN
        IF (Nv.EQ.Nt) THEN
          DO i=1,Nv
          CALL RECALC_08_W (datetime(i,:))
          CALL GSWGSE_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        ELSE IF (Nv.EQ.1) THEN
          DO i=1,Nt
          CALL RECALC_08_W (datetime(i,:))
          CALL GSWGSE_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        ELSE
          CALL RECALC_08_W (datetime(1,:))
          DO i=1,Nv
          CALL GSWGSE_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        END IF

      ELSE IF (trans=='GSEtoGSM') THEN
        IF (Nv.EQ.Nt) THEN
          DO i=1,Nv
          CALL RECALC_08_W (datetime(i,:))
          CALL GSWGSE_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          END DO
        ELSE IF (Nv.EQ.1) THEN
          DO i=1,Nt
          CALL RECALC_08_W (datetime(i,:))
          CALL GSWGSE_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   indatav(1,1),indatav(1,2),indatav(1,3),-1)
          END DO
        ELSE
          CALL RECALC_08_W (datetime(1,:))
          DO i=1,Nv
          CALL GSWGSE_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          END DO
        END IF



      ELSE IF (trans=='MAGtoSM') THEN
        IF (Nv.EQ.Nt) THEN
          DO i=1,Nv
          CALL RECALC_08_W (datetime(i,:))
          CALL MAGSM_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        ELSE IF (Nv.EQ.1) THEN
          DO i=1,Nt
          CALL RECALC_08_W (datetime(i,:))
          CALL MAGSM_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        ELSE
          CALL RECALC_08_W (datetime(1,:))
          DO i=1,Nv
          CALL MAGSM_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        END IF

      ELSE IF (trans=='SMtoMAG') THEN
        IF (Nv.EQ.Nt) THEN
          DO i=1,Nv
          CALL RECALC_08_W (datetime(i,:))
          CALL MAGSM_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          END DO
        ELSE IF (Nv.EQ.1) THEN
          DO i=1,Nt
          CALL RECALC_08_W (datetime(i,:))
          CALL MAGSM_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   indatav(1,1),indatav(1,2),indatav(1,3),-1)
          END DO
        ELSE
          CALL RECALC_08_W (datetime(1,:))
          DO i=1,Nv
          CALL MAGSM_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          END DO
        END IF


C these transformations rely on multiple calls.

      ELSE IF (trans=='MAGtoGSM') THEN
        IF (Nv.EQ.Nt) THEN
          DO i=1,Nv
          CALL RECALC_08_W (datetime(i,:))
          CALL MAGSM_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   xtmp, ytmp, ztmp, 1)
          CALL SMGSW_08(xtmp, ytmp, ztmp,
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        ELSE IF (Nv.EQ.1) THEN
          DO i=1,Nt
          CALL RECALC_08_W (datetime(i,:))
          CALL MAGSM_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   xtmp, ytmp, ztmp, 1)
          CALL SMGSW_08(xtmp, ytmp, ztmp,
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        ELSE
          CALL RECALC_08_W (datetime(1,:))
          DO i=1,Nv
          CALL MAGSM_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   xtmp, ytmp, ztmp, 1)
          CALL SMGSW_08(xtmp, ytmp, ztmp,
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        END IF

      ELSE IF (trans=='GSMtoMAG') THEN
        IF (Nv.EQ.Nt) THEN
          DO i=1,Nv
          CALL RECALC_08_W (datetime(i,:))
          CALL SMGSW_08(xtmp, ytmp, ztmp,
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          CALL MAGSM_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   xtmp, ytmp, ztmp,-1)
          END DO
        ELSE IF (Nv.EQ.1) THEN
          DO i=1,Nt
          CALL SMGSW_08(xtmp, ytmp, ztmp,
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          CALL MAGSM_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   xtmp, ytmp, ztmp,-1)
          END DO
        ELSE
          CALL RECALC_08_W (datetime(1,:))
          DO i=1,Nv
          CALL SMGSW_08(xtmp, ytmp, ztmp,
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          CALL MAGSM_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   xtmp, ytmp, ztmp,-1)
          END DO
        END IF



      ELSE IF (trans=='GEItoMAG') THEN
        IF (Nv.EQ.Nt) THEN
          DO i=1,Nv
          CALL RECALC_08_W (datetime(i,:))
          CALL GEIGEO_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   xtmp, ytmp, ztmp, 1)
          CALL GEOMAG_08(xtmp, ytmp, ztmp,
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        ELSE IF (Nv.EQ.1) THEN
          DO i=1,Nt
          CALL RECALC_08_W (datetime(i,:))
          CALL GEIGEO_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   xtmp, ytmp, ztmp, 1)
          CALL GEOMAG_08(xtmp, ytmp, ztmp,
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        ELSE
          CALL RECALC_08_W (datetime(1,:))
          DO i=1,Nv
          CALL GEIGEO_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   xtmp, ytmp, ztmp, 1)
          CALL GEOMAG_08(xtmp, ytmp, ztmp,
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        END IF
      ELSE IF (trans=='MAGtoGEI') THEN
        IF (Nv.EQ.Nt) THEN
          DO i=1,Nv
          CALL RECALC_08_W (datetime(i,:))
          CALL GEOMAG_08(xtmp, ytmp, ztmp,
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          CALL GEIGEO_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   xtmp, ytmp, ztmp,-1)
          END DO
        ELSE IF (Nv.EQ.1) THEN
          DO i=1,Nt
          CALL RECALC_08_W (datetime(i,:))
          CALL GEOMAG_08(xtmp, ytmp, ztmp,
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          CALL GEIGEO_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   xtmp, ytmp, ztmp,-1)
          END DO
        ELSE
          CALL RECALC_08_W (datetime(1,:))
          DO i=1,Nv
          CALL GEOMAG_08(xtmp, ytmp, ztmp,
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          CALL GEIGEO_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   xtmp, ytmp, ztmp,-1)
          END DO
        END IF


      ELSE IF (trans=='GEOtoSM') THEN
        IF (Nv.EQ.Nt) THEN
          DO i=1,Nv
          CALL RECALC_08_W (datetime(i,:))
          CALL GEOMAG_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   xtmp, ytmp, ztmp, 1)
          CALL MAGSM_08(xtmp, ytmp, ztmp,
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        ELSE IF (Nv.EQ.1) THEN
          DO i=1,Nt
          CALL RECALC_08_W (datetime(i,:))
          CALL GEOMAG_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   xtmp, ytmp, ztmp, 1)
          CALL MAGSM_08(xtmp, ytmp, ztmp,
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        ELSE
          CALL RECALC_08_W (datetime(1,:))
          DO i=1,Nv
          CALL GEIGEO_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   xtmp, ytmp, ztmp, 1)
          CALL MAGSM_08(xtmp, ytmp, ztmp,
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        END IF
      ELSE IF (trans=='SMtoGEO') THEN
        IF (Nv.EQ.Nt) THEN
          DO i=1,Nv
          CALL RECALC_08_W (datetime(i,:))
          CALL MAGSM_08(xtmp, ytmp, ztmp,
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          CALL GEOMAG_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   xtmp, ytmp, ztmp,-1)
          END DO
        ELSE IF (Nv.EQ.1) THEN
          DO i=1,Nt
          CALL RECALC_08_W (datetime(i,:))
          CALL MAGSM_08(xtmp, ytmp, ztmp,
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          CALL GEOMAG_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   xtmp, ytmp, ztmp,-1)
          END DO
        ELSE
          CALL RECALC_08_W (datetime(1,:))
          DO i=1,Nv
          CALL MAGSM_08(xtmp, ytmp, ztmp,
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          CALL GEOMAG_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   xtmp, ytmp, ztmp,-1)
          END DO
        END IF




      ELSE IF (trans=='SMtoGSE') THEN
        IF (Nv.EQ.Nt) THEN
          DO i=1,Nv
          CALL RECALC_08_W (datetime(i,:))
          CALL SMGSW_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   xtmp, ytmp, ztmp, 1)
          CALL GSWGSE_08(xtmp, ytmp, ztmp,
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        ELSE IF (Nv.EQ.1) THEN
          DO i=1,Nt
          CALL RECALC_08_W (datetime(i,:))
          CALL SMGSW_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   xtmp, ytmp, ztmp, 1)
          CALL GSWGSE_08(xtmp, ytmp, ztmp,
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        ELSE
          CALL RECALC_08_W (datetime(1,:))
          DO i=1,Nv
          CALL SMGSW_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   xtmp, ytmp, ztmp, 1)
          CALL GSWGSE_08(xtmp, ytmp, ztmp,
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        END IF
      ELSE IF (trans=='GSEtoSM') THEN
        IF (Nv.EQ.Nt) THEN
          DO i=1,Nv
          CALL RECALC_08_W (datetime(i,:))
          CALL GSWGSE_08(xtmp, ytmp, ztmp,
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          CALL SMGSW_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   xtmp, ytmp, ztmp,-1)
          END DO
        ELSE IF (Nv.EQ.1) THEN
          DO i=1,Nt
          CALL RECALC_08_W (datetime(i,:))
          CALL GSWGSE_08(xtmp, ytmp, ztmp,
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          CALL SMGSW_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   xtmp, ytmp, ztmp,-1)
          END DO
        ELSE
          CALL RECALC_08_W (datetime(1,:))
          DO i=1,Nv
          CALL GSWGSE_08(xtmp, ytmp, ztmp,
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          CALL SMGSW_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   xtmp, ytmp, ztmp,-1)
          END DO
        END IF




      ELSE IF (trans=='GEOtoGSE') THEN
        IF (Nv.EQ.Nt) THEN
          DO i=1,Nv
          CALL RECALC_08_W (datetime(i,:))
          CALL GEOGSW_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   xtmp, ytmp, ztmp, 1)
          CALL GSWGSE_08(xtmp, ytmp, ztmp,
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        ELSE IF (Nv.EQ.1) THEN
          DO i=1,Nt
          CALL RECALC_08_W (datetime(i,:))
          CALL GEOGSW_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   xtmp, ytmp, ztmp, 1)
          CALL GSWGSE_08(xtmp, ytmp, ztmp,
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        ELSE
          CALL RECALC_08_W (datetime(1,:))
          DO i=1,Nv
          CALL GEOGSW_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   xtmp, ytmp, ztmp, 1)
          CALL GSWGSE_08(xtmp, ytmp, ztmp,
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        END IF
      ELSE IF (trans=='GSEtoGEO') THEN
        IF (Nv.EQ.Nt) THEN
          DO i=1,Nv
          CALL RECALC_08_W (datetime(i,:))
          CALL GSWGSE_08(xtmp, ytmp, ztmp,
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          CALL GEOGSW_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   xtmp, ytmp, ztmp,-1)
          END DO
        ELSE IF (Nv.EQ.1) THEN
          DO i=1,Nt
          CALL RECALC_08_W (datetime(i,:))
          CALL GSWGSE_08(xtmp, ytmp, ztmp,
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          CALL GEOGSW_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   xtmp, ytmp, ztmp,-1)
          END DO
        ELSE
          CALL RECALC_08_W (datetime(1,:))
          DO i=1,Nv
          CALL GSWGSE_08(xtmp, ytmp, ztmp,
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          CALL GEOGSW_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   xtmp, ytmp, ztmp,-1)
          END DO
        END IF





      ELSE IF (trans=='GEItoGSM') THEN
        IF (Nv.EQ.Nt) THEN
          DO i=1,Nv
          CALL RECALC_08_W (datetime(i,:))
          CALL GEIGEO_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   xtmp, ytmp, ztmp, 1)
          CALL GEOGSW_08(xtmp, ytmp, ztmp,
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        ELSE IF (Nv.EQ.1) THEN
          DO i=1,Nt
          CALL RECALC_08_W (datetime(i,:))
          CALL GEIGEO_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   xtmp, ytmp, ztmp, 1)
          CALL GEOGSW_08(xtmp, ytmp, ztmp,
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        ELSE
          CALL RECALC_08_W (datetime(1,:))
          DO i=1,Nv
          CALL GEIGEO_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   xtmp, ytmp, ztmp, 1)
          CALL GEOGSW_08(xtmp, ytmp, ztmp,
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        END IF
      ELSE IF (trans=='GSMtoGEI') THEN
        IF (Nv.EQ.Nt) THEN
          DO i=1,Nv
          CALL RECALC_08_W (datetime(i,:))
          CALL GEOGSW_08(xtmp, ytmp, ztmp,
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          CALL GEIGEO_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   xtmp, ytmp, ztmp,-1)
          END DO
        ELSE IF (Nv.EQ.1) THEN
          DO i=1,Nt
          CALL RECALC_08_W (datetime(i,:))
          CALL GEOGSW_08(xtmp, ytmp, ztmp,
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          CALL GEIGEO_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   xtmp, ytmp, ztmp,-1)
          END DO
        ELSE
          CALL RECALC_08_W (datetime(1,:))
          DO i=1,Nv
          CALL GEOGSW_08(xtmp, ytmp, ztmp,
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          CALL GEIGEO_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   xtmp, ytmp, ztmp,-1)
          END DO
        END IF

C start of three subroutine calls for transformations

      ELSE IF (trans=='GEItoSM') THEN
        IF (Nv.EQ.Nt) THEN
          DO i=1,Nv
          CALL RECALC_08_W (datetime(i,:))
          CALL GEIGEO_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   xtmp, ytmp, ztmp, 1)
          CALL GEOMAG_08(xtmp, ytmp, ztmp,
     *                   xtmp2, ytmp2, ztmp2,1)
          CALL MAGSM_08(xtmp2, ytmp2, ztmp2,
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        ELSE IF (Nv.EQ.1) THEN
          DO i=1,Nt
          CALL GEIGEO_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   xtmp, ytmp, ztmp, 1)
          CALL GEOMAG_08(xtmp, ytmp, ztmp,
     *                   xtmp2, ytmp2, ztmp2,1)
          CALL MAGSM_08(xtmp2, ytmp2, ztmp2,
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        ELSE
          CALL RECALC_08_W (datetime(1,:))
          DO i=1,Nv
          CALL GEIGEO_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   xtmp, ytmp, ztmp, 1)
          CALL GEOMAG_08(xtmp, ytmp, ztmp,
     *                   xtmp2, ytmp2, ztmp2,1)
          CALL MAGSM_08(xtmp2, ytmp2, ztmp2,
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        END IF
      ELSE IF (trans=='SMtoGEI') THEN
        IF (Nv.EQ.Nt) THEN
          DO i=1,Nv
          CALL RECALC_08_W (datetime(i,:))
          CALL MAGSM_08(xtmp, ytmp, ztmp,
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          CALL GEOMAG_08(xtmp2, ytmp2, ztmp2,
     *                   xtmp, ytmp, ztmp,-1)
          CALL GEIGEO_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   xtmp2, ytmp2, ztmp2,-1)
          END DO
        ELSE IF (Nv.EQ.1) THEN
          DO i=1,Nt
          CALL RECALC_08_W (datetime(i,:))
          CALL MAGSM_08(xtmp, ytmp, ztmp,
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          CALL GEOMAG_08(xtmp2, ytmp2, ztmp2,
     *                   xtmp, ytmp, ztmp,-1)
          CALL GEIGEO_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   xtmp2, ytmp2, ztmp2,-1)
          END DO
        ELSE
          CALL RECALC_08_W (datetime(1,:))
          DO i=1,Nv
          CALL MAGSM_08(xtmp, ytmp, ztmp,
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          CALL GEOMAG_08(xtmp2, ytmp2, ztmp2,
     *                   xtmp, ytmp, ztmp,-1)
          CALL GEIGEO_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   xtmp2, ytmp2, ztmp2,-1)
          END DO
        END IF




      ELSE IF (trans=='GEItoGSE') THEN
        IF (Nv.EQ.Nt) THEN
          DO i=1,Nv
          CALL RECALC_08_W (datetime(i,:))
          CALL GEIGEO_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   xtmp, ytmp, ztmp, 1)
          CALL GEOGSW_08(xtmp, ytmp, ztmp,
     *                   xtmp2, ytmp2, ztmp2,1)
          CALL GSWGSE_08(xtmp2, ytmp2, ztmp2,
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        ELSE IF (Nv.EQ.1) THEN
          DO i=1,Nt
          CALL GEIGEO_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   xtmp, ytmp, ztmp, 1)
          CALL GEOGSW_08(xtmp, ytmp, ztmp,
     *                   xtmp2, ytmp2, ztmp2,1)
          CALL GSWGSE_08(xtmp2, ytmp2, ztmp2,
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        ELSE
          CALL RECALC_08_W (datetime(1,:))
          DO i=1,Nv
          CALL GEIGEO_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   xtmp, ytmp, ztmp, 1)
          CALL GEOGSW_08(xtmp, ytmp, ztmp,
     *                   xtmp2, ytmp2, ztmp2,1)
          CALL GSWGSE_08(xtmp2, ytmp2, ztmp2,
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        END IF
      ELSE IF (trans=='GSEtoGEI') THEN
        IF (Nv.EQ.Nt) THEN
          DO i=1,Nv
          CALL RECALC_08_W (datetime(i,:))
          CALL GSWGSE_08(xtmp, ytmp, ztmp,
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          CALL GEOGSW_08(xtmp2, ytmp2, ztmp2,
     *                   xtmp, ytmp, ztmp,-1)
          CALL GEIGEO_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   xtmp2, ytmp2, ztmp2,-1)
          END DO
        ELSE IF (Nv.EQ.1) THEN
          DO i=1,Nt
          CALL RECALC_08_W (datetime(i,:))
          CALL GSWGSE_08(xtmp, ytmp, ztmp,
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          CALL GEOGSW_08(xtmp2, ytmp2, ztmp2,
     *                   xtmp, ytmp, ztmp,-1)
          CALL GEIGEO_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   xtmp2, ytmp2, ztmp2,-1)
          END DO
        ELSE
          CALL RECALC_08_W (datetime(1,:))
          DO i=1,Nv
          CALL GSWGSE_08(xtmp, ytmp, ztmp,
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          CALL GEOGSW_08(xtmp2, ytmp2, ztmp2,
     *                   xtmp, ytmp, ztmp,-1)
          CALL GEIGEO_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   xtmp2, ytmp2, ztmp2,-1)
          END DO
        END IF



      ELSE IF (trans=='MAGtoGSE') THEN
        IF (Nv.EQ.Nt) THEN
          DO i=1,Nv
          CALL RECALC_08_W (datetime(i,:))
          CALL MAGSM_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   xtmp, ytmp, ztmp, 1)
          CALL SMGSW_08(xtmp, ytmp, ztmp,
     *                   xtmp2, ytmp2, ztmp2,1)
          CALL GSWGSE_08(xtmp2, ytmp2, ztmp2,
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        ELSE IF (Nv.EQ.1) THEN
          DO i=1,Nt
          CALL MAGSM_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   xtmp, ytmp, ztmp, 1)
          CALL SMGSW_08(xtmp, ytmp, ztmp,
     *                   xtmp2, ytmp2, ztmp2,1)
          CALL GSWGSE_08(xtmp2, ytmp2, ztmp2,
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        ELSE
          CALL RECALC_08_W (datetime(1,:))
          DO i=1,Nv
          CALL MAGSM_08(indatav(i,1),indatav(i,2),indatav(i,3),
     *                   xtmp, ytmp, ztmp, 1)
          CALL SMGSW_08(xtmp, ytmp, ztmp,
     *                   xtmp2, ytmp2, ztmp2,1)
          CALL GSWGSE_08(xtmp2, ytmp2, ztmp2,
     *                   outdatav(i,1),outdatav(i,2),outdatav(i,3),1)
          END DO
        END IF
      ELSE IF (trans=='GSEtoMAG') THEN
        IF (Nv.EQ.Nt) THEN
          DO i=1,Nv
          CALL RECALC_08_W (datetime(i,:))
          CALL GSWGSE_08(xtmp, ytmp, ztmp,
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          CALL SMGSW_08(xtmp2, ytmp2, ztmp2,
     *                   xtmp, ytmp, ztmp,-1)
          CALL MAGSM_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   xtmp2, ytmp2, ztmp2,-1)
          END DO
        ELSE IF (Nv.EQ.1) THEN
          DO i=1,Nt
          CALL RECALC_08_W (datetime(i,:))
          CALL GSWGSE_08(xtmp, ytmp, ztmp,
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          CALL SMGSW_08(xtmp2, ytmp2, ztmp2,
     *                   xtmp, ytmp, ztmp,-1)
          CALL MAGSM_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   xtmp2, ytmp2, ztmp2,-1)
          END DO
        ELSE
          CALL RECALC_08_W (datetime(1,:))
          DO i=1,Nv
          CALL GSWGSE_08(xtmp, ytmp, ztmp,
     *                   indatav(i,1),indatav(i,2),indatav(i,3),-1)
          CALL SMGSW_08(xtmp2, ytmp2, ztmp2,
     *                   xtmp, ytmp, ztmp,-1)
          CALL MAGSM_08(outdatav(i,1),outdatav(i,2),outdatav(i,3),
     *                   xtmp2, ytmp2, ztmp2,-1)
          END DO
        END IF

      END IF

      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

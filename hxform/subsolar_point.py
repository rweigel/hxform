def _nint(a):
  """Replicates Fortran's NINT behavior (rounds halfway cases away from zero)."""
  # https://stackoverflow.com/questions/45884588/how-does-fortran-convert-a-real-number-to-integer
  if a >= 0:
    return int(a + 0.5)
  else:
    return int(a - 0.5)


def _mead(IYR, IDAY, SECS, lib_transform='geopack_08_dp'):

  # TODO: Implement.
  # Russel 1971 Appendix 2.
  """
  G.D. Mead (private communication) has written a simple subroutine to calculate the position of the Sun in GEI coordinates. It
  is accurate for years 1901 through 2099, to within 0.006 deg. The input is the year, day of year and seconds of the day in UT.
  The output is Greenwich Mean Sideral Time in degrees, the ecliptic longitude, apparent right ascension and declination of the
  Sun in degrees. The listing of this program follows. We note that the cartesian coordinates of the vector from the Earth to the
  Sun are:
  X = cos(SRASN) cos(SDEC)
  Y = sin(SRASN) cos(SDEC)
  Z = sin(SDEC)
  SUBROUTINE SUN(IYR, IDAY, SECS, GST, SLONG, SRASN, SDEC)
  C PROGRAM TO CALCULATE SIDEREAL TIME AND POSITION OF THE SUN.
  C GOOD FOR YEARS 1901 THROUGH 2099. ACCURACY 0.006 DEGREE.
  C INPUT IS IYR, IDAY (INTEGERS), AND SECS, DEFINING UN. TIME.
  C OUTPUT IS GREENWICH MEAN SIDEREAL TIME (GST) IN DEGREES,
  C LONGITUDE ALONG ECLIPTIC (SLONG), AND APPARENT RIGHT ASCENSION
  C AND DECLINATION (SRASN, SDEC) OF THE SUN, ALL IN DEGREES
  """

  from math import pi, sin, cos, atan, tan, atan2, sqrt
  def cotan(x):
    return 1.0 / tan(x)
  def DMOD(x, y):
    return x - y * int(x / y)


  #DATA RAD /57.29578/
  RAD = 180.0 / pi

  #DOUBLE PRECISION DJ, FDAY
  #IF(IYR. LT. 1901. OR. IYR. GT. 2099) RETURN
  if IYR < 1901 or IYR > 2099:
    raise ValueError('1901 <= year <= 2099 required')

  #FDAY = SECS/86400
  FDAY = SECS/86400

  #DJ = 365* (IYR-1900) + (IYR-1901)/4 + IDAY + FDAY -0.5D0
  DJ = 365* (IYR-1900) + (IYR-1901)/4 + IDAY + FDAY -0.5

  #T = DJ / 36525
  T = DJ / 36525

  #VL = DMOD (279.696678 + 0.9856473354*DJ, 360.D0)
  VL = DMOD (279.696678 + 0.9856473354*DJ, 360.)

  #GST = DMOD (279.690983 + 0.9856473354*DJ + 360.*FDAY + 180., 360.D0)
  GST = DMOD (279.690983 + 0.9856473354*DJ + 360.*FDAY + 180., 360.)

  #G = DMOD (358.475845 + 0.985600267*DJ, 360.D0) / RAD
  G = DMOD (358.475845 + 0.985600267*DJ, 360.) / RAD

  #SLONG = VL + (1.91946 -0.004789*T)*SIN(G) + 0.020094*SIN (2.*G)
  SLONG = VL + (1.91946 -0.004789*T)*sin(G) + 0.020094*sin (2.*G)

  #OBLIQ = (23.45229 -0.0130125*T) / RAD
  OBLIQ = (23.45229 -0.0130125*T) / RAD

  #SLP = (SLONG -0.005686) / RAD
  SLP = (SLONG -0.005686) / RAD

  #SIND = SIN (OBLIQ)*SIN (SLP)
  SIND = sin (OBLIQ)*sin (SLP)

  #COSD = SQRT(1.-SIND**2)
  COSD = sqrt(1.-SIND**2)

  #SDEC = RAD * ATAN (SIND/COSD)
  SDEC = RAD * atan (SIND/COSD)

  #SRASN = 180. -RAD*ATAN2
  #(COTAN (OBLIQ)*SIND/COSD, -COS (SLP)/COSD)

  SRASN = 180. -RAD*atan2(cotan (OBLIQ)*SIND/COSD, -cos (SLP)/COSD)

  # Added
  sbsllon = SRASN - GST

  return SDEC, sbsllon


def _laundal(year, doy, ut):

  from math import floor, sin, cos, pi, atan2, asin

  # Appendix C of https://link.springer.com/article/10.1007/s11214-016-0275-y#appendices
  # Reference frame of output not stated in Appendix, but from paragraph that
  # references appendix, it appears to be GEO.
  yr = year - 2000
  nleap = floor((year - 1601)/4.)
  nleap = nleap - 99
  if year <= 1900:
    ncent = floor((year - 1601)/100.)
    ncent = 3 - ncent
    nleap = nleap + ncent
  l0 = -79.549+(-.238699*(yr-4*nleap)+ 3.08514e-2*nleap)
  g0 = -2.472+(-.2558905*(yr-4*nleap)- 3.79617e-2*nleap)
  df = (ut/86400. - 1.5) + doy
  lf = .9856474*df
  gf = .9856003*df
  l = l0 + lf
  g = g0 + gf
  grad = g*pi/180.
  lmbda = l + 1.915*sin(grad) + .020*sin(2.*grad)
  lmrad = lmbda*pi/180.
  sinlm = sin(lmrad)
  n = df + 365.*yr + nleap
  epsilon = 23.439 - 4.0e-7*n
  epsrad = epsilon*pi/180.
  alpha = atan2(cos(epsrad)*sinlm, cos(lmrad)) * 180./pi
  delta = asin(sin(epsrad)*sinlm) * 180./pi
  sbslcolat = 90 - delta
  etdeg = l - alpha
  nrot = round(etdeg/360.)
  etdeg = etdeg - 360.*nrot
  aptime = ut/240. + etdeg
  sbsllon = 180. - aptime
  nrot = round(sbsllon/360.)
  sbsllon = sbsllon - 360.*nrot

  return 90-sbslcolat, sbsllon


def _apex(iyr, iday, ihr, imn, sec):
  # Subroutine subsol in
  # https://github.com/NCAR/apex_fortran/blob/master/apex.f90
  # subsolar geographic latitude and longitude given the date and time (Universal Time)
  from math import sin, cos, atan2, asin

  dtr = 0.0174532925199432957692369076847
  rtd = 57.2957795130823208767981548147

  minyear = 1601
  maxyear = 2100

  sbsllat=0.
  sbsllon=0.

  yr = iyr-2000
  nleap = (iyr-1601)/4
  nleap = nleap - 99
  if (iyr <= 1900):
    if (iyr < minyear):
      raise ValueError('subsolr invalid before 1601:  input year = ',iyr)
    ncent = (iyr-minyear)/100
    ncent = 3 - ncent
    nleap = nleap + ncent
  if (iyr > maxyear):
    raise ValueError('subsolr invalid after 2100:  input year = ',iyr)

  l0 = -79.549 + (-.238699*(yr-4*nleap) + 3.08514e-2*nleap)
  g0 = -2.472 + (-.2558905*(yr-4*nleap) - 3.79617e-2*nleap)
  ut = float(ihr*3600 + imn*60) + sec
  df = (ut/86400. - 1.5) + iday
  lf = .9856474*df
  gf = .9856003*df
  l = l0 + lf
  g = g0 + gf
  grad = g*dtr
  lambda_ = l + 1.915*sin(grad) + .020*sin(2.*grad)
  lamrad = lambda_*dtr
  sinlam = sin(lamrad)
  n = df + 365.*yr + float(nleap)
  epsilon = 23.439 - 4.e-7*n
  epsrad = epsilon*dtr
  alpha = atan2(cos(epsrad)*sinlam,cos(lamrad))*rtd
  delta = asin(sin(epsrad)*sinlam)*rtd
  sbsllat = delta
  etdeg = l - alpha
  nrot = _nint(etdeg/360.)
  etdeg = etdeg - float(360*nrot)
  aptime = ut/240. + etdeg
  sbsllon = 180. - aptime
  nrot = _nint(sbsllon/360.)
  sbsllon = sbsllon - float(360*nrot)

  return sbsllat, sbsllon


def _tiegcm(IYR, IDAY, IHR, IMN, SEC):
  from math import sin, cos, atan2, asin
  # https://fpi.hao.ucar.edu/hao/aim3/marsal/tiegcm/maggrd/src_code/subsol.F
  # Find subsolar geographic latitude and longitude given the date and time (Universal Time)
  D2R = 0.0174532925199432957692369076847
  R2D = 57.2957795130823208767981548147
  YR = IYR - 2000
  NLEAP = (IYR-1601)/4
  NLEAP = NLEAP - 99
  if (IYR <= 1900):
    if (IYR <= 1600):
      raise ValueError('SUBSOLR INVALID BEFORE 1601')

    NCENT = (IYR-1601)/100
    NCENT = 3 - NCENT
    NLEAP = NLEAP + NCENT

  if (IYR >= 2101):
    raise ValueError('SUBSOLR INVALID AFTER 2100')

  L0 = -79.549 + (-.238699*(YR-4*NLEAP) + 3.08514E-2*NLEAP)
  G0 = -2.472 + (-.2558905*(YR-4*NLEAP) - 3.79617E-2*NLEAP)
  UT = float(IHR*3600 + IMN*60) + SEC
  DF = (UT/86400. - 1.5) + IDAY
  LF = .9856474*DF
  GF = .9856003*DF
  L = L0 + LF
  G = G0 + GF
  GRAD = G*D2R
  LAMBDA = L + 1.915*sin(GRAD) + .020*sin(2.*GRAD)
  LAMRAD = LAMBDA*D2R
  SINLAM = sin(LAMRAD)
  N = DF + float(365*YR + NLEAP)
  EPSILON = 23.439 - 4.E-7*N
  EPSRAD = EPSILON*D2R
  ALPHA = atan2(cos(EPSRAD)*SINLAM, cos(LAMRAD))*R2D
  DELTA = asin(sin(EPSRAD)*SINLAM)*R2D
  SBSLLAT = DELTA
  ETDEG = L - ALPHA
  NROT = _nint(ETDEG/360.)
  ETDEG = ETDEG - float(360*NROT)
  APTIME = UT/240. + ETDEG
  SBSLLON = 180. - APTIME
  NROT = _nint(SBSLLON/360.)
  SBSLLON = SBSLLON - float(360*NROT)

  return SBSLLAT, SBSLLON


def _hxform(time, frame='GEO', lib='geopack_08_dp'):
  """Compute the subsolar point at a given time in a specified frame.
  This is done by transforming [1, 0, 0] from GSM to the desired frame.
  """

  import numpy
  import hxform
  outer_type = type(time)

  if lib not in hxform.libs():
    raise ValueError(f"Library '{lib}' not in available libraries: {hxform.libs()}")

  subsol_pt = hxform.transform([1., 0., 0.], time, 'GSM', frame, lib=lib)

  if isinstance(time, outer_type):
    return outer_type(subsol_pt)
  else:
    return numpy.array(subsol_pt)


def subsolar_point(t, frame='GEO', lib='geopack_08_dp', lib_transform='geopack_08_dp'):

  # TODO: Also use https://docs.astropy.org/en/latest/api/astropy.coordinates.get_sun.html
  import numpy as np

  import hxform

  libs_alt = ['tiegcm', 'apex', 'laundal', 'mead']
  libs = [*libs_alt, *hxform.libs()]
  if lib not in libs:
    raise ValueError(f"Library '{lib}' not in available libraries: {libs}")

  if lib not in libs_alt:
    return _hxform(t, frame=frame, lib=lib)

  # rtp: radius, theta (colat), phi (lon)
  from hxform.transform import _t_prep
  if lib in libs_alt:
    outer_type = type(t)

    ts, _ = _t_prep(t)

    rtp = np.full((len(ts), 3), np.nan)
    rtp[:, 0] = 1.0  # radius

    for i, t in enumerate(ts):
      yr = t[0] # year
      doy = hxform.timelib.ymd2doy(t[0:3]) # day of year
      hms = t[3:6] # hours, minutes, seconds
      sod = t[3]*3600 + t[4]*60 + t[5] # second of day

      if lib == 'tiegcm':
        lat, lon = _tiegcm(yr, doy, *hms)

      if lib == 'apex':
        lat, lon = _apex(yr, doy, *hms)

      if lib == 'laundal':
        lat, lon = _laundal(yr, doy, sod)

      if lib == 'mead':
        lat, lon = _mead(yr, doy, sod)

      rtp[i, 1] = lat
      rtp[i, 2] = lon

    xyz = hxform.sph2car(rtp)

    if lib == 'mead':
      if frame != 'GEI':
        xyz = hxform.transform(xyz, t, 'GEI', frame, lib=lib_transform)
    else:
      if frame != 'GEO':
        xyz = hxform.transform(xyz, t, 'GEO', frame, lib=lib_transform)

    if outer_type in (list, tuple, str):
      return outer_type(xyz[0])

    return xyz


if __name__ == '__main__':

  import hxform
  t = [2018, 7, 1, 12, 0, 0]

  # year
  yr = t[0]

  # day of year
  doy = hxform.timelib.ymd2doy(t[0:3])

  # hours, minutes, seconds
  hms = t[3:6]

  # seconds after midnight (second of day)
  sod = t[3]*3600 + t[4]*60 + t[5]

  subsolar_pt = _mead(yr, doy, sod, lib_transform='geopack_08_dp')
  print(subsolar_pt)
  # (23.07090356829239, 0.9811996644603767)

  # https://github.com/NCAR/apex_fortran/blob/master/test.out
  # subsol inputs: iyr,iday,ihr,imn,sec= 2018  182   12    0   0. UT
  # subsol returns: sbsllat,sbsllon =   23.09    0.97 deg

  subsolar_pt = _tiegcm(yr, doy, *hms)
  print(subsolar_pt)
  # (23.07096778320955, 0.9804334139973321)

  subsolar_pt = _apex(yr, doy, *hms)
  print(subsolar_pt)
  # (23.07096778320955, 0.9804334139973321)

  subsolar_pt = _laundal(2018, doy, sod)
  print(subsolar_pt)
  # (23.08837872675872, 0.9684353956029668)

  subsolar_pt = _hxform(t, frame='GEO')
  rtp = hxform.car2sph(subsolar_pt)
  print(tuple(rtp[1:3]))
  # (23.088330367250137, 0.969196427178192)

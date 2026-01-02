def _nint(a):
  """Replicates Fortran's NINT behavior (rounds halfway cases away from zero)."""
  # https://stackoverflow.com/questions/45884588/how-does-fortran-convert-a-real-number-to-integer
  if a >= 0:
    return int(a + 0.5)
  else:
    return int(a - 0.5)


def _laundal(year, doy, ut):

  from math import floor, sin, cos, pi, atan2, asin

  # Appendix C of https://link.springer.com/article/10.1007/s11214-016-0275-y#appendices
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


def _hxform(time, frame='MAG', lib='geopack_08_dp'):
  """Compute the subsolar point at a given time in a specified frame.
  Usage:
  ------
  subsol_pt = hxform.subsolar_point(time, frame='MAG', lib='geopack_08_dp')
  """

  import numpy
  import hxform

  if lib not in hxform.libs():
    raise ValueError(f"Library '{lib}' not in available libraries: {hxform.libs()}")

  subsol_pt = hxform.transform([1., 0., 0.], time, 'GSM', frame, lib=lib)

  if isinstance(time, (list, tuple, str)):
    return subsol_pt
  else:
    return numpy.array(subsol_pt)


def subsolar_point(t, frame='MAG', lib='geopack_08_dp', lib_transform='geopack_08_dp'):

  import numpy as np

  import hxform

  libs_alt = ['tiegcm', 'apex', 'laundal']
  libs = [*libs_alt, *hxform.libs()]
  if lib not in libs:
    raise ValueError(f"Library '{lib}' not in available libraries: {libs}")

  if lib in libs_alt:
    outer_type = type(t)

    from hxform.transform import _t_prep
    ts, _ = _t_prep(t)
    rtp = np.full((len(ts), 3), np.nan)
    rtp[:, 0] = 1.0  # radius

    for i, t in enumerate(ts):
      # year
      yr = t[0]
      # day of year
      doy = hxform.timelib.doy(t[0:3]) 
      # hours, minutes, seconds
      hms = t[3:6]
      # seconds after midnight
      sf = t[3]*3600 + t[4]*60 + t[5]

      if lib == 'tiegcm':
        lat, lon = _tiegcm(yr, doy, *hms)
      if lib == 'apex':
        lat, lon = _apex(yr, doy, *hms)
      if lib == 'laundal':
        lat, lon = _laundal(yr, doy, sf)

      rtp[i, 1] = lat
      rtp[i, 2] = lon

    xyz = hxform.sph2car(rtp)
    if frame != 'GEO':
      xyz = hxform.transform(xyz, t, 'GEO', frame, lib=lib_transform)

    if outer_type in (list, tuple, str):
      return outer_type(xyz[0])

    return xyz

  return _hxform(t, frame=frame, lib=lib)


if __name__ == '__main__':

  import hxform
  t = [2018, 7, 1, 12, 0, 0]

  # year
  yr = t[0]

  # day of year
  doy = hxform.timelib.doy(t[0:3]) 

  # hours, minutes, seconds
  hms = t[3:6]

  # seconds after midnight
  sf = t[3]*3600 + t[4]*60 + t[5]

  # https://github.com/NCAR/apex_fortran/blob/master/test.out
  # subsol inputs: iyr,iday,ihr,imn,sec= 2018  182   12    0   0. UT
  # subsol returns: sbsllat,sbsllon =   23.09    0.97 deg

  subsolar_pt = _tiegcm(yr, doy, *hms)
  print(subsolar_pt)
  # (23.07096778320955, 0.9804334139973321)

  subsolar_pt = _apex(yr, doy, *hms)
  print(subsolar_pt)
  # (23.07096778320955, 0.9804334139973321)

  subsolar_pt = _laundal(2018, doy, sf)
  print(subsolar_pt)
  # (23.08837872675872, 0.9684353956029668)

  subsolar_pt = _hxform(t, frame='GEO')
  rtp = hxform.car2sph(subsolar_pt)
  print(rtp[1:3])
  # [23.088330367250137, 0.969196427178192]

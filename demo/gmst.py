def is_leap_year(iyear):
  """Return True if a year is a leap year."""

  if not isinstance(iyear, int):
    raise ValueError(f"year must be an integer, got {type(iyear)}")

  if (iyear % 4 == 0 and iyear % 100 != 0) or (iyear % 400 == 0):
    return True
  else:
    return False


def julian_centuries_since_j2000(iyear, iday, secs, method='basic'):
  """Julian centuries since J2000.0 (1 Jan 2000 12:00 UT)."""

  checkargs(iyear, iday, secs)

  return julian_days_since_j2000(iyear, iday, secs, method=method)/36525.0


def julian_days_since_j2000(iyear, iday, secs, method='basic'):
  """Julian days since J2000.0 (2000-01-01T12:00 UT) given year, doy, and secs in Gregorian and Proleptic calendar."""

  checkargs(iyear, iday, secs)

  if method == 'basic':
    # A straightforward implementation of the definition of Julian days since J2000.0.
    days = 0
    if iyear >= 2000:
      iyears = range(2000, iyear)
      sign = 1
    else:
      iyears = range(iyear, 2000)
      sign = -1

    for y in iyears:
      if is_leap_year(y):
        days += sign*366
      else:
        days += sign*365

    return days + (iday - 1) + secs/86400.0 - 0.5

  if method == 'standard':
    import math
    year_jd = iyear - 1
    month_jd = 13
    day_jd = 1
    century = year_jd//100
    gregorian = 2 - century + century//4
    jd = math.floor(365.25*(year_jd + 4716))
    jd += math.floor(30.6001*(month_jd + 1))
    jd += day_jd + gregorian - 1524.5
    jd += (iday - 1) + secs/86400.0
    mjd = jd - 2400000.5

    return mjd - 51544.5


def julian_days(iyear, iday, secs, calendar='julian'):
  """Julian days given year, doy, and secs in Julian or Gregorian calendar."""

  checkargs(iyear, iday, secs)

  if calendar == 'julian':
    return julian_days_since_j2000(iyear, iday, secs, method='basic') + 2451545.0
  elif calendar == 'gregorian':
    return julian_days_since_j2000(iyear, iday, secs, method='standard') + 2451545.0
  else:
    raise ValueError(f"calendar must be 'julian' or 'gregorian', got {calendar}")


def julian_days_since_j2000_test():
  max_diff = 0.0
  for iyear in range(-4800, -4000):
    for iday in range(1, 367):
      if iday == 366 and not is_leap_year(iyear):
        continue
      for secs in [0.0, 43200.0]:
        days_basic = julian_days_since_j2000(iyear, iday, secs, method='basic')
        days_jd = julian_days_since_j2000(iyear, iday, secs, method='standard')
        max_diff = max(max_diff, abs(days_basic - days_jd))
        if days_basic != days_jd:
          print(f"Difference for {iyear}-{iday:03d} {secs:.3f} sec: basic={days_basic:.6f}, standard={days_jd:.6f}")

  if max_diff == 0.0:
    print("No differences found between basic and standard methods.")
  else:
    print(f"Maximum difference between basic and standard methods: {max_diff:.6e} days = {max_diff*86400:.6e} seconds")


def mjd(iyear, iday, secs):
  checkargs(iyear, iday, secs)

  d0 = julian_days_since_j2000(iyear, iday, secs)

  mjd = d0 + 51544.5

  return mjd


def checkargs(iyear, iday, secs):
  if not isinstance(iyear, int):
    raise ValueError(f"year must be an integer, got {type(iyear)}")

  if not isinstance(iday, int):
    raise ValueError(f"iday must be an integer, got {type(iday)}")

  if not isinstance(secs, (int, float)):
    raise ValueError(f"secs must be a number, got {type(secs)}")

  if not secs >= 0.0 and secs < 86400.0:
    raise ValueError(f"secs must be in the range [0, 86400), got {secs}")

  if iday < 1 or iday > 366:
    raise ValueError(f"iday must be in the range [1, 366], got {iday}")

  if not is_leap_year(iyear) and iday == 366:
    raise ValueError(f"iday must be in the range [1, 365] for non-leap years, got {iday} for year {iyear}")


def deg2timestr(deg):
  """Convert degrees in the range [0, 360) to HH:MM:SS.ssssss."""

  if not isinstance(deg, (int, float)):
    raise ValueError(f"deg must be a number, got {type(deg)}")
  if deg < 0.0 or deg >= 360.0:
    raise ValueError(f"deg must be in the range [0, 360), got {deg}")

  hours = int(deg/15.0)
  minutes = int((deg/15.0 - hours)*60.0)
  fseconds = ((deg/15.0 - hours)*60.0 - minutes)*60.0
  seconds = int(fseconds)
  microseconds = int((fseconds - seconds)*1e6)
  ts = f"{hours:02d}:{minutes:02d}:{int(seconds):02d}.{microseconds:06d}"

  return ts


def timestr2deg(ts):
  """Convert time string in the form HH:MM:SS.ssssss to degrees"""

  parts = ts.split(':')
  hours = int(parts[0])
  minutes = int(parts[1])
  seconds = float(parts[2])

  if hours < 0 or hours >= 24:
    raise ValueError(f"hours must be in the range [0, 24), got {hours}")
  if minutes < 0 or minutes >= 60:
    raise ValueError(f"minutes must be in the range [0, 60), got {minutes}")
  if seconds < 0.0 or seconds >= 60.0:
    raise ValueError(f"seconds must be in the range [0, 60), got {seconds}")

  deg = 15*(hours + minutes/60.0 + seconds/3600.0)

  return deg


def gmstR71(iyear, iday, secs):
  """Greenwich Mean Sidereal Time in degrees using Russell 1971 formula."""

  # http://jsoc.stanford.edu/~jsoc/keywords/Chris_Russel/Geophysical%20Coordinate%20Transformations.htm
  """
  The following is from
  http://jsoc.stanford.edu/~jsoc/keywords/Chris_Russel/Geophysical%20Coordinate%20Transformations.htm
  which is a copy of a page previous posted at on Russell's website.

  Appendix 2. The Calculation of the Position of the Sun

  G.D. Mead (private communication) has written a simple subroutine to calculate
  the position of the Sun in GEI coordinates. It is accurate for years 1901
  through 2099, to within 0.006 deg. The input is the year, day of year and
  seconds of the day in UT. The output is Greenwich Mean Sideral Time in degrees,
  the ecliptic longitude, apparent right ascension and declination of the Sun
  in degrees. The listing of this program follows. We note that the cartesian
  coordinates of the vector from the Earth to the Sun are:

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
    DATA RAD /57.29578/
    DOUBLE PRECISION DJ, FDAY
    IF(IYR. LT. 1901. OR. IYR. GT. 2099) RETURN
    FDAY = SECS/86400
    DJ = 365* (IYR-1900) + (IYR-1901)/4 + IDAY + FDAY -0.5D0 
    T = DJ / 36525
    VL = DMOD (279.696678 + 0.9856473354*DJ, 360.D0)
    GST = DMOD (279.690983 + 0.9856473354*DJ + 360.*FDAY + 180., 360.D0)
    G = DMOD (358.475845 + 0.985600267*DJ, 360.D0) / RAD
    SLONG = VL + (1.91946 -0.004789*T)*SIN(G) + 0.020094*SIN (2.*G)
    OBLIQ = (23.45229 -0.0130125*T) / RAD
    SLP = (SLONG -0.005686) / RAD
    SIND = SIN (OBLIQ)*SIN (SLP)
    COSD = SQRT(1.-SIND**2)
    SDEC = RAD * ATAN (SIND/COSD)
    SRASN = 180. -RAD*ATAN2
    (COTAN (OBLIQ)*SIND/COSD, -COS (SLP)/COSD)
    RETURN
    END
  """

  checkargs(iyear, iday, secs)
  if iyear < 1901 or iyear > 3000:
    raise ValueError(f"year must be in the range [1901, 2099], got {iyear}")

  fday = secs/86400.0
  dj = 365*(iyear - 1900) + (iyear - 1901)//4 + iday + fday - 0.5

  gst_deg = (279.690983 + 0.9856473354*dj + 360.0*fday + 180.0) % 360.0

  # Geopack-2008 version of this is
  # DATA RAD/57.295779513D0/
  # GST=DMOD(279.690983D0 + .9856473354D0*DJ + 360.D0*FDAY + 180.D0, 360.D0)/
  # * RAD
  # (the * before RAD is a line continuation character in Fortran, not a multiplication operator)

  """
  Newcomb (https://apps.dtic.mil/sti/tr/pdf/ADA548492.pdf Eqn 11):
  Long of Sun: 279° 41' 48''.04 + 129602768''.13 T + 1.089''T^2
  where T is the time reckoned in Julian centuries of 36525 days since
  January 0, 1900, 12 h UT

  279° 41' 48.04 ≈ 279.6966777777777778 degrees
  129602768.13/36525/3600 ≈ 0.9856473353867213 degrees/day

  The difference in the leading constant in Russell 1971 and Newcomb's formula is
  279.6966777777777778 - 279.690983 = 0.005694777777762283 degrees

  This seems to correspond to the aberration of Δλ = -20.4898'' = 0.00569161111 degrees
  due to light-time correction for the Sun's apparent position (Meeus 1998, pg 133 "... affected
  by aberration (-20''.5) and measured ...").
  """

  return gst_deg


def newcomb(iyear, iday, secs, order=2):
  """Longitude of the Sun in degrees using Newcomb's formula."""
  checkargs(iyear, iday, secs)

  # Equation 11 of McCarthy 2011
  #Long of Sun: 279° 41' 48''.04 + 129602768''.13 T + 1.089''T^2
  #where T is the time reckoned in Julian centuries of 36525 days since
  #January 0, 1900, 12 h UT

  # 36525.0 = julian_days_since_j2000(2000, 1, 86400/2) - julian_days_since_j2000(1899, 365, 86400/2)# 
  dj = julian_days_since_j2000(iyear, iday, secs) + 36525.0

  T = dj/36525.0

  a = 279 + 41/60.0 + 48.04/3600.0
  b = (129602768.13)/3600.0
  c = 1.08936525/3600.0
  if order == 1:
    c = 0.0

  long = (a + b*T + c*T**2) % 360.0

  # J1900.0 = 1899 December 31, 12:00 TT = JD 2415020.0
  # J2000.0 = 2000 January 1, 12:00 TT = JD 2451545.0
  # Difference is 1 Julian Century

  # To convert to J2000.0, replace T with T + 1, which gives
  # L = (a + b*(T+1) + c*(T+1)^2) % 360.0
  # L = (a + b + c^2 + (b + 2*c)*T + c*T^2) % 360.0

  # (a + b + c^2) % 360.0 = 280.46560286934255
  # (b + 2*c) = 36000.76953020292
  # long = (280.46560286934255 + 36000.76953020292*T + c*T**2) % 360.0

  # Meeus 1998, Eqn 12.4 is
  # T = (JD - 2451 545.0)/36525
  #theta0 = 280.460 618 37 + 360.985 647 366 29 (JD - 2451 545.0) + 0.000 387 933*T^2 - T^3/38 710 000
  return long


def gmstF02(iyear, iday, secs):
  """Greenwich Mean Sidereal Time in degrees using Fränz and Harper 2002 formula.

  Fränz and Harper 2002, Section 2
  ------------------------------------------------------------------------------
  The Julian Day Number (JD) starts at Greenwich mean noon 4713 Jan. 1, B.C.
  [S.2.26]. The epoch day number is defined in this paper as the fractional
  number of days of 86400 s from the epoch:

  d0 = (JD - 2451545.0)                                                   (1)

  Formulae from S. and A. use Julian centuries (T0) from J2000.0. One Julian
  century has 36525 days, one Julian year has 365.25 days, s.t.[S.T3.222.2]

  T0 = d0/36525.0   and   y0 = d0/365.25                                  (2)
  ------------------------------------------------------------------------------

  Fränz and Harper 2002, Section 3.3.1
  ------------------------------------------------------------------------------
  For the precision needed in this paper, we may neglect the difference between
  TU and T0, such that (Meeus, 2000):

  GMST ≈ 280◦.46061837
       + 360◦.98564736629*d0
       +   0◦.0003875*T0^2
       -   0◦.000000026*T0^3
  ------------------------------------------------------------------------------

  Note that 1/38710000 (= 2.5833118057349522e-08) is used in Meeus 1998 Eqn 12.4
  (see notes in gmstM98) and 2.6e-08 is used in Fränz and Harper 2002 below (and
  Siedelman 1992). Also, 0.000387933 is used in Meeus 1998 Eqn 12.4
  (and Siedelman 1992) but 0.0003875 is used in Fränz and Harper 2002.

  Finally, Fränz and Harper 2002 cite "Meeus, J, 2000, Astronomical Algorithms
  2nd Edition" but but I am only able to find 1998 associated with a 2nd edition.
  """

  checkargs(iyear, iday, secs)

  T0 = julian_centuries_since_j2000(iyear, iday, secs)
  d0 = julian_days_since_j2000(iyear, iday, secs)

  # Fränz and Harper 2002, Eqn (20)
  theta = 280.46061837 + 360.98564736629*d0 + 0.0003875*T0**2 - (2.6e-8)*T0**3
  """
  Meeus 1998, Eqn 12.4
  theta = 280.46061837 + 360.98564736629*d0 + 0.000387933*T0**2 - (1/38710000)*T0**3
  """
  #print(0.0003875*T0**2)
  #print(0.000387933*T0**2)

  theta = theta % 360.0

  return theta


def gmstM98(iyear, iday, secs):
  """
  Meeus 1998, Chapter 12
  ------------------------------------------------------------------------------
  We shall denote by Theta0 the sidereal time at Greenwich at 0h UT of a given
  date, and by theta0 the sidereal time at Greenwich for any given instant UT.

  The sidereal time at the meridian of Greenwich, at 0h Universal Time of a given
  date, can be obtained as follows.

  Calculate the JD corresponding to that date at 0h UT (Chapter 7). Thus, this
  is a number ending on .5. Then find T by

  (12.1)  T = (JD - 2451 545.0)/36525

  The mean sidereal time at Greenwich at 0h UT is then given by the following
  expression which was adopted in 1982 by the International Astronomical Union:

  (12.2) Theta0 = 6h41m50s.54841 + 8640 184s.812 866*T + 0s.093 104*T^2 - 0s.000 0062*T^3

  Expressed in degrees and decimals, this formula can be written

  (12.3) Theta0 = 100.460 618 37 + 36000.770 053 608*T + 0.000 387 933*T^2 - T^3/38 710 000

  Important: the formulae (12.2) and (12.3) are valid only for those values of T
  which correspond to 0h UT of a date. All other values would give incorrect
  results. To obtain the sidereal time theta0 at Greenwich for any instant UT of a
  given date, multiply that instant by 1.00273790935 and add the result to
  the sidereal time Theta0 at Oh UT.

  To obtain the sidereal time theta0 at Greenwich for any instant UT of a given date,
  multiply that instant by 1.002 737909 35 and add the result to the sidereal time
  Theta0 at 0h UT.

  The mean sidereal time at Greenwich, expressed in degrees, can also be found
  directly for any instant as follows. If JD is the Julian Day corresponding to
  that instant in UT (not necessarily 0h), find T by formula (12.1), and then

  (12.4) theta0 = 280.460 618 37 + 360.985 647 366 29 (JD - 2451 545.0) + 0.000 387 933*T^2 - T^3/38 710 000
  ------------------------------------------------------------------------------
  """

  checkargs(iyear, iday, secs)

  T0 = julian_centuries_since_j2000(iyear, iday, secs)
  d0 = julian_days_since_j2000(iyear, iday, secs)

  theta = 280.46061837 + 360.98564736629*d0 + 0.000387933*T0**2 - (1/38710000)*T0**3
  theta = theta % 360.0

  return theta


def gmstH92(iyear, iday, secs):
  """Compute Greenwich Mean Sidereal Time in degrees using Hapgood 1992 formula.

  Hapgood 1992, Section 4.1
  ------------------------------------------------------------------------------
  All of the fundamental transformations defined in the following sections are
  time dependent. To maintain a uniform style, time is there specified by modified
  Julian date (MJD), which is the time measured in days from 00:00 UT on 17
  November 1858 (Julian date 2400000.5). In this paper we use only the integer
  part of MJD, i.e. the value at 00:00 UT on the day of interest. For some
  applications it is also necessary to give the time within the day as
  Universal Time in hours (UT).
  ------------------------------------------------------------------------------

  Hapgood 1992, Section 4.1
  ------------------------------------------------------------------------------
  The rotation angle theta is the Greenwich mean sidereal time. This can be
  calculated using the following formula (U.S. Naval Observatory, 1989):

  theta = 100.461 + 360000.770*T0 + 15.04107 UT,                          (3)

  where

  T0 = (MJD - 51544.5)/36525.0

  Note that T0 is the time in Julian centuries (36525 days) from 12:00 UT on
  1 January 2000 (known as epoch 2000.0) to the previous midnight.
  ------------------------------------------------------------------------------

  Hapgood 1995, Section 3
  ------------------------------------------------------------------------------
  The Russell (1971) paper is more complex. GST and $\\lambda_{\\odot}$ are given
  by the Fortran subroutine in his Appendix 2. The statements in that routine
  are functionally equivalent to the equations in Hapgood (1992) and thus yield
  values in the mean epoch-of-date.
  ...

  However, the values of GST and $\\lambda_{\\odot}$ predicted by Russell agree
  with those of Hapgood (1992) to within $0.01^\\circ$.
  ------------------------------------------------------------------------------

  Note that equation (3) for theta in Hapgood (1992) is related to Theta0 in 
  Meeus 1998 (see notes in gmstM198). Meeus 1998 Eqn 12.3 is

  Theta0 = 100.46061837 + 36000.770053608*T + 0.000387933*T^2 - T^3/38710000.

  which is followed by "...To obtain the sidereal time theta0 at Greenwich for
  any instant UT of a given date, multiply that instant by 1.00273790935 and add
  the result to the sidereal time Theta0 at 0h UT."

  Equation (3) above has the first two numbers rounded to three decimal places
  and the T^2 and T^3 terms are omitted. The 15.04107 comes from rounding
  15*1.00273790935 = 15.0410686402499998.
  """

  checkargs(iyear, iday, secs)

  T0 = (mjd(iyear, iday, 0.0) - 51544.5)/36525.0

  theta = 100.461 + 36000.770*T0 + 15.04107*(secs/3600.0)
  theta = theta % 360.0

  return theta


if __name__ == "__main__":
  #julian_days_since_j2000_test()
  #exit()
  from utilrsw.xprint import xprint

  """
  Entering data into https://aa.usno.navy.mil/data/siderealtime gives

  Date (UT1)	Time (UT1)  GST (Mean)  GST (Apparent)  LST (Mean)    LST (Apparent) EoE

  https://aa.usno.navy.mil/calculated/siderealtime?date=2000-01-01&time=00:00:00.000&intv_mag=1.00&intv_unit=1&reps=1&lat=0.0000&lon=0.0000&label=&submit=Get+Data
  2000-Jan-1	00:00:00.0	06:39:52.2717	06:39:51.4197	06:39:52.2717 06:39:51.4197  -0.8520

  https://aa.usno.navy.mil/calculated/siderealtime?date=2000-01-01&time=12%3A00%3A00.000&intv_mag=1.00&intv_unit=1&reps=1&lat=0.0000&lon=0.0000&label=&submit=Get+Data
  2000-Jan-1	12:00:00.0	18:41:50.5494	18:41:49.6974	18:41:50.5494	18:41:49.6974	-0.8520

  https://aa.usno.navy.mil/calculated/siderealtime?date=2020-01-01&time=00:00:00.000&intv_mag=1.00&intv_unit=1&reps=1&lat=0.0000&lon=0.0000&label=&submit=Get+Data
  2020-Jan-1	00:00:00.0	06:40:29.2343	06:40:28.2256	06:40:29.2343	06:40:28.2256	-1.0087

  where GST = Greenwich Sidereal Time, LST = Local Sidereal Time, and 
  EoE = Equation of the Equinoxes = GST (Apparent) - GST (Mean) = LST (Apparent) - LST (Mean)
  """

  tests = [
    (2000, 1, 0, "06:39:52.2717"),
    (2000, 1, 43200, "18:41:50.5494"),
    (2020, 1, 0, "06:40:29.2343")
  ]

  for year, iday, secs, usno_str in tests:
    print(80*"-")
    xprint(f"{year}-{iday:03d} {secs:.3f} sec")

    usno = timestr2deg(usno_str)
    R71 = gmstR71(year, iday, secs)
    H92 = gmstH92(year, iday, secs)
    F02 = gmstF02(year, iday, secs)
    M98 = gmstM98(year, iday, secs)

    xprint(f"Newcomb long of Sun (1st order only)  {newcomb(year, iday, secs, order=1):9.6f} deg")
    xprint(f"Newcomb long of Sun (full equation)   {newcomb(year, iday, secs, order=2):9.6f} deg")
    xprint("")

    #ref = "USNO"
    #ref_val = timestr2deg(usno_str)

    ref = "R71"
    ref_val = R71

    xprint(f"USNO         GMST = {usno_str}   | {usno:.6f} deg | {usno - ref_val:9.6f} deg diff from {ref} value")
    xprint(f"Russell 1971 GMST = {deg2timestr(R71)} | {R71:.6f} deg | {R71 - ref_val:9.6f} deg diff from {ref} value")
    xprint(f"Hapgood 1992 GMST = {deg2timestr(H92)} | {H92:.6f} deg | {H92 - ref_val:9.6f} deg diff from {ref} value")
    xprint(f"F & H 2002   GMST = {deg2timestr(F02)} | {F02:.6f} deg | {F02 - ref_val:9.6f} deg diff from {ref} value")
    xprint(f"Meeus 1998   GMST = {deg2timestr(M98)} | {M98:.6f} deg | {M98 - ref_val:9.6f} deg diff from {ref} value")

    xprint(f"\nF & H 2002 - Meeus 1998 {F02 - M98:.16e}")

    print(80*"-")
    print("")


  if False:
    year = 2200
    iday = 1
    secs = 0.0

    print("Days since J2000.0 = ", julian_days_since_j2000(year, iday, secs))
    print("Days since J2000.0 = ", 36525.0*julian_centuries_since_j2000(year, iday, secs))
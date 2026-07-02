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

  if iday < 1 or iday > 366:
    raise ValueError(f"iday must be in the range [1, 366], got {iday}")
  if not is_leap_year(iyear) and iday == 366:
    raise ValueError(f"iday must be in the range [1, 365] for non-leap years, got {iday} for year {iyear}")
  if not secs >= 0.0 and secs < 86400.0:
    raise ValueError(f"secs must be in the range [0, 86400), got {secs}")


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

  checkargs(iyear, iday, secs)
  if iyear < 1901 or iyear > 3000:
    raise ValueError(f"year must be in the range [1901, 2099], got {iyear}")

  fday = secs/86400.0
  dj = 365*(iyear - 1900) + (iyear - 1901)//4 + iday + fday - 0.5

  gst_deg = (279.690983 + 0.9856473354*dj + 360.0*fday + 180.0) % 360.0

  """
  On iyear = 1900 iday = 1, fday = 0.5 => dj = 0, so
  gst_deg = (279.690983 + 360.0*0.5 + 180.0) % 360.0
  gst_deg = (279.690983 + 360.0) % 360.0

  On iyear = 1901 iday = 1, fday = 0.5 => dj = 365, so
  0.9856473354*dj = 359.761277421 and

  gst_deg = (279.690983 + 359.761277421 + 360) % 360
  gst_deg = (639.452260421 + 360) % 360

  Newcomb 1898 pg 9 gives

  𝜏 = 18^h 38^m 45^s.836 + 8 640 184^s.542T + 0.0929T^2

  where T is # the number of Julian centuries of 36525 days since
  1900, Jan 0, Greenwich Mean Noon.

  This can be written as

  𝜏 = 279.69098333333335° + (8640184.542)*15/3600 T
    = 279.69098333333335° + 36000.768925 T

  Note that 36000.768925/36525 = 0.9856473353867213, which matches the second
  term in Russell's gst_deg equation to the precision given.

  When T = 0
  𝜏 = 279.69098333333335°, which matches the R71 case (even though equation
  is said to not be valid for iyear = 1900).

  On 1901, Jan 1, Greenwich Mean Noon, T = 365/36525,
  𝜏 = 279.69098333333335° + 36000.768925 * (365/36525)
  𝜏 = 639.452 260 749 487°

  Compare this with the R71 value of
     (639.452 260 421 + 360) % 360
  the difference is due to rounding.
  """

  return gst_deg



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

  Fran̈z and Harper 2002, Section 1, last par
  ------------------------------------------------------------------------------
  We also cite the formulae and methods given by Hapgood (1992) for geocentric
  systems, which are based on the Astronomical Almanac for Computers (1988)
  which is no longer updated by the Nautical Almanac Offices. The formulae used
  by Hapgood (1992) are first order approximations of the third order formulae
  given in Expl.Suppl. (1961).
  ------------------------------------------------------------------------------

  Note that 1/38710000 = 2.5833118057349522e-08 is used in Meeus 1998 Eqn 12.4
  (see notes in gmstM98) and 2.6e-08 is used in Fränz and Harper 2002 (and
  Seidelmann 1992). Also, 0.000387933*T^2 used in Meeus 1998 Eqn 12.4
  (and Seidelmann 1992) but 0.0003875*T0^2 is used in Fränz and Harper 2002.
  Finally, Fränz and Harper 2002 cite "Meeus, J, 2000, Astronomical Algorithms
  2nd Edition" but but I am only  able to find 1998 associated with a 2nd edition.
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
  the sidereal time Theta0 at 0h UT.

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

  theta = 100.461 + 36000.770*T0 + 15.04107 UT,

  where

  T0 = (MJD - 51544.5)/36525.0                                             (3)

  Note that T0 is the time in Julian centuries (36525 days) from 12:00 UT on
  1 January 2000 (known as epoch 2000.0) to the previous midnight.
  ------------------------------------------------------------------------------

  Note 1:
    U.S. Naval Observatory, 1989 does not give this formula explicitly.
    It gives (page B2, eqn. 1)

    GMST = 6h.69737456 + 2400.051336 T_0 + 0.0000258622 T_0^2 + 1.002737909 UT

    with definitions on page B3:

    To and T are time intervals in Julian centuries from J2000.0:
      To = (JDo - 2451545.0) / 36525    T = (JD - 2451545.0) / 36525;
    UT is the universal time in hours;
    JDo and JD are the Julian dates at 0h UT and at an arbitrary time of the day,
    respectively;

    U.S. Naval Observatory, 1989 and Hapgood 1992 are equivalent because
      360*6.69737456/24 = 100.4606184000000

      1.002737909*15 = 15.041068634999998

      MJD - 51544.5 = JD - 2451545.0

      15*2400.051336 = 36000.77004000000

    and Hapgood 1992 is

    theta = 100.461 + 36000.770*T0 + 15.04107 UT

  Note 2:
    Equation (3) for theta is related to Theta0 in Meeus 1998 (see notes in
    gmstM198). Meeus 1998 Eqn 12.3 is

    Theta0 = 100.46061837 + 36000.770053608*T + 0.000387933*T^2 - T^3/38710000.

    which is followed by "... To obtain the sidereal time theta0 at Greenwich for 
    any instant UT of a given date, multiply that instant by 1.00273790935 and add
    the result to the sidereal time Theta0 at 0h UT."

    Equation (3) above has the first two numbers rounded to the nearest 0.001
    degree and the T^2 and T^3 terms are omitted. The 15.04107 comes from rounding
    15*1.00273790935 = 15.0410686402499998.
  """

  checkargs(iyear, iday, secs)

  T0 = (mjd(iyear, iday, 0.0) - 51544.5)/36525.0

  theta = 100.461 + 36000.770*T0 + 15.04107*(secs/3600.0)
  theta = theta % 360.0

  return theta


def IAU2000():
  pass
  # https://www.celestialprogramming.com/snippets/greenwichMeanSiderealTime.html


if __name__ == "__main__":
  """
  Entering data into https://aa.usno.navy.mil/data/siderealtime gives

  https://aa.usno.navy.mil/calculated/siderealtime?date=1899-12-31&time=12%3A00%3A00.000&intv_mag=1.00&intv_unit=1&reps=1&lat=0.0000&lon=0.0000&label=&submit=Get+Data
  Date (UT1)	Time (UT1)  GST (Mean)  GST (Apparent)  LST (Mean)    LST (Apparent) EoE
  1899-Dec-31	12:00:00.0	18:38:45.8477	18:38:46.9084	18:38:45.8477	18:38:46.9084	1.0607

  https://aa.usno.navy.mil/calculated/siderealtime?date=1901-12-31&time=12%3A00%3A00.000&intv_mag=1.00&intv_unit=1&reps=1&lat=0.0000&lon=0.0000&label=&submit=Get+Data
  Date (UT1)	Time (UT1)  GST (Mean)  GST (Apparent)  LST (Mean)    LST (Apparent) EoE
  1901-Dec-31	12:00:00.0	18:36:51.2622	18:36:51.9875	18:36:51.2622	18:36:51.9875	0.7253

  https://aa.usno.navy.mil/calculated/siderealtime?date=1970-01-01&time=12:00:00.000&intv_mag=1.00&intv_unit=1&reps=1&lat=0.0000&lon=0.0000&label=&submit=Get+Data
  Date (UT1)	Time (UT1)  GST (Mean)  GST (Apparent)  LST (Mean)    LST (Apparent) EoE
  1970-Jan-1	12:00:00.0	18:42:53.3971	18:42:53.6719	18:42:53.3971	18:42:53.6719	0.2748

  https://aa.usno.navy.mil/calculated/siderealtime?date=1980-01-01&time=12:00:00.000&intv_mag=1.00&intv_unit=1&reps=1&lat=0.0000&lon=0.0000&label=&submit=Get+Data
  Date (UT1)	Time (UT1)  GST (Mean)  GST (Apparent)  LST (Mean)    LST (Apparent) EoE
  1980-Jan-1	12:00:00.0	18:41:13.5942	18:41:13.1176	18:41:13.5942	18:41:13.1176	-0.4766

  https://aa.usno.navy.mil/calculated/siderealtime?date=1990-01-01&time=12:00:00.000&intv_mag=1.00&intv_unit=1&reps=1&lat=0.0000&lon=0.0000&label=&submit=Get+Data
  Date (UT1)	Time (UT1)  GST (Mean)  GST (Apparent)  LST (Mean)    LST (Apparent) EoE
  1990-Jan-1	12:00:00.0	18:43:30.3485	18:43:31.0725	18:43:30.3485	18:43:31.0725	0.7239

  https://aa.usno.navy.mil/calculated/siderealtime?date=2000-01-01&time=12:00:00.000&intv_mag=1.00&intv_unit=1&reps=1&lat=0.0000&lon=0.0000&label=&submit=Get+Data
  Date (UT1)	Time (UT1)  GST (Mean)  GST (Apparent)  LST (Mean)    LST (Apparent) EoE
  2000-Jan-1	12:00:00.0	18:41:50.5494	18:41:49.6974	18:41:50.5494	18:41:49.6974	-0.8520

  https://aa.usno.navy.mil/calculated/siderealtime?date=2010-01-01&time=12:00:00.000&intv_mag=1.00&intv_unit=1&reps=1&lat=0.0000&lon=0.0000&label=&submit=Get+Data
  2010-Jan-1	12:00:00.0	18:44:07.3074	18:44:08.3189	18:44:07.3074	18:44:08.3189	1.0115

  https://aa.usno.navy.mil/calculated/siderealtime?date=2020-01-01&time=12:00:00.000&intv_mag=1.00&intv_unit=1&reps=1&lat=0.0000&lon=0.0000&label=&submit=Get+Data
  Date (UT1)	Time (UT1)  GST (Mean)  GST (Apparent)  LST (Mean)    LST (Apparent) EoE
  2020-Jan-1	12:00:00.0	18:42:27.5120	18:42:26.5019	18:42:27.5120	18:42:26.5019	-1.0101

  where GST = Greenwich Sidereal Time, LST = Local Sidereal Time, and 
  EoE = Equation of the Equinoxes = GST (Apparent) - GST (Mean) = LST (Apparent) - LST (Mean)
  """
  usno_data = {
    #(1899, 365, 43200.0): "18:38:45.8477", # Too early for R71
    (1901, 365, 43200.0): "18:36:51.2622",
    (1970, 1, 43200.0): "18:42:53.3971",
    (1980, 1, 43200.0): "18:41:13.5942",
    (1990, 1, 43200.0): "18:43:30.3485",
    (2000, 1, 43200.0): "18:41:50.5494",
    (2010, 1, 43200.0): "18:44:07.3074",
    (2020, 1, 43200.0): "18:42:27.5120",
  }

  for key in usno_data:

    gmstUSNO_str = gmstUSNO_str = usno_data[key]
    year, iday, secs = key

    _gmstUSNO = timestr2deg(gmstUSNO_str)
    _gmstR71 = gmstR71(year, iday, secs)
    _gmstH92 = gmstH92(year, iday, secs)
    _gmstF02 = gmstF02(year, iday, secs)
    _gmstM98 = gmstM98(year, iday, secs)

    #ref = "USNO"
    #ref_val = timestr2deg(gmstUSNO_str)
    ref = "R71"
    ref_val = _gmstR71

    print(40*"-")
    print(f"{year}-{iday:03d} {secs/3600.0:.2f} h")
    print(40*"-")
    print(f"USNO         GMST = {gmstUSNO_str}   | {_gmstUSNO:.6f} deg | {_gmstUSNO - ref_val:9.6f} deg diff from {ref} value")
    print(f"Russell 1971 GMST = {deg2timestr(_gmstR71)} | {_gmstR71:.6f} deg | {_gmstR71 - ref_val:9.6f} deg diff from {ref} value")
    print(f"Hapgood 1992 GMST = {deg2timestr(_gmstH92)} | {_gmstH92:.6f} deg | {_gmstH92 - ref_val:9.6f} deg diff from {ref} value")
    print(f"F & H 2002   GMST = {deg2timestr(_gmstF02)} | {_gmstF02:.6f} deg | {_gmstF02 - ref_val:9.6f} deg diff from {ref} value")
    print(f"Meeus 1998   GMST = {deg2timestr(_gmstM98)} | {_gmstM98:.6f} deg | {_gmstM98 - ref_val:9.6f} deg diff from {ref} value")
    print(40*"-")

  if False:
    year = 2200
    iday = 1
    secs = 0.0

    print("Days since J2000.0 = ", julian_days_since_j2000(year, iday, secs))
    print("Days since J2000.0 = ", 36525.0*julian_centuries_since_j2000(year, iday, secs))
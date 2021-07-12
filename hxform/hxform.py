import numpy as np
from datetime import datetime
import time as tm # angel addon for checking timing.py

def tpad(time, length=7):
    """Pad list with 3 or more elements with zeros.

    Example:
    --------
    >>> from hxform import hxform as hx
    >>> print(hx.tpad([2000,1,1]))                 # [2000, 1, 1, 0, 0, 0, 0]
    >>> print(hx.tpad([2000,1,1], length=4))       # [2000, 1, 1, 0]
    >>> print(hx.tpad([2000,1,1,2,3,4], length=3)) # [2000, 1, 1]
    """

    in_type = type(time)

    # TODO: Check that time is valid
    time = np.array(time)

    assert(len(time) > 2)

    if len(time.shape) == 1:
        if len(time) > length:
            time = time[0:length]
        else:
            pad = length - len(time)
            time = np.pad(time, (0,pad), 'constant', constant_values=0)
    else:
        if len(time[0]) > length:
            time = time[:,0:length]
        else:
            pad = length - len(time)
            time = np.pad(time, ((0,0),(0,pad)), 'constant', constant_values=0)



    if in_type == np.ndarray:
        return time
    elif in_type == tuple:
        return tuple(map(tuple,time))
    else:
        return list(map(list,time))


def  is_leap_year(year):
    if isinstance(year, int):
        if year % 100 == 0:
            return year % 400 == 0
        return year % 4 == 0
    else:
        year = np.array(year)
        leap400 = np.where(year%400==0, True,False)
        leap100 = np.where(year%100==0, True, False)
        leap4 = np.where(year%4==0, True, False)
        cor1 = np.where(leap100, False, leap4)
        return np.where(leap400, True, cor1)

def doy(date):
    """
    doy([2021,7,5])  ->  186
    doy([[2000,9,30],[1900,9,30],[1980,5,5]]) -> [274, 273. 126]
    """
    date = np.array(date)

    if len(date.shape) == 1:
        year, month, day = date[0], date[1], date[2]

        if is_leap_year(year): K = 1
        else: K = 2
    else:
        year, month, day = date[:,0], date[:,1], date[:,2]

        K = np.where(is_leap_year(year),1, 2)

    N = np.fix((275* month) / 9.0) - K * np.fix((month + 9) / 12.0) + day - 30

    return N.astype(int)

def to_doy(t):
    """Convert from [y, m, d, h, min, sec] to [y, doy, h, min, sec].

    Example
    -------
    >>> to_doy([2000,2,1,9,9,9]) # [2000,32,9,9,9]
    """
    in_type = type(t)

    t = np.array(t)

    if len(t.shape) == 1:
        pad = 6 - len(t)
        t = np.pad(t, (0,pad), 'constant', constant_values=0)
    else:
        pad = 6 - len(t[0])
        t = np.pad(t, ((0,0),(0,pad)), 'constant', constant_values=0)

    if len(t.shape) == 1:
        day_of_year = datetime(*t).timetuple().tm_yday
        doy_list = np.array([t[0], day_of_year, t[3], t[4], t[5]])
    else:
        day_of_year = doy(t[:,:3])
        t = np.column_stack((t[:,0], day_of_year, t[:,3], t[:,4], t[:,5]))


    if in_type == np.ndarray:
        return t
    elif in_type == tuple:
        return tuple(t.tolist()) #tuple(map(tuple,t))
    else:
        return t.tolist()



def transform(v, time, csys_in, csys_out, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Transfrom between coordinates systems using Geopack or SpacePy.

    Parameters
    ----------
    v : array-like

        (Nv, 3) float np.array

        np.array of three floats

        list of three floats

        list containing lists of three floats

        list of 3-element np.arrays

    time : array-like
           list of 3+ ints
           list containing lists of 3+ ints
           np.array of 3+ ints
           (Nt, 3) float np.array, where Nt = 1 or Nt = Nv

           The 3+ ints are [year, month, day, [hours, [minutes, [seconds]]]]
           Zeros are used for any missing optional value.

    csys_in : str
              One of MAG, GEI, GEO, GSE, GSM, SM

    csys_out : str
               One of MAG, GEI, GEO, GSE, GSM, SM

    ctype_in : str
               'car' (default) or 'sph'
               For spherical coordinates, `v` should be in r, latitude, longitude,
               with angles in degrees.

    ctype_out : str
               'car' (default) or 'sph'

    lib : str
          'geopack_08_dp' (default) or 'spacepy'

    Returns
    -------
    array-like with dimensions matching either `time` (if `Nt` != 1 and `Nv` = 1) or
    `v` (if `Nv` =! 1 and `Nt` = 1). If `Nv` and `Nt` != 1, dimensions are same as `v`.

    Return type will match that of `v`. Note that if a list of 3-element np.arrays are
    passed, execution time will be larger. Use `np.ndarrays` for `v` and `time` for fastest
    execution time.

    Examples
    --------
    >>> from hxform import hxform as hx
    >>> t1 = [2000, 1, 1, 0, 0, 0] # or np.array([2000, 1, 1, 0, 0, 0])
    >>> v1 = [0, 0, 1]             # or np.array([0, 0, 1])
    >>> # All of the following are equivalent and return a list with three floats
    >>> from hxform import hxform as hx
    >>> hx.transform(v1, time1, 'GSM', 'GSE')
    >>> hx.transform(v1, time1, 'GSM', 'GSE', ctype_in='car')
    >>> hx.transform(v1, time1, 'GSM', 'GSE', ctype_in='car', ctype_out='car')

    The following 3 calls return a list with two lists of 3 elements

    1. Transform two vectors at same time t1

        >>> from hxform import hxform as hx
        >>> hx.transform([v1, v1], t1, 'GSM', 'GSE')

    2. Transform one vector at two different times

        >>> from hxform import hxform as hx
        >>> hx.transform(v1, [t1, t1], 'GSM', 'GSE')

    3. Transform two vectors, each at different times

        >>> from hxform import hxform as hx
        >>> hx.transform([v1, v1], [t1, t1], 'GSM', 'GSE')
    """

    assert(lib in ['cxform', 'geopack_08_dp', 'spacepy'])
    in_type = type(v)

    list_of_arrays = False
    if isinstance(v[0], np.ndarray) and isinstance(v, list):
        list_of_arrays = True

    v = np.array(v, dtype=np.double)
    time = np.array(time, dtype=np.int32)

    if len(time.shape) > 1 and len(v.shape) > 1:
        if time.shape[0] != v.shape[0]:
            raise ValueError("time and v cannot be different lengths")

    if lib == 'cxform':
        import os
        import glob
        import ctypes
        this_script = os.path.join(os.path.dirname(__file__),
                                "../hxform/cxform_wrapper*")

        for lib_file in glob.glob(this_script):
            # The name of the .so or .dll file will not be the same on all
            # systems, so we need to find it. (For example, on one system
            # it is cxform_wrapper.cpython-37m-darwin.so.)
            # TODO: Find a better way to do this.
            break
        lib_path = os.path.join(lib_file)
        lib_obj = ctypes.cdll.LoadLibrary(lib_path)

        assert(len(time.shape) <= 2 and len(v.shape) <= 2)

        if len(v.shape) == 1:
            v = np.array([v])

        if len(time.shape) == 1:
            time = np.array([time])

        Nt = time.shape[0]

        if time.shape[1] < 3:
            raise ValueError("At least year, month, and day must be given for time.")

        if ctype_in == 'sph':
            v[:,0], v[:,1], v[:,2] = StoC(v[:,0], v[:,1], v[:,2])

        vp = np.full(v.shape, np.nan)

        nz = time.shape[1]
        if nz != 6:
            # Pad time. TODO: Do this in wrapper so extra memory is not needed.
            tmp = np.zeros((time.shape[0], 6-nz), dtype=np.int32)
            time = np.concatenate((time, tmp), 1)

        ret = lib_obj.cxform_wrapper(
                ctypes.c_void_p(v.ctypes.data),
                ctypes.c_void_p(time.ctypes.data),
                ctypes.c_char_p(str.encode(csys_in)),
                ctypes.c_char_p(str.encode(csys_out)),
                ctypes.c_void_p(vp.ctypes.data),
                ctypes.c_int(v.shape[0]),
                ctypes.c_int(int(Nt))
            )

        if ctype_out == 'sph':
            vp[:,0], vp[:,1], vp[:,2] = CtoS(vp[:,0], vp[:,1], vp[:,2])

    if lib == 'geopack_08_dp':
        import hxform.geopack_08_dp as geopack_08_dp
        trans = csys_in + 'to' + csys_out

        if len(v.shape) == 1:
            v = np.array([v])

        if len(time.shape) == 1:
            time = np.array([time])

        dtime = np.array(to_doy(time))

        if v.shape[0] <= time.shape[0]:
            outsize = time.shape[0]
        else:
            outsize = v.shape[0]

        if ctype_in == 'sph':
            v[:,0], v[:,1], v[:,2] = StoC(v[:,0], v[:,1], v[:,2])

        vp = geopack_08_dp.transform(v, trans, dtime, outsize)

        if ctype_out == 'sph':
            vp[:,0], vp[:,1], vp[:,2] = CtoS(vp[:,0], vp[:,1], vp[:,2])


    if lib == 'spacepy':
        try:
            # SpacePy is not installed when hxform is installed due to
            # frequent install failures and so the default is to not use it.
            import spacepy.coordinates as sc
            from spacepy.time import Ticktock
            import numpy.matlib
        except ImportError as error:
            print(error.__class__.__name__ + ": " + error.message)
        except Exception as exception:
            print(exception, False)
            print(exception.__class__.__name__ + ": " + exception.message)

        if len(time.shape) == 1 and len(v.shape) > 1:
            time = numpy.matlib.repmat(time, v.shape[0], 1)
        if len(v.shape) == 1 and len(time.shape) > 1:
            v = numpy.matlib.repmat(v, time.shape[0], 1)
        if len(v.shape) == 1:
            v = np.array([v])

        cvals = sc.Coords(v, csys_in, ctype_in)

        if len(time.shape) == 1:
            # SpacePy requires time values to be strings with second precision
            t_str = '%04d-%02d-%02dT%02d:%02d:%02d' % tuple(tpad(time, length=6))
        else:
            t_str = []
            for i in range(time.shape[0]):
                t_str.append('%04d-%02d-%02dT%02d:%02d:%02d' % tuple(tpad(time[i,:], length=6)))
            t_str = np.array(t_str)

        cvals.ticks = Ticktock(t_str, 'ISO')
        newcoord = cvals.convert(csys_out, ctype_out)

        vp = newcoord.data

    if len(time.shape) == 1 and len(v.shape) == 1:
        vp = vp[0, :]

    if in_type == np.ndarray:
        return vp
    else:
        if list_of_arrays is True:
            vp2 = []
            for i in range(vp.shape[0]):
                vp2.append(vp[i])
            return vp2
        else:
            return vp.tolist()


def MAGtoGEI(v, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Equivalent to transform(v, time, 'MAG', 'GEI', ...)"""
    return transform(v, time, 'MAG', 'GEI', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def MAGtoGEO(v, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Equivalent to transform(v, time, 'MAG', 'GEO', ...)"""
    return transform(v, time, 'MAG', 'GEO', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def MAGtoGSE(v, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Equivalent to transform(v, time, 'MAG', 'GSE', ...)"""
    return transform(v, time, 'MAG', 'GSE', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def MAGtoGSM(v, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Equivalent to transform(v, time, 'MAG', 'GSM', ...)"""
    return transform(v, time, 'MAG', 'GSM', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def MAGtoSM(v, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Equivalent to transform(v, time, 'MAG', 'SM', ...)"""
    return transform(v, time, 'MAG', 'SM', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)


def GEItoMAG(v, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Equivalent to transform(v, time, 'GEI', 'MAG', ...)"""
    return transform(v, time, 'GEI', 'MAG', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def GEItoGEO(v, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Equivalent to transform(v, time, 'GEI', 'GEO', ...)"""
    return transform(v, time, 'GEI', 'GEO', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def GEItoGSE(v, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Equivalent to transform(v, time, 'GEI', 'GSE', ...)"""
    return transform(v, time, 'GEI', 'GSE', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def GEItoGSM(v, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Equivalent to transform(v, time, 'GEI', 'GSM', ...)"""
    return transform(v, time, 'GEI', 'GSM', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def GEItoSM(v, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Equivalent to transform(v, time, 'GEI', 'SM', ...)"""
    return transform(v, time, 'GEI', 'SM', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)


def GEOtoMAG(v, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Equivalent to transform(v, time, 'GEO', 'MAG', ...)"""
    return transform(v, time, 'GEO', 'MAG', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def GEOtoGEI(v, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Equivalent to transform(v, time, 'GEO', 'GEI', ...)"""
    return transform(v, time, 'GEO', 'GEI', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def GEOtoGSE(v, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Equivalent to transform(v, time, 'GEO', 'GSE', ...)"""
    return transform(v, time, 'GEO', 'GSE', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def GEOtoGSM(v, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Equivalent to transform(v, time, 'GEO', 'GSM', ...)"""
    return transform(v, time, 'GEO', 'GSM', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def GEOtoSM(v, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Equivalent to transform(v, time, 'GEO', 'SM', ...)"""
    return transform(v, time, 'GEO', 'SM', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)


def GSEtoMAG(v, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Equivalent to transform(v, time, 'GSE', 'MAG', ...)"""
    return transform(v, time, 'GSE','MAG', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def GSEtoGEI(v, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Equivalent to transform(v, time, 'GSE', 'GEI', ...)"""
    return transform(v, time, 'GSE','GEI', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def GSEtoGEO(v, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Equivalent to transform(v, time, 'GSE', 'GEO', ...)"""
    return transform(v, time, 'GSE','GEO', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def GSEtoGSM(v, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Equivalent to transform(v, time, 'GSE', 'GEO', ...)"""
    return transform(v, time, 'GSE','GSM', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def GSEtoSM(v, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Equivalent to transform(v, time, 'GSE', 'SM', ...)"""
    return transform(v, time, 'GSE','SM', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)


def GSMtoMAG(v, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Equivalent to transform(v, time, 'GSM', 'MAG', ...)"""
    return transform(v, time, 'GSM', 'MAG', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def GSMtoGEI(v, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Equivalent to transform(v, time, 'GSM', 'GEI', ...)"""
    return transform(v, time, 'GSM', 'GEI', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def GSMtoGEO(v, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Equivalent to transform(v, time, 'GSM', 'GEO', ...)"""
    return transform(v, time, 'GSM', 'GEO', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def GSMtoGSE(v, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Equivalent to transform(v, time, 'GSM', 'GSE', ...)"""
    return transform(v, time, 'GSM', 'GSE', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def GSMtoSM(v, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Equivalent to transform(v, time, 'GSM', 'SM', ...)"""
    return transform(v, time, 'GSM', 'SM', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)


def SMtoMAG(v, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Equivalent to transform(v, time, 'SM', 'MAG', ...)"""
    return transform(v, time, 'SM', 'MAG', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def SMtoGEI(v, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Equivalent to transform(v, time, 'SM', 'GEI', ...)"""
    return transform(v, time, 'SM', 'GEI', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def SMtoGEO(v, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Equivalent to transform(v, time, 'SM', 'GEO', ...)"""
    return transform(v, time, 'SM', 'GEO', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def SMtoGSE(v, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Equivalent to transform(v, time, 'SM', 'GSE', ...)"""
    return transform(v, time, 'SM', 'GSE', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def SMtoGSM(v, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Equivalent to transform(v, time, 'SM', 'GSM', ...)"""
    return transform(v, time, 'SM', 'GSM', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)


def CtoS(x, y, z):
    """Convert from cartesian to r, colatitude [degrees], longitude [degrees]."""
    r = np.sqrt(np.power(x, 2) + np.power(y, 2) + np.power(z, 2))
    assert np.all(r > 0), 'radius must be greater than zero.'
    theta = 90.0 - (180.0/np.pi)*np.arccos(z/r)
    phi = (180.0/np.pi)*np.arctan2(y, x)

    return r, theta, phi

def StoC(r, colat, long):
    """Convert r, colatitude [degrees], longitude [degrees] to cartesian."""
    assert np.all(r > 0), 'radius must be greater than zero.'
    x = r*np.cos((np.pi/180.0)*long)*np.cos((np.pi/180.0)*colat)
    y = r*np.sin((np.pi/180.0)*long)*np.cos((np.pi/180.0)*colat)
    z = r*np.cos(np.pi/2.0 - (np.pi/180.0)*colat)

    return x, y, z


def UTtoHMS(UT, **kwargs):
    """Convert universal time in fractional hours into integer hour, minutes, seconds.

    Example
    -------
    >>> from hxform import hxform as hx
    >>> print(hx.UTtoHMS(12))              # [12, 0, 0]
    >>> print(hx.UTtoHMS(24))              # [0, 0, 0]
    >>> print(hx.UTtoHMS(24, keep24=True)) # [24, 0, 0]
    """

    keep24 = False
    if 'keep24' in kwargs:
        keep24 = kwargs['keep24']

    if UT > 24 or UT < 0:
        raise ValueError('Required: 0 <= UT <= 24.')

    hours = int(UT)
    minutes = int((UT-hours)*60.)
    seconds = int(round((UT-hours-minutes/60.)*3600.))
    if seconds == 60:
        seconds = 0
        minutes = minutes + 1
    if minutes == 60:
        minutes = 0
        hours = hours + 1

    if hours == 24 and keep24 == False:
        return [0, 0, 0]

    return [hours, minutes, seconds]


def MAGtoMLT(pos, time, csys='sph', lib='geopack_08_dp'):
    """Compute magnetic local time given a UT and MAG position or longitude.

    Uses equation 93 in Laundal and Richmond, 2016 (10.1007/s11214-016-0275-y)

    Usage:
    ------
    >>> from hxform import hxform as hx
    >>> mlt = hx.MAGtoMLT(MAGlong, time)
    >>> mlt = hx.MAGtoMLT([MAGlong1, Mlong2, ...], time)

    >>> mlt = hx.MAGtoMLT([MAGx, MAGy, MAGz], time, csys='car')
    >>> mlt = hx.MAGtoMLT([[MAGx1, MAGy1, MAGz1],...], time, csys='car')

    Returns:
    --------
    mlt: float or array-like

    Examples:
    --------
    >>> from hxform import hxform as hx
    >>> mlt = hx.MAGtoMLT(0., [2000, 1, 1, 0, 0, 0])
    >>> print(mlt) # 18.869936573301775

    >>> from hxform import hxform as hx
    >>> mlt = hx.MAGtoMLT([0., 0.], [2000, 1, 1, 0, 0, 0])
    >>> print(mlt) # [18.86993657 18.86993657]

    >>> from hxform import hxform as hx
    >>> mlt = hx.MAGtoMLT([-1., 0., 0.], [2000, 1, 1, 0, 0, 0], csys='car')
    >>> print(mlt) # 6.869936573301775

    >>> from hxform import hxform as hx
    >>> mlt = hx.MAGtoMLT([[-1., 0., 0.],[-1., 0., 0.]], [2000, 1, 1, 0, 0, 0], csys='car')
    >>> print(mlt) # [6.86993657 6.86993657]
"""

    assert(csys == 'car' or csys == 'sph')

    pos = np.array(pos)
    time = np.array(time)

    if not isinstance(pos, float):
        pos = np.array(pos)

    if csys == 'sph':
        phi = pos*np.pi/180.
    else:
        if pos.shape == (3, ):
            phi = np.arctan2(pos[1], pos[0])
        else:
            phi = np.arctan2(pos[:, 1], pos[:, 0])

    subsol_pt = transform(np.array([1., 0., 0.]), time, 'GSM', 'MAG', lib=lib)

    if len(subsol_pt.shape) == 1:
        phi_cds = np.arctan2(subsol_pt[1], subsol_pt[0])
    else:
        phi_cds = np.arctan2(subsol_pt[:, 1], subsol_pt[:, 0])

    delta = phi - phi_cds # note np.array([a1, a2, ...])+b == np.array([a1+b, a2+b, ...])

    if isinstance(delta, float):
        delta = np.array([delta])

    idx = np.where(delta > np.pi)
    delta[idx] = delta[idx] - 2.*np.pi
    idx = np.where(delta <= -np.pi)
    delta[idx] = delta[idx] + 2.*np.pi

    if delta.size == 1:
        delta = delta[0]

    MLT = 12. + delta*24./(2.*np.pi)
    return MLT

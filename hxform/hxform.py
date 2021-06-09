import os
import sys
import numpy as np

def tpad(time, length=7):

    # TODO: Check that time is valid
    time = list(time)

    assert(len(time) > 2)

    if len(time) > length:
        time = time[0:length]
    else:
        pad = length - len(time)
        time = time + pad*[0]

    # TODO: If time was tuple, return tuple.
    #       If time was np.array, return np.array.
    return tuple(time)

def transform(v, time, csys_in, csys_out, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    """Transfrom between coordinates systems using Geopack or SpacePy.

    Parameters
    ----------
    v : array-like
        list of three floats
        list containing lists of three floats
        np.array of three floats
        (Nv, 3) float np.array

    time : array-like
        list of 3+ ints
        list containing lists of 3+ ints
        np.array of 3+ ints
        (Nt, 3) float np.array, where Nt = 1 or Nt = Nv

        The 3+ ints are [year, month, day, [hours, [minutes, [seconds]]]]
        Zeros are used for any missing optional value.

    csys_in : str

    csys_out : str

    ctype_in : str
               'car' (default) or 'sph'

    ctype_out : str
               'car' (default) or 'sph'

    lib : str
          'geopack_08_dp' (default) or 'spacepy'

    Returns
    -------
    array-like with structure matching either `time` (if `Nt` != 1) or
    `v` (if `Nv` =! 1). If `Nv` and `Nt` != 1, size is same as `v`.

    Examples
    --------
    t1 = [2000, 1, 1, 0, 0, 0]
    t2 = [2000, 1, 1, 2, 0, 0]
    v1 = [0, 0, 1]
    v2 = [0, 1, 0]

    # All equivalent and return a list with three floats
    transform([0, 1, 1], time1)
    transform([0, 1, 1], time1, ctype_in='car')
    transform([0, 1, 1], time1, ctype_in='car', ctype_out='car')

    # The following 3 calls return a list with two lists of 3 elements
    # 1. Transform two vectors at same time t1
    transform([v1, v2], t1)

    # 2. Transform two vectors, each at different times
    transform([v1, v2], [t1, t2])

    # 3. Transform one vector at two different times
    transform(v1, [t1, t2])

    # 4. Transform one vector at three times
    transform(np.array([v1, v2]), t1)

    """

    in_type = type(v)
    v = np.array(v)
    t = np.array(time)

    if len(t.shape) > 1 and len(v.shape) > 1:
        if t.shape[0] != v.shape[0]:
            raise ValueError("t and v cannot be different lengths")

    if lib == 'geopack_08_dp':
        import hxform.geopack_08_dp as geopack_08_dp
        trans = csys_in + 'to' + csys_out
        if len(t.shape) == 1:
            if len(v.shape) == 1:
                v = np.array([v])
            #print(v)
            dtime = np.array(tpad(t, length=5), dtype=np.int32)
            ret = np.column_stack(geopack_08_dp.transform(v[:,0], v[:,1], v[:,2], trans, dtime))
        else:
            if len(v.shape) == 1:
                ret = np.full((t.shape[0], 3), np.nan)
                for i in range(0, t.shape[0]):
                    dtime = np.array(tpad(t[i,0:5], length=5), dtype=np.int32)
                    ret[i,:] = np.column_stack(geopack_08_dp.transform(v[0], v[1], v[2], trans, dtime))
            else:
                ret = np.full((v.shape[0], 3), np.nan)
                for i in range(0, t.shape[0]):
                    dtime = np.array(tpad(t[i,0:5], length=5), dtype=np.int32)
                    ret[i,:] = np.column_stack(geopack_08_dp.transform(v[i,0], v[i,1], v[i,2], trans, dtime))

    if lib == 'spacepy':
        try:
            # SpacePy is not installed when hxform is installed due to
            # frequent install failures and the default is to not use it.
            import spacepy.coordinates as sc
            from spacepy.time import Ticktock
            import numpy.matlib
        except ImportError as error:
            print(error.__class__.__name__ + ": " + error.message)
        except Exception as exception:
            print(exception, False)
            print(exception.__class__.__name__ + ": " + exception.message)

        if len(t.shape) == 1 and len(v.shape) > 1:
            t = numpy.matlib.repmat(t, v.shape[0], 1)
        if len(v.shape) == 1 and len(t.shape) > 1:
            v = numpy.matlib.repmat(v, t.shape[0], 1)
        if len(v.shape) == 1:
            v = np.array([v])

        #print(v)
        cvals = sc.Coords(v, csys_in, ctype_in)

        if len(t.shape) == 1:
            # SpacePy requires time values to be strings with second precision
            t_str = '%04d-%02d-%02dT%02d:%02d:%02d' % tpad(t, length=6)
        else:
            t_str = []
            for i in range(t.shape[0]):
                t_str.append('%04d-%02d-%02dT%02d:%02d:%02d' % tpad(t[i,:], length=6))
            t_str = np.array(t_str)

        cvals.ticks = Ticktock(t_str, 'ISO')
        newcoord = cvals.convert(csys_out, ctype_out)

        ret = newcoord.data

    if len(t.shape) == 1 and len(v.shape) == 1:
        ret = ret[0, :]

    if in_type == np.ndarray:
        return ret
    else:
        return ret.tolist()

def MAGtoGSM(v_MAG, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    return transform(v_MAG, time, 'MAG', 'GSM', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def GSMtoMAG(v_GSM, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    return transform(v_GSM, time, 'GSM', 'MAG', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def GEOtoGSM(v_GEO, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    return transform(v_GEO, time, 'GEO', 'GSM', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def GEOtoMAG(v_GEO, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    return transform(v_GEO, time, 'GEO', 'MAG', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def MAGtoGEO(v_MAG, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    return transform(v_MAG, time, 'MAG', 'GEO', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def MAGtoSM(v_MAG, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    return transform(v_MAG, time, 'MAG', 'SM', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def MAGtoGEI(v_MAG, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    return transform(v_MAG, time, 'MAG', 'GEI', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def GEOtoGEI(v_GEO, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    return transform(v_GEO, time, 'GEO', 'GEI', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def GEOtoSM(v_GEO, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    return transform(v_GEO, time, 'GEO', 'SM', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def GSEtoGSM(v_GSE, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    return transform(v_GSE, time, 'GSE', 'GSM', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def GSMtoGSE(v_GSM, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    return transform(v_GSM, time, 'GSM', 'GSE', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def SMtoGSM(v_SM, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    return transform(v_SM, time, 'SM', 'GSM', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def GSMtoSM(v_SM, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    return transform(v_SM, time, 'GSM', 'SM', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def GEItoGEO(v_GEI, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp'):
    return transform(v_GEI, time, 'GEI','GEO', ctype_in=ctype_in, ctype_out=ctype_out, lib=lib)

def StoC(r, theta, phi):
    """Convert from spherical to cartesian coordinates

    r, theta, phi: array-like

    """
    x = r*np.cos(phi)*np.sin(theta)
    y = r*np.sin(phi)*np.sin(theta)
    z = r*np.cos(theta)
    return x, y, z


def UTtoHMS(UT, **kwargs):
    """Convert universal time in fractional hours into integer hour, minutes, seconds

    from cxtransform import UTtoHMS

    print(UTtoHMS(12))              # [12, 0, 0]

    print(UTtoHMS(24))              # [0, 0, 0]

    print(UTtoHMS(24, keep24=True)) # [24, 0, 0]
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


def MAGtoMLT(pos, time, csys='sph', lib='geopack_08_dp', debug=False):
    """Compute magnetic local time given a UT and MAG position or longitude

    Uses equation 93 in https://arxiv.org/abs/1611.10321

    Usage:
    ------
    import cxtransform as cx
    mlt = cx.MAGtoMLT(MAGlong, time)
    mlt = cx.MAGtoMLT([MAGlong1, Mlong2, ...], time)

    mlt = cx.MAGtoMLT([MAGx, MAGy, MAGz], time, csys='car')
    mlt = cx.MAGtoMLT([[MAGx1, MAGy1, MAGz1],...], time, csys='car')

    Returns:
    --------
    mlt: float or array-like

    Examples:
    --------
    import cxtransform as cx

    mlt = cx.MAGtoMLT(0., [2000, 1, 1, 0, 0, 0])
    print(mlt) # 18.869936573301775

    mlt = cx.MAGtoMLT([0., 0.], [2000, 1, 1, 0, 0, 0])
    print(mlt) # [18.86993657 18.86993657]

    mlt = cx.MAGtoMLT([-1., 0., 0.], [2000, 1, 1, 0, 0, 0], csys='car')
    print(mlt) # 6.869936573301775

    mlt = cx.MAGtoMLT([[-1., 0., 0.],[-1., 0., 0.]], [2000, 1, 1, 0, 0, 0], csys='car')
    print(mlt) # [6.86993657 6.86993657]
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

    if debug:
        print('phi =' + str(phi))

    subsol_pt = transform(np.array([1, 0, 0]), time, 'GSM', 'MAG', lib=lib)

    #import pdb;pdb.set_trace()
    if len(subsol_pt.shape) == 1:
        phi_cds = np.arctan2(subsol_pt[1], subsol_pt[0])
    else:
        phi_cds = np.arctan2(subsol_pt[:, 1], subsol_pt[:, 0])

    if debug:
        print('phi_cds =' + str(phi_cds))

    delta = phi - phi_cds # note np.array([a1, a2, ...])+b == np.array([a1+b, a2+b, ...])

    if debug:
        print('delta =' + str(delta))

    if isinstance(delta,float):
        delta = np.array([delta])

    idx = np.where(delta > np.pi)
    delta[idx] = delta[idx] - 2.*np.pi
    idx = np.where(delta <= -np.pi)
    delta[idx] = delta[idx] + 2.*np.pi

    if delta.size == 1:
        delta = delta[0]

    MLT = 12. + delta*24./(2.*np.pi)
    return MLT

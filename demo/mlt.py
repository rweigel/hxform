import hxform
lib = 'geopack_08_dp'
t = '2000-01-01T00:00:00'

mlong = 0.0
mlt = hxform.mlt(mlong, t, lib=lib)
hxform.xprint(mlt) # 18.8708344

mlt = hxform.mlt([mlong, mlong], t, lib=lib)
hxform.xprint(mlt) # [18.8708344 18.8708344]

mlt = hxform.mlt([-1., 0., 0.], t, lib=lib, csys='car')
hxform.xprint(mlt) # 6.8708344027652934

mlt = hxform.mlt([[-1., 0., 0.], [-1., 0., 0.]], [2000, 1, 1, 0, 0, 0], lib=lib, csys='car')
hxform.xprint(mlt) # [6.8708344 6.8708344]

libs = hxform.libs()
hxform.xprint(f'MLT at 0 deg MAG long at {t}Z')
for lib in libs:
  try:
    mlt = hxform.mlt(0., [2000, 1, 1, 0, 0, 0], lib=lib)
    hxform.xprint(f'{lib:13s}: {mlt}')
  except Exception as e:
    hxform.xprint(f'{lib:13s}: Error: {e}')
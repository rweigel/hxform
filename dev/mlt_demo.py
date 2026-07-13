import hxform

t = '2000-01-01T00:00:00'
mlong = 0.0

libs = hxform.libs()
hxform.xprint(f'MLT at {mlong} deg MAG long at {t}Z')
for lib in libs:
  try:
    mlt = hxform.mlt(0., [2000, 1, 1, 0, 0, 0], lib=lib)
    hxform.xprint(f'{lib:14s}: {mlt}')
  except Exception as e:
    hxform.xprint(f'{lib:14s}: Error: {e}')
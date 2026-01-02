import hxform
t = [2018, 7, 1, 12, 0, 0]

libs_alt = ['tiegcm', 'apex', 'laundal']
libs = [*libs_alt, *hxform.libs()]
for lib in libs:
  try:
    print(f'Library: {lib}')
    xyz = hxform.subsolar_point(t, frame='GEO', lib=lib)
    rtp = hxform.car2sph(xyz)
  except Exception as e:
    print(f'  Error: {e}')
    continue

  hxform.xprint(f'  x, y, z: {xyz[0]:.6f}, {xyz[1]:.6f}, {xyz[2]:.6f}')
  hxform.xprint(f'  Lat, Lon: {rtp[1]:.6f}, {rtp[2]:.6f}')

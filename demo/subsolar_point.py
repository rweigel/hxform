import hxform

def subsolars(t, frame='GEO', libs=None, libs_exclude=None):
  #t = [2018, 7, 1, 12, 0, 0]

  if libs is None:
    libs_alt = ['tiegcm', 'apex', 'laundal', 'mead']
    libs = [*libs_alt, *hxform.libs()]

  # Filter out excluded libraries
  if libs_exclude is not None:
    libs = [lib for lib in libs if lib not in libs_exclude]

  thetas = []
  phis = []
  for lib in libs:

    try:
      print(f'Library: {lib}')
      xyz = hxform.subsolar_point(t, frame=frame, lib=lib)
      rtp = hxform.car2sph(xyz) # radius, theta (colat), phi (lon)
      thetas.append(rtp[1])
      phis.append(rtp[2])
      hxform.xprint(f'  x, y, z: {xyz[0]:.6f}, {xyz[1]:.6f}, {xyz[2]:.6f}')
      hxform.xprint(f'  Lat, Lon: {rtp[1]:.6f}, {rtp[2]:.6f}')
    except Exception as e:
      thetas.append(None)
      phis.append(None)
      print(f'  Error: {e}')

  return libs, thetas, phis

import datetime
import pandas
import utilrsw

to = [2010, 1, 1, 0, 0, 0]
tf = [2010, 1, 2, 0, 0, 0]
delta = {'hours': 1}
t = utilrsw.time.ints_list(to, tf, delta)

times = []
thetas = []
phis = []
for ti in t:
  times.append(datetime.datetime(*ti))
  libs, theta, phi = subsolars(ti, frame='GEO', libs_exclude=['sscweb'])
  thetas.append(theta)
  phis.append(phi)

thetas = pandas.DataFrame(thetas, columns=libs, index=times)
print(thetas)
phis = pandas.DataFrame(phis, columns=libs, index=times)
print(phis)

import matplotlib.pyplot as plt

plt.figure(figsize=(10, 6))
for lib in libs:
  plt.plot(times, thetas[lib], marker='o', label=lib)
plt.xlabel('Time')
plt.title('Subsolar colatitude in GEO frame')
plt.legend()
plt.grid(True)
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()
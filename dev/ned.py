def ned(v_cart, x_cart):

  import numpy

  v_sph = get_spherical_vector_components(v_cart, x_cart)
  v_north = -v_sph[:,1]
  v_east  = +v_sph[:,2]
  v_down  = -v_sph[:,0]

  return numpy.column_stack([v_north, v_east, v_down])


def get_spherical_vector_components(v_cart, x_cart):

  import numpy as np
  x = x_cart[:,0]
  y = x_cart[:,1]
  z = x_cart[:,2]
  vx = v_cart[:,0]
  vy = v_cart[:,1]
  vz = v_cart[:,2]

  r = np.sqrt(x**2 + y**2 + z**2)
  L = np.sqrt(x**2 + y**2)

  v_r = (x*vx + y*vy + z*vz)/r
  v_theta = ((x*vx + y*vy)*z - (L**2)*vz)/(r*L)
  v_phi = (-y*vx + x*vy)/L

  return np.column_stack([v_r, v_theta, v_phi])



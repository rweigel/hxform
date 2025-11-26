# Test car2sph and sph2car

import numpy as np

import hxform


def test_basic():
    car0 = np.array([[1., 1., 1.], [2., 2., 2.]])
    sph = hxform.car2sph(car0[:,0], car0[:,1], car0[:,2])
    car = hxform.sph2car(sph[0], sph[1], sph[2])

    for i in range(3):
      np.testing.assert_allclose(car0[:,i], car[i], atol=1e-15, rtol=0.0)

def test_components2matrix():
  out1 = np.array([[1.0, 0.0, 0.0]])
  out2 = np.array([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])

  # Three inputs
  # Case 1. in docstring
  sph = hxform.components2matrix(1.0, 0.0, 0.0)
  assert(np.all(sph == out1))

  # Case 2. in docstring
  sph = hxform.components2matrix([1.0], [0.0], [0.0])
  assert(np.all(sph == out1))

  sph = hxform.components2matrix([1.0, 2.0], [0.0, 0.0], [0.0, 0.0])
  assert(np.all(sph[0, :] == out2[0, :]))
  assert(np.all(sph[1, :] == out2[1, :]))

  sph = hxform.components2matrix((1.0, 0.0, 0.0))
  assert(np.all(sph == out1))

  # Case 3. in docstring
  sph = hxform.components2matrix((1.0), (0.0), (0.0))
  assert(np.all(sph == out1))

  sph = hxform.components2matrix((1.0, 2.0), (0.0, 0.0), (0.0, 0.0))
  assert(np.all(sph[0, :] == out2[0, :]))
  assert(np.all(sph[1, :] == out2[1, :]))

  # One input
  # Case 4. in docstring
  sph = hxform.components2matrix([1.0, 0.0, 0.0])
  assert(np.all(sph == out1))

  # Case 5. in docstring
  sph = hxform.components2matrix((1.0, 0.0, 0.0))
  assert(np.all(sph == out1))

  # Case 6. in docstring
  sph = hxform.components2matrix([1.0], [0.0], [0.0])
  assert(np.all(sph == out1))

  sph = hxform.components2matrix((1.0), (0.0), (0.0))
  assert(np.all(sph == out1))

  # Case 7. in docstring
  sph = hxform.components2matrix([1.0, 2.0], [0.0, 0.0], [0.0, 0.0])
  assert(np.all(sph[0, :] == out2[0, :]))
  assert(np.all(sph[1, :] == out2[1, :]))

  sph = hxform.components2matrix((1.0, 2.0), (0.0, 0.0), (0.0, 0.0))
  assert(np.all(sph[0, :] == out2[0, :]))
  assert(np.all(sph[1, :] == out2[1, :]))

  # Case 8. in docstring
  sph = hxform.components2matrix(np.array([1.0, 0.0, 0.0]))
  assert(np.all(sph == out1))

  # Case 9. in docstring
  out1c = out1.copy()
  out1c.shape = (1, 3)
  sph = hxform.components2matrix(out1c)
  assert(np.all(sph == out1))

  # Case 10. in docstring
  sph = hxform.components2matrix(out2)
  assert(np.all(sph == out2))


def test_simple_values():

    # x = 1., y = 0., z = 0.
    sph = hxform.car2sph(1.0, 0.0, 0.0)
    np.testing.assert_allclose(sph[0], 1.0, atol=1e-15, rtol=0.0)
    np.testing.assert_allclose(sph[1], 0.0, atol=1e-15, rtol=0.0)
    np.testing.assert_allclose(sph[2], 0.0, atol=1e-15, rtol=0.0)

    sph = hxform.sph2car(1.0, 0.0, 0.0)
    np.testing.assert_allclose(sph[0], 1.0, atol=1e-15, rtol=0.0)
    np.testing.assert_allclose(sph[1], 0.0, atol=1e-15, rtol=0.0)
    np.testing.assert_allclose(sph[2], 0.0, atol=1e-15, rtol=0.0)


    # x = 0., y = 1., z = 0.
    sph = hxform.car2sph(0.0, 1.0, 0.0)
    np.testing.assert_allclose(sph[0], 1.0, atol=1e-15, rtol=0.0)
    np.testing.assert_allclose(sph[1], 0.0, atol=1e-15, rtol=0.0)
    np.testing.assert_allclose(sph[2], 90.0, atol=1e-15, rtol=0.0)

    car = hxform.sph2car(1.0, 0.0, 90.0)
    np.testing.assert_allclose(0., car[0], atol=1e-15, rtol=0.0)
    np.testing.assert_allclose(1., car[1], atol=1e-15, rtol=0.0)
    np.testing.assert_allclose(0., car[2], atol=1e-15, rtol=0.0)


    # x = 0., y = 0., z = 1.
    sph = hxform.car2sph(0.0, 0.0, 1.0)
    np.testing.assert_allclose(sph[0], 1.0, atol=1e-15, rtol=0.0)
    np.testing.assert_allclose(sph[1], 90.0, atol=1e-15, rtol=0.0)
    np.testing.assert_allclose(sph[2], 0.0, atol=1e-15, rtol=0.0)

    car = hxform.sph2car(1.0, 90.0, 0.0)
    np.testing.assert_allclose(0.0, car[0], atol=1e-15, rtol=0.0)
    np.testing.assert_allclose(0.0, car[1], atol=1e-15, rtol=0.0)
    np.testing.assert_allclose(1.0, car[2], atol=1e-15, rtol=0.0)


def test_r_equal_zero():
  import pytest
  car0 = np.array([0., 0., 0.])
  with pytest.raises(AssertionError):
    hxform.car2sph(*car0)


def test_random_values():

    rng = np.random.default_rng(1)
    car0 = -1.0 + 2.0*rng.random(size=(1000, 3))
    sph = hxform.car2sph(car0[:,0], car0[:,1], car0[:,2])
    car = hxform.sph2car(sph[0], sph[1], sph[2])

    for i in range(3):
      np.testing.assert_allclose(car0[:,i], car[i], atol=1e-12, rtol=0.0)


if __name__ == '__main__':
  # Could use: https://stackoverflow.com/a/28643737 to call all functions.
  test_components2matrix()
  exit()
  test_basic()
  test_simple_values()
  test_r_equal_zero()
  test_random_values()

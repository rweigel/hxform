# Test car2sph and sph2car

import numpy as np

import hxform


def test_basic():

  kwargs = {'atol': 1e-15, 'rtol': 0.0}
  sph = hxform.car2sph(0, 0, 1)
  np.testing.assert_allclose(sph[0], 1.0, **kwargs)
  np.testing.assert_allclose(sph[1], 90.0, **kwargs)
  np.testing.assert_allclose(sph[2], 0.0, **kwargs)

  car = hxform.sph2car(1.0, 90.0, 0.0)
  np.testing.assert_allclose(car[0], 0.0, **kwargs)
  np.testing.assert_allclose(car[1], 0.0, **kwargs)
  np.testing.assert_allclose(car[2], 1.0, **kwargs)


def assert_raises(exc_type, func, args, match=None, **kwargs):
  """Call a function and assert that it raises an exception.

  Args:
      exc_type (Exception): _expected exception type
      func (function): _function to call_
      args (list): _call function using func(*args)_
      match (str, optional): _string to match in exception message_. Defaults to None.

  Raises:
      AssertionError: _Expected assertion to raise exc_type_
      AssertionError: _Exception message does not match_
      AssertionError: _Exception not raised_

  Examples:
      # Assert that ValueError is raised when calling myfunc(arg1, arg2)
      assert_raises(ValueError, myfunc, [arg1, arg2])

      # Assert that ValueError is raised with message containing 'invalid value'
      # when calling myfunc(arg1, arg2)
      assert_raises(ValueError, myfunc, [arg1, arg2], match='invalid value')
  """
  try:
    func(*args, **kwargs)
  except Exception as e:
    if not isinstance(e, exc_type):
      emsg = f"Expected assertion to raise {exc_type}, got {type(e)}"
      raise AssertionError(emsg)
    if match and match not in str(e):
      raise AssertionError(f"Exception message '{str(e)}' does not contain '{match}'")
    return e
  raise AssertionError(f"{exc_type} not raised")


def test_matrix2components():

  comps0 = (1.0, 0.0, 0.0)

  # hxform.components2matrix(x, y, z) where x, y, z are numbers
  mat = hxform.components2matrix(*comps0)
  comps = hxform.matrix2components(*comps0, mat)
  assert(comps == comps0)

  # hxform.components2matrix([x, y, z]) where x, y, z are numbers
  mat = hxform.components2matrix(list(comps0))
  comps = hxform.matrix2components(list(comps0), mat)
  assert(comps == comps0)

  # hxform.components2matrix((x, y, z)) where x, y, z are numbers
  mat = hxform.components2matrix(comps0)
  comps = hxform.matrix2components(comps0, mat)
  assert(comps == comps0)


  comps0l = ([1.0], [0.0], [0.0])
  comps0t = ((1.0, ), (0.0, ), (0.0, ))

  # hxform.components2matrix(x, y, z) where x, y, z are lists
  mat = hxform.components2matrix(*comps0l)
  comps = hxform.matrix2components(*comps0l, mat)
  assert(comps == comps0)

  # hxform.components2matrix([x, y, z]) where x, y, z are lists
  mat = hxform.components2matrix(comps0l)
  comps = hxform.matrix2components(comps0l, mat)
  assert(comps == comps0)

  # hxform.components2matrix([x, y, z]) where x, y, z are tuples
  mat = hxform.components2matrix(list(comps0t))
  comps = hxform.matrix2components(list(comps0t), mat)
  assert(comps == comps0)

  # hxform.components2matrix((x, y, z)) where x, y, z are tuples
  mat = hxform.components2matrix(comps0t)
  comps = hxform.matrix2components(comps0t, mat)
  assert(comps == comps0)

  # hxform.components2matrix((x, y, z)) where x, y, z are lists
  mat = hxform.components2matrix(tuple(comps0t))
  comps = hxform.matrix2components(tuple(comps0t), mat)
  assert(comps == comps0)


  comps0l = ([1.0, 2.0], [0.0, 0.0], [0.0, 0.0])
  comps0t = ((1.0, 2.0), (0.0, 0.0), (0.0, 0.0))

  # hxform.components2matrix(x, y, z) where x, y, z are lists
  mat = hxform.components2matrix(*comps0l)
  comps = hxform.matrix2components(*comps0l, mat)
  assert(comps == comps0l)

  # hxform.components2matrix(x, y, z) where x, y, z are tuples
  mat = hxform.components2matrix(*comps0t)
  comps = hxform.matrix2components(*comps0t, mat)
  assert(comps == comps0t)

  # hxform.components2matrix([x, y, z]) where x, y, z are lists
  mat = hxform.components2matrix(comps0l)
  comps = hxform.matrix2components(comps0l, mat)
  assert(comps == comps0)

  # hxform.components2matrix([x, y, z]) where x, y, z are tuples
  mat = hxform.components2matrix(list(comps0t))
  comps = hxform.matrix2components(list(comps0t), mat)
  assert(comps == comps0)

  # hxform.components2matrix((x, y, z)) where x, y, z are tuples
  mat = hxform.components2matrix(comps0t)
  comps = hxform.matrix2components(comps0t, mat)
  assert(comps == comps0)

  # hxform.components2matrix((x, y, z)) where x, y, z are lists
  mat = hxform.components2matrix(tuple(comps0t))
  comps = hxform.matrix2components(tuple(comps0t), mat)
  assert(comps == comps0)


  comps0l = [np.array(1), np.array(0), np.array(0)]
  comps0t = (np.array(1), np.array(0), np.array(0))
  # hxform.components2matrix(x, y, z) where x, y, z are scalar numpy.ndarrays
  mat = hxform.components2matrix(*comps0l)
  comps = hxform.matrix2components(*comps0l, mat)
  assert(comps == comps0l)

  # hxform.components2matrix([x, y, z]) where x, y, z are scalar numpy.ndarrays
  mat = hxform.components2matrix(comps0l)
  comps = hxform.matrix2components(comps0l, mat)
  assert(comps == comps0)

  # hxform.components2matrix((x, y, z)) where x, y, z are scalar numpy.ndarrays
  mat = hxform.components2matrix(comps0t)
  comps = hxform.matrix2components(comps0t, mat)
  assert(comps == comps0)


  comps0l = (np.array([1.0, 2.0]), np.array([0.0, 0.0]), np.array([0.0, 0.0]))
  comps0t = (np.array([1.0, 2.0]), np.array([0.0, 0.0]), np.array([0.0, 0.0]))

  # hxform.components2matrix(x, y, z) where x, y, z are shape = (2, ) numpy.ndarrays
  mat = hxform.components2matrix(*comps0l)
  comps = hxform.matrix2components(*comps0l, mat)
  for i in range(3):
    assert(np.all(comps[i] == comps0l[i]))

  # hxform.components2matrix([x, y, z]) where x, y, z are shape = (2, ) numpy.ndarrays
  mat = hxform.components2matrix(comps0l)
  comps = hxform.matrix2components(comps0l, mat)
  for i in range(3):
    assert(np.all(comps[i] == comps0l[i]))

  # hxform.components2matrix((x, y, z)) where x, y, z are shape = (2, ) numpy.ndarrays
  mat = hxform.components2matrix(comps0t)
  comps = hxform.matrix2components(comps0t, mat)
  for i in range(3):
    assert(np.all(comps[i] == comps0t[i]))




def test_components2matrix():

  comps1 = np.array([[1.0, 0.0, 0.0]])
  comps2 = np.array([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])

  # Three inputs

  # Case 1. in docstring
  mat = hxform.components2matrix(1.0, 0.0, 0.0)
  assert(mat.shape == (1, 3))
  assert(np.all(mat == comps1))

  match = "Number of arguments must be 1 or 3"
  # hxform.components2matrix(1, 2)
  assert_raises(ValueError, hxform.components2matrix, [1, 2], match=match)
  # hxform.components2matrix(1, 2, 3, 4)
  assert_raises(ValueError, hxform.components2matrix, [1, 2, 3, 4], match=match)
  # hxform.components2matrix([1, 2, 3], 1)
  assert_raises(ValueError, hxform.components2matrix, [[1, 2, 3], 1], match=match)

  match = 'If the first input is a number, all inputs must numbers'
  # hxform.components2matrix(1, [2], 3)
  assert_raises(ValueError, hxform.components2matrix, [1, [2], 3], match=match)
  # hxform.components2matrix(1, (2,), 3)
  assert_raises(ValueError, hxform.components2matrix, [1, (2,), 3], match=match)
  match = 'If the first input is not numeric, all inputs have same type'
  # hxform.components2matrix([1], 2, 3)
  assert_raises(ValueError, hxform.components2matrix, [[1], 2, 3], match=match)

  match = 'All inputs must be number, list, tuple, or numpy.ndarray'
  # hxform.components2matrix({}, {}, {})
  assert_raises(ValueError, hxform.components2matrix, [{}, {}, {}], match=match)

  # Case 2. in docstring
  match = 'All arguments must have same length'
  # hxform.components2matrix([1, 2], [0], [0])
  assert_raises(ValueError, hxform.components2matrix, [[1, 2], [0], [0]], match=match)
  # hxform.components2matrix((1, 2), (0, ), (0, ))
  assert_raises(ValueError, hxform.components2matrix, [(1, 2), (0, ), (0, )], match=match)

  mat = hxform.components2matrix([1.0], [0.0], [0.0])
  assert(mat.shape == (1, 3))
  assert(np.all(mat == comps1))

  mat = hxform.components2matrix([1.0, 2.0], [0.0, 0.0], [0.0, 0.0])
  assert(mat.shape == (2, 3))
  assert(np.all(mat[0, :] == comps2[0, :]))
  assert(np.all(mat[1, :] == comps2[1, :]))

  # Case 3. in docstring
  mat = hxform.components2matrix((1.0), (0.0), (0.0))
  assert(mat.shape == (1, 3))
  assert(np.all(mat == comps1))

  mat = hxform.components2matrix((1.0, 2.0), (0.0, 0.0), (0.0, 0.0))
  assert(mat.shape == (2, 3))
  assert(np.all(mat[0, :] == comps2[0, :]))
  assert(np.all(mat[1, :] == comps2[1, :]))

  # Case 4. in docstring
  mat = hxform.components2matrix(np.array(1), np.array(0), np.array(0))
  assert(mat.shape == (1, 3))
  assert(np.all(mat == comps1))

  mat = hxform.components2matrix(np.array([1, 2]), np.array([0, 0]), np.array([0, 0]))
  assert(mat.shape == (2, 3))
  assert(np.all(mat == comps2))

  match = 'All numpy.ndarrays must be 1D'
  arg = [np.array([[1, 2], [3, 4]]), np.array([0, 0]), np.array([0, 0])]
  # hxform.components2matrix(arg[0])
  assert_raises(ValueError, hxform.components2matrix, arg, match=match)

  match = 'All numpy.ndarrays must have same shape[0]'
  arg = [np.array([1, 2]), np.array([1, 2]), np.array([1, 2, 3])]
  # hxform.components2matrix(arg[0])
  assert_raises(ValueError, hxform.components2matrix, arg, match=match)


  # One input

  match = 'Input list/tuple must have three elements'
  # hxform.components2matrix([1.0, 0.0])
  assert_raises(ValueError, hxform.components2matrix, [[1.0, 0.0]], match=match)

  # Case 5. in docstring
  mat = hxform.components2matrix([1.0, 0.0, 0.0])
  assert(mat.shape == (1, 3))
  assert(np.all(mat == comps1))

  # Case 6. in docstring
  mat = hxform.components2matrix((1.0, 0.0, 0.0))
  assert(mat.shape == (1, 3))
  assert(np.all(mat == comps1))

  # Case 7. in docstring
  mat = hxform.components2matrix([[1, 0, 0], [2, 0, 0]])
  assert(mat.shape == (2, 3))
  assert(np.all(mat == comps2))

  # Mis-match of number types allowed. Coerced to common type by numpy.array()
  mat = hxform.components2matrix([[1.0, 0, 0], [2, 0, 0]])
  assert(mat.shape == (2, 3))
  assert(np.all(mat == comps2))

  # Case 8. in docstring
  mat = hxform.components2matrix(((1, 0, 0), (2, 0, 0)))
  assert(mat.shape == (2, 3))
  assert(np.all(mat == comps2))

  # Case 7 and 8 error cases
  match = 'If the first input is not numeric, all inputs must have same type'
  assert_raises(ValueError, hxform.components2matrix, [[[1, 0, 0], (2, 0, 0)]], match=match)

  # Case 9. in docstring
  mat = hxform.components2matrix(np.array([1.0, 0.0, 0.0]))
  assert(np.all(mat == comps1))

  match = 'Input 1-D numpy.ndarray must have three elements'
  # hxform.components2matrix(np.array([1.0, 0.0]))
  assert_raises(ValueError, hxform.components2matrix, [np.array([1.0, 0.0])], match=match)

  # Case 10. in docstring
  comps1c = comps1.copy()
  comps1c.shape = (1, 3)
  mat = hxform.components2matrix(comps1c)
  assert(mat.shape == (1, 3))
  assert(np.all(mat == comps1))

  match = 'Input 1-D numpy.ndarray must have three elements'
  # hxform.components2matrix(np.array([1.0, 0.0]))
  assert_raises(ValueError, hxform.components2matrix, [np.array([1.0, 0.0])], match=match)

  match = 'Input numpy.ndarray must be 1D or 2D'
  # hxform.components2matrix(np.array([[[1.0, 0.0, 0.0]]]))
  assert_raises(ValueError, hxform.components2matrix, [np.array([[[1.0, 0.0, 0.0]]])], match=match)

  # Case 11. in docstring
  mat = hxform.components2matrix(comps2)
  assert(mat.shape == (2, 3))
  assert(np.all(mat == comps2))

  match = 'Input numpy.ndarray must have three columns'
  arg = [np.array([[1.0, 0.0, 0.0, 0.0]])]
  # hxform.components2matrix(arg)
  assert_raises(ValueError, hxform.components2matrix, arg, match=match)

  arg = [np.array([[1.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0]])]
  # hxform.components2matrix(arg)
  assert_raises(ValueError, hxform.components2matrix, arg, match=match)


def test_simple_values():

    # x = 1., y = 0., z = 0.
    sph = hxform.car2sph(1.0, 0.0, 0.0)
    np.testing.assert_allclose(sph[0], 1.0, atol=1e-15, rtol=0.0)
    np.testing.assert_allclose(sph[1], 0.0, atol=1e-15, rtol=0.0)
    np.testing.assert_allclose(sph[2], 0.0, atol=1e-15, rtol=0.0)

    car = hxform.sph2car(1.0, 0.0, 0.0)
    np.testing.assert_allclose(car[0], 1.0, atol=1e-15, rtol=0.0)
    np.testing.assert_allclose(car[1], 0.0, atol=1e-15, rtol=0.0)
    np.testing.assert_allclose(car[2], 0.0, atol=1e-15, rtol=0.0)

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
  #test_matrix2components()
  #test_basic()
  #test_simple_values()
  #test_r_equal_zero()
  #test_random_values()

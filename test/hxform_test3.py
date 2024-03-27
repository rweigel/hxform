# Test CtoS and StoC

import pytest
import numpy as np
from hxform import hxform as hx

def test_basic():
    car0 = np.array([[1., 1., 1.], [2., 2., 2.]])
    sph = hx.CtoS(car0[:,0], car0[:,1], car0[:,2])
    car = hx.StoC(sph[0], sph[1], sph[2])

    for i in range(3):
        np.testing.assert_allclose(car0[:,i], car[i], atol=0.0, rtol=1e-10)

def test_simple_values():

    # x = 1., y = 0., z = 0.
    sph = hx.CtoS(1.0, 0.0, 0.0)
    np.testing.assert_allclose(sph[0], 1.0, atol=1e-15, rtol=0.0)
    np.testing.assert_allclose(sph[1], 0.0, atol=1e-15, rtol=0.0)
    np.testing.assert_allclose(sph[2], 0.0, atol=1e-15, rtol=0.0)

    sph = hx.StoC(1.0, 0.0, 0.0)
    np.testing.assert_allclose(sph[0], 1.0, atol=1e-15, rtol=0.0)
    np.testing.assert_allclose(sph[1], 0.0, atol=1e-15, rtol=0.0)
    np.testing.assert_allclose(sph[2], 0.0, atol=1e-15, rtol=0.0)


    # x = 0., y = 1., z = 0.
    sph = hx.CtoS(0.0, 1.0, 0.0)
    np.testing.assert_allclose(sph[0], 1.0, atol=1e-15, rtol=0.0)
    np.testing.assert_allclose(sph[1], 0.0, atol=1e-15, rtol=0.0)
    np.testing.assert_allclose(sph[2], 90.0, atol=1e-15, rtol=0.0)

    car = hx.StoC(1.0, 0.0, 90.0)
    np.testing.assert_allclose(0., car[0], atol=1e-15, rtol=0.0)
    np.testing.assert_allclose(1., car[1], atol=1e-15, rtol=0.0)
    np.testing.assert_allclose(0., car[2], atol=1e-15, rtol=0.0)


    # x = 0., y = 0., z = 1.
    sph = hx.CtoS(0.0, 0.0, 1.0)
    np.testing.assert_allclose(sph[0], 1.0, atol=1e-15, rtol=0.0)
    np.testing.assert_allclose(sph[1], 90.0, atol=1e-15, rtol=0.0)
    np.testing.assert_allclose(sph[2], 0.0, atol=1e-15, rtol=0.0)

    car = hx.StoC(1.0, 90.0, 0.0)
    np.testing.assert_allclose(0.0, car[0], atol=1e-15, rtol=0.0)
    np.testing.assert_allclose(0.0, car[1], atol=1e-15, rtol=0.0)
    np.testing.assert_allclose(1.0, car[2], atol=1e-15, rtol=0.0)


def test_r_equal_zero():
    car0 = np.array([0., 0., 0.])
    with pytest.raises(AssertionError):
        sph = hx.CtoS(*car0)


def test_random_values():

    rng = np.random.default_rng(1)
    car0 = -1.0 + 2.0*rng.random(size=(1000, 3))
    sph = hx.CtoS(car0[:,0], car0[:,1], car0[:,2])
    car = hx.StoC(sph[0], sph[1], sph[2])

    for i in range(3):
        np.testing.assert_allclose(car0[:,i], car[i], atol=1e-12, rtol=0.0)


if __name__ == '__main__':
    # Could use: https://stackoverflow.com/a/28643737 to call all functions.
    test_basic()
    test_simple_values()
    test_r_equal_zero()
    test_random_values()

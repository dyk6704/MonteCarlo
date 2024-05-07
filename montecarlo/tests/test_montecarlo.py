"""
Unit and regression test for the montecarlo package.
"""

# Import package, test suite, and other packages as needed
import montecarlo
import sys
import numpy as np
import pytest


def test_montecarlo_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "montecarlo" in sys.modules

def test_Bistring():
    N = 9
    bs = montecarlo.BitString(N)
    assert(len(bs) == N)

    bs.flip_site(3)
    bs.flip_site(4)
    bs.flip_site(4)
    bs.flip_site(8)

    assert(all(bs.config == [0, 0, 0, 1, 0, 0, 0, 0, 1]))

    bs2 = montecarlo.BitString(12)
    bs2.set_config([1,0,0,0,0,1,0,0,1,0,1,0])
    assert(all(bs2.config == [1,0,0,0,0,1,0,0,1,0,1,0]))

    assert(bs2.on() == 4)
    assert(bs2.off() == 8)

    assert(bs2.int() == 2122)

    bs3 = montecarlo.BitString(21)
    bs3.set_int_config(5000)

    assert(all(bs3.config == [0,0,0,0,0,0,0,0,1,0,0,1,1,1,0,0,0,1,0,0,0]))


def test_average_values():
    N = 6
    Jval = 2.0
    T = 2.0

    mu = [.1 for i in range(N)]
    J = []
    for i in range(N):
        J.append([((i+1) % N, Jval), ((i-1)%N, Jval)])

    ham = montecarlo.IsingHamiltonian(J=J, mu=mu)

    E, M, HC, MS = ham.compute_average_values(T)


    assert(np.isclose(E,  -11.95991923))
    assert(np.isclose(M,   -0.00000000))
    assert(np.isclose(HC,   0.31925472))
    assert(np.isclose(MS,   0.01202961))


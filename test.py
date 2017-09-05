#!/usr/bin/env python

from __future__ import print_function

import numpy as np
from pymopac import get_energy
from pymopac import get_gradient

if __name__ == "__main__":

    atomtypes = ['O', 'H', 'H']
    coordinates = np.array(
            [[  0.000000,  0.000000,  0.000000], 
             [ -0.926021, -0.036279,  0.354716], 
             [  0.354364,  0.731748,  0.541187]])

    e = get_energy(coordinates, atomtypes, method="PM6-D3")
    print("E =", e)

    e, g = get_gradient(coordinates, atomtypes, method="PM6-D3")
    print("E =", e)
    print("G =", g)


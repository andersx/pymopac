from __future__ import print_function
import os
import uuid
import numpy as np

# Default MOPAC2016.exe path
__EXECUTABLE__ = "/opt/mopac/MOPAC2016.exe"

# Default method for PyMOPAC
__DEFAULT_METHOD__ = "PM7"

COMPONENT_INDEX = dict()
COMPONENT_INDEX["X"] = 0
COMPONENT_INDEX["Y"] = 1
COMPONENT_INDEX["Z"] = 2

def execute_mopac(filename):

    cmd = "LD_LIBRARY_PATH=/opt/mopac:$LD_LIBRARY_PATH OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 %s %s 2> /dev/null" % \
            (__EXECUTABLE__, filename)

    os.system(cmd)


def write_energy(coords, species, filename, method=__DEFAULT_METHOD__):
    """ Writes an input file to calculate the
        OM2 energy and gradient with MNDO99.

        Coords: Numpy array of coordinates in angstrom
        Species: list of names of each elements
        filename: Filename of the inputfile
    """

    content = "%s NOSYM 1SCF\n\n" % method

    for i in range(len(species)):

        content += "\n%-2s %20.10f 0 %20.10f 0 %20.10f 0 " % \
            (species[i], coords[i,0], coords[i,1], coords[i,2])

    f = open(filename, "w")
    f.write(content)
    f.close()


def read_energy(filename):
    """ Reads the output file from a energy+gradient
        calculation and returns the energy in [kcal/mol]
        and gradient in [kcal/(mol*angstrom)]

        filename: Filename of the outputfile
    """

    f = open(filename)
    lines = f.readlines()
    f.close()
    start = -1

    energy = 1e20

    for i, line in enumerate(lines):
        if "FINAL HEAT OF FORMATION" in line:
            tokens = line.split()
            energy = float(tokens[5])
            break
    
    if (energy > 1e19):
        print("MOPAC calculation failed with bad energy")
        exit()


    return energy


def write_gradient(coords, species, filename, method="PM6"):
    """ Writes an input file to calculate the
        OM2 energy and gradient with MNDO99.

        Coords: Numpy array of coordinates in angstrom
        Species: list of names of each elements
        filename: Filename of the inputfile
    """

    content = "%s NOSYM 1SCF GRADIENTS\n\n" % method

    for i in range(len(species)):

        content += "\n%-2s %20.10f 1 %20.10f 1 %20.10f 1 " % \
            (species[i], coords[i,0], coords[i,1], coords[i,2])

    f = open(filename, "w")
    f.write(content)
    f.close()

def read_gradient(filename):
    """ Reads the output file from a energy+gradient
        calculation and returns the energy in [kcal/mol]
        and gradient in [kcal/(mol*angstrom)]

        filename: Filename of the outputfile
    """

    f = open(filename)
    lines = f.readlines()
    f.close()
    start = -1

    energy = 1e20
    atoms = -1

    for i, line in enumerate(lines):

        if "Empirical Formula" in line:
            tokens = line.split()
            atoms = int(tokens[-2])

        if "FINAL HEAT OF FORMATION" in line:
            tokens = line.split()
            energy = float(tokens[5])
            break
    
    if (energy > 1e19):
        print("QM calculation failed with bad energy")
        exit()


    for i, line in enumerate(lines):
        if "PARAMETER     ATOM    TYPE            VALUE       GRADIENT" in line:
            start = i + 1
            break

    if (start == -1) or (atoms == -1):
        print("QM calculation failed without gradient")
        exit()

    grad = np.zeros((atoms,3))


    for line in lines[start:]:
        tokens = line.split()
        if len(tokens) != 8:
            break

        atom_id = int(tokens[1]) - 1
        component = COMPONENT_INDEX[tokens[4]]

        value = float(tokens[6])

        grad[atom_id, component] = value

    
    return energy, grad

def get_energy(coords, species, method=__DEFAULT_METHOD__):
    """ Calculates the OM2 energy and gradient using MOPAC2016.

        Coords: Numpy array of coordinates in angstrom
        Species: list of names of each elements
        method: selected method (default="PM7")
    """

    prefix = uuid.uuid4()

    inp_filename = "%s.mop" % prefix
    out_filename = "%s.out" % prefix
    arc_filename = "%s.arc" % prefix

    write_energy(coords, species, inp_filename, method=method)

    execute_mopac(inp_filename)

    energy = read_energy(out_filename)

    os.remove(inp_filename)
    os.remove(out_filename)
    os.remove(arc_filename)

    return energy

def get_gradient(coords, species, method=__DEFAULT_METHOD__):
    """ Calculates the OM2 energy and gradient using MOPAC2016.

        Coords: Numpy array of coordinates in angstrom
        Species: list of names of each elements
        method: selected method (default="PM7")
    """

    prefix = uuid.uuid4()

    inp_filename = "%s.mop" % prefix
    out_filename = "%s.out" % prefix
    arc_filename = "%s.arc" % prefix

    write_gradient(coords, species, inp_filename, method=method)

    execute_mopac(inp_filename)

    energy, gradient = read_gradient(out_filename)

    os.remove(inp_filename)
    os.remove(out_filename)
    os.remove(arc_filename)

    return energy, gradient

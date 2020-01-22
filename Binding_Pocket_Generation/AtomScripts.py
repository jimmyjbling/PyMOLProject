from Binding_Pocket_Generation import FindPyMol
import sys

try:
    from pymol import finish_launching, cmd, selector, stored
    import pymol2
    import numpy as np
except ImportError:
    FindPyMol.pymol_not_found()
    sys.exit()


def select_user_selection(sele, name="sele"):
    cmd.select(name, selector.process(sele))


def fetch_pdb(code, waters=1):
    cmd.fetch(code)
    if waters:
        cmd.remove('solvent')
    return code


def select_alpha_carbons(sele, name='alpha_c'):
    cmd.select(name, selector.process(sele + " and n. ca"))

    return name


def get_atom_cords(sele):
    stored.xyz = []

    cmd.iterate_state(1, sele, "stored.xyz.append([x,y,z])")

    list_xyz = []
    for atom in stored.xyz:
        list_xyz.append(atom)

    return list_xyz


def generate_center_point(sele):
    def get_average(my_list):
        my_sum = 0
        for x in my_list:
            my_sum = my_sum + x

        return my_sum / len(my_list)
    cords = get_atom_cords(sele)

    x_values = [x[0] for x in cords]
    y_values = [y[1] for y in cords]
    z_values = [z[2] for z in cords]

    x_average = get_average(x_values)
    y_average = get_average(y_values)
    z_average = get_average(z_values)

    center_point = [x_average, y_average, z_average]

    return center_point


def get_atoms(sele='polymer'):
    polymer = cmd.get_model(sele)  # this gets the protein's data from pdb
    atoms = [x for x in polymer.atom]  # gets all atom in protein. These are pyChem atom objects
    return atoms


def get_distance(point1, point2):
    x1 = point1[0]
    x2 = point2[0]
    y1 = point1[1]
    y2 = point2[1]
    z1 = point1[2]
    z2 = point2[2]
    distance = ((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2) ** 0.5
    return distance

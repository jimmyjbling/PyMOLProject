from Binding_Pocket_Generation import FindPyMol
import sys

try:
    from pymol import finish_launching, cmd, selector, stored
    import pymol2
    import numpy as np
except ImportError:
    FindPyMol.pymol_not_found()
    sys.exit()


def get_atoms(sele='polymer'):
    solvent_radius = 3.0  # this is vanderwalls of solvent, this case its water
    polymer = cmd.get_model(sele)  # this gets the protein's data from pdb
    atoms = [x for x in polymer.atom]  # gets all atom in protein. These are pyChem atom objects
    natoms = len(atoms)  # counts atoms
    radius = [x.vdw for x in atoms]  # gets vanderwalls radius of each atom.

    # the following will generate 1000 random points on the surface of a unit sphere
    rand1 = np.random.rand(1000)
    rand2 = np.random.rand(1000)

    ab_points = np.empty((1000, 3))

    z = rand1  # generates random plane of z
    theta = 2 * np.pi * rand2  # generates random angle theta

    h = np.sqrt(1 - (z**2))  # gets unknown side length

    x = h * np.cos(theta)  # gets x coordinate from angle and side length
    y = h * np.sin(theta)  # gets y coordinate from angle and side length

    # sets points
    ab_points[:, 0] = x
    ab_points[:, 1] = y
    ab_points[:, 2] = z

    # gets locations of atoms relative from each other
    coords = [x.coord for x in atoms]

    def get_point_distance(point1, point2):
        x1 = point1[0]
        x2 = point2[0]
        y1 = point1[1]
        y2 = point2[1]
        z1 = point1[2]
        z2 = point2[2]
        distance = ((x2 - x1)**2 + (y2-y1)**2 + (z2-z1)**2)**0.5
        # print(distance)
        # print(point1, point2)
        return distance

    # nearby_atoms = np.zeros((natoms, natoms))

    cut_off = 5.0 * max(radius)  # shortens program run time by ignoring atoms far way
    atoms_with_sasa = []  # hold atoms that are physically making up the SASA of the protein

    """
    SASA stands for solvent accessible surface area. This means rather than taking just the vanderwalls radius of 
    the atoms into account, it making the surface that is the limit for how far the solvent molecules can travel
    in the protein. 
    """

    for i in range(natoms):
        atom_rad = radius[i] + solvent_radius
        check_points = (ab_points * atom_rad) + coords[i]
        non_impacted_points = check_points.tolist()
        for j in range(natoms):
            dist = get_point_distance(coords[i], coords[j])
            if (dist <= cut_off) and (dist != 0):
                tmp = non_impacted_points.copy()
                for point in non_impacted_points:
                    if get_point_distance(point, coords[j]) <= radius[j]:
                        tmp.remove(point)
                non_impacted_points = tmp.copy()
        if len(non_impacted_points) > 0:
            atoms_with_sasa.append(atoms[i])

    with open("test.txt", "w") as file:
        for atom in atoms_with_sasa:
            file.write(atom.id)
    print([atom.id for atom in atoms_with_sasa])


if __name__ == '__main__':
    #finish_launching(['pymol', '-q'])
    cmd.fetch("1azm")
    print(get_atoms())

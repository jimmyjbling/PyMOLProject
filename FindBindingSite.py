"""
SASA stands for solvent accessible surface area. This means rather than taking just the vanderwalls radius of
the atoms into account, it making the surface that is the limit for how far the solvent molecules can travel
in the protein.
"""

import FindPyMol
import sys
import AtomScripts
import DockingBox

try:
    from pymol import finish_launching, cmd, selector, stored
    import pymol2
    import numpy as np
except ImportError:
    FindPyMol.pymol_not_found()
    sys.exit()


def shrake_and_rupley(atoms, box, npoints=500):
    """
    Uses the common shrake and rupley algroithum to find the SASA of the protien inside the box area of the protien
    Will make sphere around each atom with npoints and check for other atom impacts
    :param atoms: List of atoms an their info and chempy objects
    :param box: A box object defining a rough area of the binding site
    :param npoints: Set the number of points on each sphere, higher is more accurate but slower
    :return: NONE
    """
    box_atoms = box.atoms_in_box
    box_radius = [x.vdw for x in box_atoms]
    box_coords = [x.coord for x in box_atoms]
    box_num_atoms = len(box_atoms)

    solvent_radius = 1.4  # this is vanderwalls of solvent, this case its water

    # the following will generate 1000 random points on the surface of a unit sphere
    rand1 = np.random.uniform(-1, 1, size=npoints)
    rand2 = np.random.uniform(-1, 1, size=npoints)

    z = rand1  # generates random plane of z
    theta = 2 * np.pi * rand2  # generates random angle theta

    h = np.sqrt(1 - (z**2))  # gets unknown side length

    x = h * np.cos(theta)  # gets x coordinate from angle and side length
    y = h * np.sin(theta)  # gets y coordinate from angle and side length

    # sets points into coord matrix
    ab_points = np.empty((npoints, 3))
    ab_points[:, 0] = x
    ab_points[:, 1] = y
    ab_points[:, 2] = z

    atoms_with_sasa = []  # hold atoms that are physically making up the SASA of the protein

    # first loop will cycle threw each atom in the passes docking box
    for i in range(box_num_atoms):
        atom_rad = box_radius[i] + solvent_radius  # adds solvent radius to atom radius
        check_points = (ab_points * atom_rad) + box_coords[i]  # sets the random sphere points to atom coords
        list_points = check_points.tolist()  # converts to list for easier use

        # finds only atoms close to the box atoms and puts into list
        nearby_atoms = [x for x in atoms if AtomScripts.get_distance(x.coord, box_coords[i]) < 6]

        impact_test_point = []  # holds number of test points impacted
        # cycles throw all the test point of each atom
        for list_point in list_points:
            impact = False
            # cycles throw each nearby atom for each point
            for nearby in nearby_atoms:
                dist = AtomScripts.get_distance(nearby.coord, list_point)  # finds distance for nearby to point
                # if the point is within the nearby radius + solvent radius it is impacted
                if dist < nearby.vdw + solvent_radius:
                    impact = True
            if impact:
                impact_test_point.append(list_point)  # adds to holder list
        # adds atom that are not 100% impacted to sasa list
        if len(impact_test_point) != npoints:
            atoms_with_sasa.append(box_atoms[i])

    # following will select (in pymol) the atoms and residues of those atom of the atoms in the sasa list
    atom_ids = [atom.id for atom in atoms_with_sasa]

    cmd.select("sasa_atoms", "ID %d" % atom_ids[0])

    for atom_id in atom_ids:
        cmd.select("sasa_atoms", "sasa_atoms or ID %d" % atom_id)

    cmd.select("sasa_residues", "byres sasa_atoms")


def get_first_layer(threshold=3):
    """
    Finds all atoms in within threshold of each ligand atom and selects them
    :param threshold: angstrom threshold of each search
    :return: NONE
    """
    ligand = AtomScripts.get_atoms('ligand')
    polymer = AtomScripts.get_atoms('polymer')

    impacted_atoms = []

    for ligand_atom in ligand:
        for polymer_atom in polymer:
            dist = AtomScripts.get_distance(ligand_atom.coord, polymer_atom.coord) - polymer_atom.vdw
            if dist < threshold:
                impacted_atoms.append(polymer_atom)

    if impacted_atoms:
        impacted_atom_ids = [y.id for y in impacted_atoms]

        cmd.select("first_layer_atoms", "ID %d" % impacted_atom_ids[0])

        for atom_id in impacted_atom_ids:
            cmd.select("first_layer_atoms", "first_layer_atoms or ID %d" % atom_id)

        cmd.select("first_layer_res", "byres first_layer_atoms")


def combo(atoms, box, npoints=500, threshold=3):
    """
    Combines methods for more accurate selection
    :param atoms: list of chempy atoms
    :param box: box object of binding site
    :param npoints: number of points for selection
    :param threshold: angrstom selection threshold
    :return: NONE
    """
    shrake_and_rupley(atoms, box, npoints)
    get_first_layer(threshold)
    cmd.select("combo", "sasa_residues and first_layer_res")

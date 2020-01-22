from Binding_Pocket_Generation import FindPyMol
import sys
from collections import Counter

try:
    from pymol import finish_launching, cmd, selector, stored
    import pymol2
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
except ImportError:
    FindPyMol.pymol_not_found()
    sys.exit()


class Box:

    def __init__(self, x_buff=9.0, y_buff=9.0, z_buff=9.0, center_point=None):
        self.atoms_in_box = None
        self.residue_atoms = None

        self.x_center = center_point[0]
        self.y_center = center_point[1]
        self.z_center = center_point[2]

        self.x_buff = x_buff
        self.y_buff = y_buff
        self.z_buff = z_buff

        self.x_axis = [self.x_center - x_buff, self.x_center + x_buff]
        self.y_axis = [self.y_center - y_buff, self.y_center + y_buff]
        self.z_axis = [self.z_center - z_buff, self.z_center + z_buff]

    def update_box(self):
        x_buff = self.x_buff
        y_buff = self.y_buff
        z_buff = self.z_buff

        self.x_axis = [self.x_center - x_buff, self.x_center + x_buff]
        self.y_axis = [self.y_center - y_buff, self.y_center + y_buff]
        self.z_axis = [self.z_center - z_buff, self.z_center + z_buff]

    def set_box_center(self, sele="ligand", center_point=None):
        if center_point is None:
            center_point = generate_center_point(sele)

        self.x_center = center_point[0]
        self.y_center = center_point[1]
        self.z_center = center_point[2]

        self.update_box()

    def set_box_buff(self, x_buff=None, y_buff=None, z_buff=None):
        if x_buff:
            self.x_buff = x_buff
        if y_buff:
            self.y_buff = y_buff
        if z_buff:
            self.z_buff = z_buff

    def get_atoms_in_box(self, polymer=None):
        if not polymer:
            polymer = cmd.get_model('polymer')

        x_min = self.x_axis[0]
        x_max = self.x_axis[1]
        y_min = self.y_axis[0]
        y_max = self.y_axis[1]
        z_min = self.z_axis[0]
        z_max = self.z_axis[1]

        atoms = [x for x in polymer.atom if x_min <= x.coord[0] <= x_max and
                 y_min <= x.coord[1] <= y_max and
                 z_min <= x.coord[2] <= z_max]

        atom_ids = [y.id for y in atoms]

        cmd.select("binding_atoms", "ID %d" % atom_ids[0])

        for atom_id in atom_ids:
            cmd.select("binding_atoms", "binding_atoms or ID %d" % atom_id)

        self.atoms_in_box = atoms

    def get_residues_inside(self, polymer=None):
        if not self.atoms_in_box:
            self.get_atoms_in_box(polymer)

        cmd.select("box_residues", "byres " + "binding_atoms")

        box_residues = cmd.get_model("box_residues")

        self.residue_atoms = [x for x in box_residues.atom]


class PyMOLWindow:
    window_tracker = 1
    PyMols = dict()

    def __init__(self, pdb_code):
        self.pdb_code = pdb_code
        self.name = "pm" + str(PyMOLWindow.window_tracker)
        PyMOLWindow.window_tracker = PyMOLWindow.window_tracker + 1
        PyMOLWindow.PyMols[self.name] = pymol2.PyMOL('custom')
        self.pymol_object = PyMOLWindow.PyMols[self.name]

    def start(self):
        PyMOLWindow.PyMols[self.name].start()

    def fetch_pdb(self, waters=1):
        self.pymol_object.cmd.fetch(self.pdb_code)
        if waters:
            self.pymol_object.cmd.remove('solvent')

    def save_image(self, filepath, x=800, y=800, dpi=-1, ray=0):
        self.pymol_object.cmd.png(filepath, x, y, dpi, ray)

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


def get_average(my_list):
    my_sum = 0
    for x in my_list:
        my_sum = my_sum + x

    return my_sum / len(my_list)


def generate_center_point(sele):
    cords = get_atom_cords(sele)

    x_values = [x[0] for x in cords]
    y_values = [y[1] for y in cords]
    z_values = [z[2] for z in cords]

    x_average = get_average(x_values)
    y_average = get_average(y_values)
    z_average = get_average(z_values)
    
    center_point = [x_average, y_average, z_average]
    
    return center_point


def find_surface_atoms(sele="polymer"):
    cmd.create("tmp", sele, zoom=0)

    # cmd.set("dot_solvent", 1, tmpName)

    cmd.get_area(selection="tmp", load_b=1)
    
    cmd.remove("tmp" + " and b < " + "1")

    cmd.select("surface_atoms", "(" + sele + ") in " + "tmp")

    cmd.delete("tmp")

    return "surface_atoms"


def find_surface_residues(sele="all"):
    sel_name = find_surface_atoms(sele)

    cmd.select("exposed_residues", "byres " + sel_name)
    
    return "exposed_residues"


def get_distance(point1, point2):
    x1 = point1[0]
    x2 = point2[0]
    y1 = point1[1]
    y2 = point2[1]
    z1 = point1[2]
    z2 = point2[2]
    distance = ((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2) ** 0.5
    return distance
     

def get_atoms(sele='polymer'):
    polymer = cmd.get_model(sele)  # this gets the protein's data from pdb
    atoms = [x for x in polymer.atom]  # gets all atom in protein. These are pyChem atom objects
    return atoms


def do_thing(atoms, box):
    box_atoms = box.atoms_in_box
    nbox_atoms = len(box_atoms)
    box_radius = [x.vdw for x in box_atoms]
    natoms = len(atoms)  # counts atoms
    radius = [x.vdw for x in atoms]  # gets vanderwalls radius of each atom.
    solvent_radius = 1.4  # this is vanderwalls of solvent, this case its water
    # the following will generate 1000 random points on the surface of a unit sphere
    rand1 = np.random.uniform(-1, 1, size=(1000))
    rand2 = np.random.uniform(-1, 1, size=(1000))

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
    box_coords = [x.coord for x in box_atoms]
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
    print(ab_points)

    ab_points = ab_points * (solvent_radius + box_radius[0])

    x=list(ab_points[:,0])
    y=list(ab_points[:,1])
    z=list(ab_points[:,2])

    print(x,y,z)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(x,y,z, c="r", marker="o")
    plt.show()
    """

    """
    SASA stands for solvent accessible surface area. This means rather than taking just the vanderwalls radius of 
    the atoms into account, it making the surface that is the limit for how far the solvent molecules can travel
    in the protein. 
    """
    for i in range(nbox_atoms):
        atom_rad = box_radius[i] + solvent_radius
        check_points = (ab_points * atom_rad) + box_coords[i]
        list_points = check_points.tolist()

        nearby_atoms = [x for x in atoms if get_distance(x.coord, box_coords[i]) < 6]

        impact_test_point = []
        for list_point in list_points:
            impact = False
            for nearby in nearby_atoms:
                dist = get_distance(nearby.coord, list_point)
                if dist < nearby.vdw + solvent_radius:
                    impact = True
            if impact:
                impact_test_point.append(list_point)
        if len(impact_test_point) != 1000:
            atoms_with_sasa.append(box_atoms[i])

    with open("test.txt", "w") as file:
        for atom in atoms_with_sasa:
            string = str(atom.id) + "\n"
            file.write(str(string))
    print([atom.id for atom in atoms_with_sasa])

    atom_ids = [atom.id for atom in atoms_with_sasa]

    cmd.select("sasa_atoms", "ID %d" % atom_ids[0])

    for atom_id in atom_ids:
        cmd.select("sasa_atoms", "sasa_atoms or ID %d" % atom_id)

    cmd.select("sasa_residues", "byres sasa_atoms")


def get_first_layer(threshold=3):
    ligand = get_atoms('ligand')
    polymer = get_atoms('polymer')

    impacted_atoms = []

    for ligand_atom in ligand:
        for polymer_atom in polymer:
            dist = get_distance(ligand_atom.coord, polymer_atom.coord) - polymer_atom.vdw
            if dist < threshold:
                impacted_atoms.append(polymer_atom)

    if impacted_atoms:
        impacted_atom_ids = [y.id for y in impacted_atoms]

        cmd.select("first_layer_atoms", "ID %d" % impacted_atom_ids[0])

        for atom_id in impacted_atom_ids:
            cmd.select("first_layer_atoms", "first_layer_atoms or ID %d" % atom_id)

        cmd.select("first_layer_res", "byres first_layer_atoms")


def get_first_layer_fit(threshold=3):
    ligand = get_atoms('ligand')
    polymer = get_atoms('polymer')

    impacted_atoms = []

    for ligand_atom in ligand:
        for polymer_atom in polymer:
            dist = get_distance(ligand_atom.coord, polymer_atom.coord) - polymer_atom.vdw
            if dist < threshold:
                impacted_atoms.append(polymer_atom)

    if impacted_atoms:
        impacted_atom_ids = [y.id for y in impacted_atoms]

        counted_impacted_atoms = Counter(impacted_atom_ids)

        fit_ids = [key for key, value in counted_impacted_atoms if value > 2]

        cmd.select("first_layer_atoms", "ID %d" % fit_ids[0])

        for atom_id in fit_ids:
            cmd.select("first_layer_atoms", "first_layer_atoms or ID %d" % atom_id)

        cmd.select("first_layer_res", "byres first_layer_atoms")


def b_load_fix():
    cmd.create("tmp", "polymer", zoom=0)

    cmd.get_area(selection="tmp", load_b=1)

    cmd.set("dot_solvent", "on")

    b_loaded_atoms = get_atoms("tmp")

    sumb = 0

    for atom in b_loaded_atoms:
        sumb = sumb + atom.b
        print(sumb)

    cmd.remove("tmp" + " and b < " + "1")

    cmd.select("surface_atoms", "polymer in tmp")

    cmd.delete("tmp")

    return "surface_atoms"


def sasa_per_residue():
    ligand = get_atoms('organic') + get_atoms('inorganic')
    polymer = get_atoms('polymer')
    ca_atoms = get_atoms('n. ca')

    stored.residues = []
    cmd.iterate_state(1, "n. ca", "stored.residues.append(resi)")

    areas = []
    for residues in stored.residues:
        areas.append(cmd.get_area("resi %s" % residues))

    print(areas)


if __name__ == '__main__':
    """
    finish_launching(['pymol', '-q'])  # open pymol in quiet mode
    pdbCode = fetch_pdb('1azm')
    # caName = select_alpha_carbons(pdbCode)
    # cords = get_atom_cords(caName)
    select_user_selection('organic', 'ligand')
    centerPointLigand = generate_center_point('ligand')
    # surfaceResidues = find_surface_residues()
    myBox = Box(x_buff=6.5,y_buff=6.5,z_buff=6.5,center_point=centerPointLigand)
    myBox.get_residues_inside()
    atoms = get_atoms()

    do_thing(atoms, myBox)

    #b_load_fix()

    #sasa_per_residue()

    get_first_layer(5)

    """
    windows = [PyMOLWindow('1azm'), PyMOLWindow('2k1d')]

    for window in windows:
        window.start()
        window.fetch_pdb()
        window.save_image("c:\\Users\\still\\downloads\\" + window.pdb_code)


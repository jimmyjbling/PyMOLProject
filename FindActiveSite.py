import FindPyMol
import sys

try:
    from pymol import finish_launching, cmd, selector, stored
    import pymol2
except ImportError:
    FindPyMol.pymol_not_found()
    sys.exit()


class Box:

    def __init__(self, x_buff=9, y_buff=9, z_buff=9, center_point=None):
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


def select_user_selection(sele, name="sele"):
    cmd.select(name, selector.process(sele))

"""
class PyMOLWindow:
    window_tracker = 1

    def __init__(self, pdb_code, mode='-q'):
        self.pdb_code = pdb_code
        self.name = "pm" + str(PyMOLWindow.window_tracker)

        finish_launching(['pymol', mode])

    def fetch_pdb(self, waters=1):
        cmd.fetch(self.pdb_code)
        if waters:
            cmd.remove('solvent')
"""


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


def get_point_distance(point1, point2):
    x1 = point1[0]
    x2 = point2[0]
    y1 = point1[1]
    y2 = point2[1]
    z1 = point1[2]
    z2 = point2[2]
    distance = ((x2**2 - x1**1) + (y2**2 - y1**2) + (z2**2 - z1**2))**0.5
    return distance
     

if __name__ == '__main__':
    finish_launching(['pymol', '-q'])  # open pymol in quiet mode
    pdbCode = fetch_pdb('1azm')
    # caName = select_alpha_carbons(pdbCode)
    # cords = get_atom_cords(caName)
    select_user_selection('organic', 'ligand')
    centerPointLigand = generate_center_point('ligand')
    # surfaceResidues = find_surface_residues()
    myBox = Box(center_point=centerPointLigand)
    myBox.get_residues_inside()

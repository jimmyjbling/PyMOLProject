import FindPyMol
import sys

try:
    from pymol import finish_launching, cmd, selector, stored
except ImportError:
    FindPyMol.pymol_not_found()
    sys.exit()


def fetch_pdb(code, waters=1):
    cmd.fetch(code)
    if waters:
        cmd.remove('solvent')
        
    return code


def select_user_selection(sele, name="sele"):
    cmd.select(name, selector.process(sele))
    
    return name


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


def set_box(sele="ligand", center_point=None):
    if center_point is None:
        center_point = generate_center_point(sele)
    
    xcenter = center_point[0]
    ycenter = center_point[1]
    zcenter = center_point[2]
    
    xbuff = 9
    ybuff = 9
    zbuff = 9
    
    x_axis = [xcenter - xbuff, xcenter + xbuff]
    y_axis = [ycenter - ybuff, ycenter + ybuff]
    z_axis = [zcenter - zbuff, zcenter + zbuff]
    
    box = [x_axis, y_axis, z_axis]
    
    return box


def get_atoms_in_box(box):
    protien = cmd.get_model('polymer')
    
    xmin = box[0][0]
    xmax = box[0][1]
    ymin = box[1][0]
    ymax = box[1][1]
    zmin = box[2][0]
    zmax = box[2][1]
    
    atoms = [x for x in protien.atom if xmin <= x.coord[0] <= xmax and
             ymin <= x.coord[1] <= ymax and
             zmin <= x.coord[2] <= zmax]
    
    atom_ids = [y.id for y in atoms]
    
    binding_atoms = cmd.get_unused_name("binding_atoms")
    
    cmd.select("binding_atoms", "ID " + binding_atoms[0])
    
    for atom_id in atom_ids:
        cmd.select("binding_atoms", "binding_atoms or ID %d" % atom_id)
    
    return "binding_atoms"


def get_residues_in_box(box):
    binding_atoms = get_atoms_in_box(box)
    
    cmd.select("box_residues", "byres " + binding_atoms)
    
    return "box_residues"
     

if __name__ == '__main__':
    finish_launching(['pymol', '-q'])  # open pymol in quiet mode
    pdbCode = fetch_pdb('1azm')
    # caName = select_alpha_carbons(pdbCode)
    # cords = get_atom_cords(caName)
    select_user_selection('organic', 'ligand')
    centerPointLigand = generate_center_point('ligand')
    # surfaceResidues = find_surface_residues()
    myBox = set_box()
    get_residues_in_box(myBox)

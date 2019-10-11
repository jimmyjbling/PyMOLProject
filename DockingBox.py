import FindPyMol
import AtomScripts
import sys

try:
    from pymol import finish_launching, cmd, selector, stored
    import pymol2
    import numpy as np
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
            center_point = AtomScripts.generate_center_point(sele)

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

import FindPyMol
import sys
import AtomScripts
import DockingBox
import FindBindingSite

try:
    from pymol import finish_launching, cmd, selector, stored
    import pymol2
    import numpy as np
except ImportError:
    FindPyMol.pymol_not_found()
    sys.exit()

if __name__ == '__main__':
    finish_launching(['pymol', '-q'])  # open pymol in quiet mode
    pdbCode = AtomScripts.fetch_pdb('1azm') #4hla looks funny?
    AtomScripts.select_user_selection('organic', 'ligand')
    centerPointLigand = AtomScripts.generate_center_point('ligand')
    myBox = DockingBox.Box(x_buff=6.5, y_buff=6.5, z_buff=6.5, center_point=centerPointLigand)
    myBox.get_residues_inside()
    atoms = AtomScripts.get_atoms()
    FindBindingSite.combo(atoms, myBox, threshold=5)

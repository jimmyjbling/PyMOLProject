"""
@author James Wellnitz
@email wellnitz.james@gmail.com
This program is designed help set up the usage oof python scripts for PyMOL
It is meant as a learning device and not used for research or commercial purposes
"""

# All of the following imports are distributed with python3 by default
import sys
import time
from Binding_Pocket_Generation import FindPyMol

# following try catch is meant to initiate a trouble shooting step if pymol cannot be imported
try:
    from pymol import finish_launching, cmd
except ImportError:
    FindPyMol.pymol_not_found()
    sys.exit(1)

if __name__ == '__main__':
    code = 0

    while True:
        code = input("Enter 4 letter name of PDB protien file:  ")
        if type(code) is not str:
            print("enter a valid code string")
        else:
            break

    finish_launching(['pymol', '-q'])  # open pymol in quiet mode
    time.sleep(2)
    cmd.fetch(code)  # gets protein code to load in
    cmd.spectrum()
    cmd.mset("1x180")
    cmd.movie.roll(1, 180, 1)
    cmd.movie.produce("movietest", mode="draw")

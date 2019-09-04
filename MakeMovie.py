"""
@author James Wellnitz
@email wellnitz.james@gmail.com
This program is designed help set up the usage oof python scripts for PyMOL
It is meant as a learning device and not used for research purposes
"""

import itertools
import sys
import time
import os
import threading


def search_for_pymol():
    search_complete = False

    def generate_spinner():
        spinner = itertools.cycle(['-', '\\', '|', '/'])
        while not search_complete:
            sys.stdout.write(next(spinner))
            sys.stdout.flush()
            time.sleep(0.3)
            sys.stdout.write('\b')

    threading.Thread(target=generate_spinner).start()

    possible_locations = []
    
    def get_pymol_paths(name, path):
        for root, dirs, files in os.walk(path):
            if name in files or name in dirs:
                possible_locations.append(os.path.join(root, name))
        return possible_locations

    get_pymol_paths("pymolhttpd.py", "c:/")

    python_locations = []
    for possible_location in possible_locations:
        if "pkgs" in possible_location or "web" in possible_location:
            possible_locations.remove(possible_location)
        else:
            index = possible_location.find("site-packages")
            python_location = possible_location[:index - 4] + "python.exe"
            python_locations.append(python_location)
    search_complete = True
    time.sleep(0.5)
    return python_locations


try:
    import pymol
except ImportError:
    print("pymol package is not associated with this python interpreter\nAttempting to locate ")
    locations = search_for_pymol()
    if locations is not None:
        print("Located python distribution with pymol. Currently running python at %s which lacks pymol. Run program "
              "with interpreter at one of the following locations instead:" % sys.executable)
        for location in locations:
            print(location)
    else:
        print("Could not find installation of pymol")
        if "anaconda" in sys.executable:
            print("Current python is an anaconda environment. Open environment manager and run following command to "
                  "install pymol into desired anaconda environment:\nconda install -c schrodinger pymol")
        else:
            print("Consider the following potential issues:\npymol might not be install on your device. Install at "
                  "https://pymol.org/2/\npymol might be installed outside the c drive. Look for pymol and the python "
                  "that comes with it outside c drive\n\nNOTE: Pymol cannot be installed with pip (easily). Either "
                  "download pymol and use its built in python or download anaconda, set up a conda python envirment "
                  "and get the pymol package there.")

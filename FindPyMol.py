import itertools
import sys
import time
import os
import threading


def search_for_pymol():
    """
    Initiates a search for a python interpreter with pymol installed
    :return: list of appropriate python.exe location in computer
    """
    search_complete = False

    def generate_spinner():
        """
        Generates a spinner so user know the program is still running during search
        :return: None
        """
        spinner = itertools.cycle(['-', '\\', '|', '/'])
        while not search_complete:
            sys.stdout.write(next(spinner))  # next spinner in cycle
            sys.stdout.flush()  # flushes buffer
            time.sleep(0.3)  # waits so it doesnt overload
            sys.stdout.write('\b')  # removes spinner text so it can be replaced

    threading.Thread(target=generate_spinner).start() # starts spinner on different thread

    possible_locations = []  # holds uncleaned pymol locations

    def get_pymol_paths(name, path):
        """
        walks the users file system starting in c drive to find file/directory and stores locations
        :param name: name of file or directory to search for
        :param path: starting directory for search
        :return: list of all full path locations with file/directory
        """
        for root, dirs, files in os.walk(path):
            if name in files or name in dirs:
                possible_locations.append(os.path.join(root, name))  # collects full path name
        return possible_locations

    get_pymol_paths("pymolhttpd.py", "c:/")  # pymolhttpd.py is a file unique to pymol package

    python_locations = []  # hold cleaned python.exe locations

    # following will select only locations of pymol with a python interpreter associated with them and change them
    # into file location of said python interpreter
    for possible_location in possible_locations:
        if "pkgs" in possible_location or "web" in possible_location:  # removes locations without interpreter
            possible_locations.remove(possible_location)
        else:
            index = possible_location.find("site-packages")
            python_location = possible_location[:index - 4] + "python.exe"  # changes location to that of python.exe
            python_locations.append(python_location)
    search_complete = True  # ends spinner thread
    time.sleep(0.5)  # wait for spinner to flush and \b in last cycle
    return python_locations


def pymol_not_found():
    """
    Activates trouble shooting for pymol import error to help use figure out why pymol is not imported correctly
    :return: None
    """
    print("pymol package is not associated with this python interpreter\nAttempting to locate ")
    locations = search_for_pymol()

    # lists locations for suitable interpreter
    if locations is not None:
        print("Located python distribution with pymol. Currently running python at %s which lacks pymol. Run program "
              "with interpreter at one of the following locations instead:" % sys.executable)
        for location in locations:
            print(location)
    # attempts to give advice for finding or installing a suitable python with pymol
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

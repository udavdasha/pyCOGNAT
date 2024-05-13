"""
This module contains class <Settings> and a method for reading an input file into it.
This is used in <COGalyser>, <Alnalyser>, <pyCOGNAT> and <NyanTranslate> software

Module version: 1.1.2
"""

import sys, os, platform

def get_program_name(directory, name):
    """
    Method constructs the exact path to a program with <name>  based on the
    operation system type
    """
    curr_name = name
    if platform.system() == "Windows":
        curr_name += ".exe"
    curr_path = os.path.join(directory, curr_name)
    return (curr_name, curr_path)

class Settings:
    def __init__(self, settings_dict, req_attributes):
        for key in settings_dict.keys():
            setattr(self, key, settings_dict[key])
        self.req_attributes = req_attributes
        self.check()
        #self.work_dir = os.path.abspath(self.work_dir)

    def check(self):
        good = True
        print ("\nChecking availability of the required settings:")
        for attribute in self.req_attributes:
            try:
                print ("[....OK....]\t%s = %s" % (attribute, getattr(self, attribute)))
            except:
                print ("[NOT FOUND]\t%s" % attribute)

    def write_file(self, filename):
        settings_file = open(filename, "w")
        settings_file.write("# This is a settings file created by a user\n")
        for attribute in self.req_attributes:
            value = getattr(self, attribute)
            settings_file.write("%s = %s\n" % (attribute, value))
        settings_file.write("#\n")
        settings_file.close()

def read_settings_file(filename, req_attributes, not_path_attr = dict()):
    """
    Method reads the settings file <filename> for the program. If the file does
    not exist, empty dictionary will be created
    Returns a <Settings> class object with obtained settings as attributes,
    and <req_attributes> list set as its required attributes.
    Dictionary <not_path_attr> with names of attributes which are not
    paths to files could be given.
    """
    set = dict() # Empty dictionary for settings
    try:
        settings_file = open(filename, "r")
        for string in settings_file:
            string = string.strip()
            if len(string) == 0:
                continue
            if string[0] == "#": # Skipping comment lines
                continue
            fields = string.split("=", 1)
            fields[0] = fields[0].strip()
            fields[1] = fields[1].strip()
            if (fields[1].isdigit()) or (fields[0] in not_path_attr):
                set[fields[0]] = fields[1]
            else:
                set[fields[0]] = os.path.abspath(fields[1])
        settings_file.close()
    except:
        print ("")
        print ("[FATAL ERROR] Settings file '%s' does not exist, consider creating one for %s!" % (filename, platform.system()))
        print ("              If your settings file has another name, pass it to the main script as an argument")
        sys.exit() #FIX (version 1.1.1): if no settings file presents, further work is impossible
    settings_object = Settings(set, req_attributes)
    return settings_object
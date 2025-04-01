#!/usr/bin/env python
import sys, os, argparse

#========================================================================================
curr_version = 1.0
parser = argparse.ArgumentParser(description = 
"This script will fix the COGNAT database by removing artificial duplicates which result from \
non-unique protein ids.  \
Current version is %s" % curr_version 
)
parser.add_argument("-i", help = "Name of the COGNAT database main folder", required = True, dest = "COGNAT_dir")
myargs = parser.parse_args()

def filter_identical_strings(in_filename, replace = False, report_file = None):
    input_file = open(in_filename)
    string_already_added = dict()
    strings = list()
    for string in input_file:
        string = string.strip()
        if len(string) == 0:
            continue
        if not string in string_already_added:
            strings.append(string)
            string_already_added[string] = True
        else:
            if report_file != None:
                report_file.write("%s\t%s\n" % (in_filename, string))
    input_file.close()

    if replace:
        output_file = open(in_filename, "w")
        for string in strings:
            output_file.write("%s\n" % string)
        output_file.close()

report_file_cogs = open("report_COGs.txt", "w")
report_file_pfam = open("report_Pfam.txt", "w")
gbff_list = os.listdir(myargs.COGNAT_dir)
n = 0
for gbff in gbff_list:
    n += 1
    print ("Proceeding %i gbff folder (out of %i)" % (n, len(gbff_list)))
    gbk_list = os.listdir(os.path.join(myargs.COGNAT_dir, gbff))
    for gbk in gbk_list:
        cog_filename = os.path.join(myargs.COGNAT_dir, gbff, gbk, "c_%s.txt" % gbk)
        filter_identical_strings(cog_filename, True, report_file_cogs)
        pfam_filename = os.path.join(myargs.COGNAT_dir, gbff, gbk, "d_%s.txt" % gbk)
        filter_identical_strings(pfam_filename, True, report_file_pfam)
report_file_cogs.close()
report_file_pfam.close()
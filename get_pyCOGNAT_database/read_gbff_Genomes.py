#!/usr/bin/env python
import sys, os, argparse, re
import shutil
import udav_base, udav_read_GenBank, udav_fasta

#========================================================================================
curr_version = "1.6.1"
#========================================================================================
parser = argparse.ArgumentParser(description = 
"This script will read genomes from input directory and do their split & read quickly. \
Current version is %s" % curr_version 
)
parser.add_argument("-i", help = "Name of directory with genome files (not splitted, *.gbff)", required = True, dest = "input_dir")
parser.add_argument("-o", help = "Prefix for the output files", required = True, dest = "output")
parser.add_argument("-r", help = "Filename with the list of record ids which should not be removed", required = False, dest = "req_records")
parser.add_argument("-a", help = "Toggle this to let all splitted files remain", action = "store_true", dest = "all_records_req")
parser.add_argument("-l", help = "List of proteins which operon structure is to be obtained (under locus_tags)", required = False, dest = "req_genes")
parser.add_argument("-s", help = "Type of script for <hmmscan> to be produced ('bat' for Windows or 'bash' for Linux)", default = "bash", dest = "script_type")
parser.add_argument("-d", help = "HMM database path to write to the <hmmscan> script", required = False, dest = "database_path")
parser.add_argument("-n", help = "Toggle this to print nucleotide sequences in FASTA format", action = "store_true", dest = "print_nucl")
parser.add_argument("-p", help = "Toggle this to print protein sequences in FASTA format", action = "store_true", dest = "print_prot")
parser.add_argument("-m", help = "Toggle this to print desired sequences in multiple separate files", action = "store_true", dest = "mult_files")
parser.add_argument("--old_locus_search", help = "Use this if table with old locuses should be printed", action = "store_true", dest = "old_locus")
parser.add_argument("--huge_genomes", help = "Use this if sequence data should not be stored in memory (should use -m option)", action = "store_true", dest = "huge_genomes")
parser.add_argument("--intergene", help = "Use this if data about intergene distances should be printed for each assembly", action = "store_true", dest = "intergene")
parser.add_argument("--non_empty_records", help = "Use this to report only records contaning genes into .records file", action = "store_true", dest = "non_empty_records")
parser.add_argument("--file_to_assembly", help = "Give here assignment between filename and real assembly ID (if the latter should not be guessed from the former)", required = False, dest = "file_to_assembly")
parser.add_argument("--pseudo", help = "Toggle this to print pseudogene nucleotide sequences in a separate file", action = "store_true", dest = "print_pseudo")
myargs = parser.parse_args()

if myargs.huge_genomes and (not myargs.mult_files):
    print ("FATAL ERROR: Option --huge_genomes used without the -m option, this is incompatible!")
    sys.exit() 

required_records = dict()
required_genes = dict()
operon_output = None
table_output = None
if myargs.req_records != None:
    required_records = udav_fasta.read_req_list (myargs.req_records, "exact")
if myargs.req_genes != None:
    required_genes = udav_fasta.read_req_list (myargs.req_genes, "exact")
    operon_output = open("%s.operone.fasta" % myargs.output, "w")
    table_output = open("%s.table.operone.txt" % myargs.output, "w")
    table_output.write("#record_id\torganism\toperon_num\tprotein_id\tlocus_tag\tproduct\tbegin\tend\tdirection\n")

file_to_assembly = dict() #FIX: version 1.6.0 add
write_mode = "w"
if myargs.file_to_assembly != None:
    file_to_assembly = udav_base.read_unique_assignment(myargs.file_to_assembly)
    write_mode = "a"
    if myargs.mult_files: print("WARNING: files will be written in 'a' mode because otherwise records could be lost!")

#-------------------------------------------------- FASTA FILES WRITING --------------
hmmscan_script = None
nucleotide_output = None
protein_output = None
pseudo_output = None # FIX: added in version 1.6.1
if myargs.mult_files: #FIX: version 1.5.3 add
    hmmscan_script = open("%s_hmmscan.%s" % (myargs.output, myargs.script_type), "w")
    path_to_database = "<PATH_TO_THE_DATABASE>"
    if myargs.database_path != None:
        path_to_database = myargs.database_path
    if myargs.script_type == "bash":
        hmmscan_script.write("#!/bin/bash\n")
        hmmscan_script.write("hmm_database=%s\n" % path_to_database)
        hmmscan_script.write("hmmscan=/home/udavdasha/_soft/COG_analyzer/soft/hmmscan\n\n")
    if myargs.script_type == "bat":
        hmmscan_script.write("@echo off\n")
        hmmscan_script.write("set hmm_database=%s\n" % path_to_database)
        hmmscan_script.write("set hmmscan=D:\\UdavBackup\\_Complete_genomes\\_bioinf_soft\\hmmer-3.1b2-cygwin64\\binaries\\hmmscan.exe\n\n")
else:
    if myargs.print_nucl:
        nucleotide_output = open("%s.nucl.fasta" % myargs.output, "w")
    if myargs.print_prot:
        protein_output = open("%s.prot.fasta" % myargs.output, "w")
if myargs.print_pseudo:
    pseudo_output = open("%s.pseudo.fasta" % myargs.output, "w")
#-------------------------------------------------------------------------------------
log_file = open("%s.get_nucl.log" % myargs.output, "w")
assembly_prot_data = dict()
records_data_file = open("%s.records" % myargs.output, "w")
records_data_file.write("#This is the output of <read_gbff_Genomes.py> script\n")
records_data_file.write("#req_mark\tinput_dir\tgbff\trecord_id\tlength_bp\tgene_number\torganism\trecord_name\n") # FIX: version 1.5.9 - name of the record is printed
assembly_prot_file = open("%s.assemblies" % myargs.output, "w") #FIX 1.5.2: now number of proteins found in each assembly is calculated separately (1.5.7: printing occurs every assembly, not in the end)
assembly_prot_file.write("#gbff\tnumber_of_parts\ttotal_length_bp\tgene_number\torganism\ttaxonomy\n")
temp_dir = "%s_temp_dir" % myargs.output
if not os.path.isdir(temp_dir):
    print ("Creating temporary directory for files: %s" % temp_dir)
    os.makedirs(temp_dir)

#files = os.listdir(myargs.input_dir)
num_files = 0
for root, subdirs, files in os.walk(myargs.input_dir):
    num_files += len(files)
log_file.write("Found %i files in the input directory %s\n" % (num_files, myargs.input_dir))
n = 0

old_locus_data = dict() #FIX: version 1.5.6

dist_direct_global = {"PROT-PROT":dict(), "PROT-RNA":dict(), "RNA-RNA":dict(), "OTHER":dict(), -1:list(), -4:list()} #FIX: version 1.5.8
dist_reverse_global = {"PROT-PROT":dict(), "PROT-RNA":dict(), "RNA-RNA":dict(), "OTHER":dict(), -1:list(), -4:list()}

for root, subdirs, files in os.walk(myargs.input_dir):
    for f in files:
        fields = f.split(".")
        curr_ext = fields[len(fields) - 1]
        if f in file_to_assembly: #FIX: version 1.6.0 (assembly could not be guessed but stated explicitly)
            curr_assembly = file_to_assembly[f]
        else:
            try:
                curr_assembly = re.match("(.+)\_genomic\.gbff$", f).group(1)
            except AttributeError:
                try:
                   curr_assembly = re.match("(.+)\.gbff$", f).group(1)
                except AttributeError:
                   print ("Error in the name format of a file: '%s'" % f)
                   raise
        n += 1 
        if (curr_ext == "gbk") or (curr_ext == "gbff") or (curr_ext == "gp") or (curr_ext == "gb"):
            #---------------------------------- FASTA FILES WRITING ------------------
            curr_prot_fasta_file = None
            curr_nucl_fasta_file = None
            if myargs.mult_files: #FIX: version 1.5.3 add
                if myargs.print_prot:
                    curr_prot_fasta_file = open("%s.fasta" % curr_assembly, write_mode)
                    if myargs.script_type == "bash":
                        hmmscan_script.write("$hmmscan --tblout %s.table --domtblout %s.table.domains -o /dev/null $hmm_database %s.fasta\n" % (curr_assembly, curr_assembly, curr_assembly))
                    if myargs.script_type == "bat":
                        hmmscan_script.write("%hmmscan% --tblout {0}.table --domtblout {0}.table.domains -o to_remove %hmm_database% {0}.fasta\n".format(curr_assembly))
                if myargs.print_nucl:
                    curr_nucl_fasta_file = open("%s.nucl.fasta" % curr_assembly, write_mode)
            #-------------------------------------------------------------------------

            curr_filename = os.path.join(root, f)
            curr_folder = os.path.join(temp_dir, f)
            if not os.path.isdir(curr_folder):
                os.makedirs(curr_folder) 
            udav_base.report_info("Working with the file %s (%i out of %i)" % (f, n, num_files), log_file)
            (c_num, f_num, name_prefix, org) = udav_read_GenBank.split_gbk(curr_filename, curr_folder, curr_folder, remove_empty = False)
            miss = 0
            dist_direct = {"PROT-PROT":dict(), "PROT-RNA":dict(), "RNA-RNA":dict(), "OTHER":dict(), -1:list(), -4:list()} #FIX: version 1.5.8
            dist_reverse = {"PROT-PROT":dict(), "PROT-RNA":dict(), "RNA-RNA":dict(), "OTHER":dict(), -1:list(), -4:list()}
            for i in range(1, f_num + c_num + 1):
                #Filenames are: <name_prefix>_<i>.gbk
                do_not_remove = myargs.all_records_req # FIX: version 1.5.5 (False by default, but can be set to True)
                filename = "%s_%i.gbk" % (name_prefix, i)
                path = os.path.join(curr_folder, filename)
                if not os.path.isfile(path): # File is missing for some reason
                    miss += 1
                    continue
                try:
                    new_genome = udav_read_GenBank.Genome_file(path)
                    new_genome.obtain_data(myargs.huge_genomes, curr_nucl_fasta_file)
                    new_genome.get_intergene_dist(dist_direct, dist_reverse) # FIX: version 1.5.8
                    new_genome.get_intergene_dist(dist_direct_global, dist_reverse_global)
                    if (new_genome.record_id in required_records) or (name_prefix in required_records):
                        do_not_remove = True
                    if ((myargs.non_empty_records) and (new_genome.gene_number != 0)) or (not myargs.non_empty_records): #FIX: version 1.5.9
                        records_data_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (do_not_remove, myargs.input_dir, f, new_genome.record_id,
                                                                                  new_genome.genome_length, new_genome.gene_number, new_genome.organism, new_genome.record_name))
                    if not curr_assembly in assembly_prot_data:
                        assembly_prot_data[curr_assembly] = [0, 0, 0, new_genome.organism, new_genome.taxonomy]
                    
                    assembly_prot_data[curr_assembly][0] += 1 # number of parts in gbff
                    assembly_prot_data[curr_assembly][1] += new_genome.genome_length # total length (bp)
                    assembly_prot_data[curr_assembly][2] += new_genome.gene_number # gene number

                    #--------------------------- FASTA FILES WRITING -----------------
                    if myargs.mult_files:
                        if myargs.print_nucl and (not myargs.huge_genomes):
                            new_genome.print_nucleotide(curr_nucl_fasta_file)
                        if myargs.print_prot:
                            new_genome.print_proteins(curr_prot_fasta_file, curr_assembly)
                    else:
                        if myargs.print_nucl: 
                            new_genome.print_nucleotide(nucleotide_output)
                        if myargs.print_prot:
                            new_genome.print_proteins(protein_output, curr_assembly)
                    if myargs.print_pseudo == True:
                        new_genome.print_pseudo(pseudo_output)                     
                    #-----------------------------------------------------------------
                    if myargs.req_genes != None:
                        new_genome.get_operons(required_genes, operon_output, table_output, "locus_tag")
                    if myargs.old_locus != None:
                        new_genome.get_old_locus_data(old_locus_data)
                except IOError:
                    log_file.write("    [WARNING]: file '%s' raised an exception (IOError)!\n" % filename)
                except TypeError:
                    log_file.write("    [WARNING]: file '%s' raised an exception (TypeError)!\n" % filename)
                if not do_not_remove:
                    os.remove(path)
            #---------------------------------- FASTA FILES WRITING ------------------
            if myargs.mult_files: #FIX: version 1.5.3 add
                if myargs.print_prot:
                    curr_prot_fasta_file.close()
                if myargs.print_nucl:
                    curr_nucl_fasta_file.close()
            #-------------------------------------------------------------------------
            #--------------------------- GENE DISTANCES DATA WRITING -----------------
            if myargs.intergene:
                udav_read_GenBank.print_dist_file("%s.dist.direct" % curr_assembly, dist_direct, -100, 100)
                udav_read_GenBank.print_dist_file("%s.dist.reverse" % curr_assembly, dist_reverse, -100, 100)
            #-------------------------------------------------------------------------

            curr_data = assembly_prot_data[curr_assembly]
            assembly_prot_file.write("%s\t%i\t%i\t%i\t%s\t%s\n" % (curr_assembly, curr_data[0], curr_data[1], curr_data[2], curr_data[3], "; ".join(curr_data[4])))
            udav_base.report_info("    DONE: %i total small files (%i with no proteins), %i missing files" % (f_num + c_num, c_num, miss), log_file)
#-------------------------------------------------- FASTA FILES WRITING --------------
if myargs.mult_files: #FIX: version 1.5.3 add
    hmmscan_script.close()
else:
    if myargs.print_nucl:
        nucleotide_output.close()
    if myargs.print_prot:
        protein_output.close()
if myargs.print_pseudo:
    pseudo_output.close()
#-------------------------------------------------------------------------------------
udav_read_GenBank.print_dist_file("%s.dist.direct" % myargs.output, dist_direct_global, -100, 100)
udav_read_GenBank.print_dist_file("%s.dist.reverse" % myargs.output, dist_reverse_global, -100, 100)

if (myargs.req_records == None) and (myargs.all_records_req == False):
    shutil.rmtree(temp_dir)
log_file.close()
records_data_file.close()
assembly_prot_file.close()

if myargs.req_genes != None:
    operon_output.close()
    table_output.close()

if myargs.old_locus: #FIX: version 1.5.6 add
    old_locus_file = open("%s.old_locus.txt" % myargs.output, "w")
    old_locus_file.write("#locus_tag\tprotein_id\told_locus_tag\trecord_id\n")
    for curr_locus in old_locus_data.keys():
        old_locus_string = old_locus_data[curr_locus][0]
        protein_id = old_locus_data[curr_locus][1]
        record_id = old_locus_data[curr_locus][2]
        if old_locus_string != None:
            old_locuses = old_locus_string.split("<>")
            for old_locus in old_locuses:
                old_locus_file.write("%s\t%s\t%s\t%s\n" % (curr_locus, protein_id, old_locus, record_id))
        else:
            old_locus_file.write("%s\t%s\t%s\t%s\n" % (curr_locus, protein_id, None, record_id))
    old_locus_file.close()
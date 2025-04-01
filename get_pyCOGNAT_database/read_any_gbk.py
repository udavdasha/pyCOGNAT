#!/usr/bin/env python
import sys, os, argparse, re
import udav_read_GenBank

#========================================================================================
curr_version = "1.0"
#========================================================================================
parser = argparse.ArgumentParser(description = 
"This script will split & read any gbk file. \
Current version is %s" % curr_version 
)
parser.add_argument("-i", help = "Name of the gbk file (not splitted)", required = True, dest = "input_gbk")
parser.add_argument("-o", help = "Prefix for the output files", required = True, dest = "output")
myargs = parser.parse_args()

def read_data_from_gbk_dir(gbk_dir, output_prefix, log_file):
    print ("Reading files from the '%s' directory" % gbk_dir)
    hmmscan_script_bash = open("%s.hmmscan.bash" % output_prefix, "w")
    hmmscan_script_bash.write("#!/bin/bash\n")
    hmmscan_script_bash.write("hmm_database=...\n")
    hmmscan_script_bash.write("hmmscan=/usr/bin/hmmscan\n\n")
    hmmscan_script_bat = open("%s.hmmscan.bat" % output_prefix, "w")
    hmmscan_script_bat.write("@echo off\n")
    hmmscan_script_bat.write("set hmm_database=..\n")
    hmmscan_script_bat.write("set hmmscan=D:\\UdavBackup\\_Complete_genomes\\_bioinf_soft\\hmmer-3.1b2-cygwin64\\binaries\\hmmscan.exe\n\n")

    files = os.listdir(gbk_dir)
    i = 0
    n = 0
    i_thresh = 100
    for single_file in files:
        i += 1
        n += 1
        if i >= i_thresh:
            i = 0
            print ("Reading %i file... (out of %i)" % (n, len(files)))
        path = os.path.join(gbk_dir, single_file)
        curr_source = ".".join(single_file.split(".")[0:-1])
        try:
            new_genome = udav_read_GenBank.Genome_file(path)
            curr_nucl_fasta_file = open("%s.nucl.fasta" % curr_source, "w")
            new_genome.obtain_data(False, curr_nucl_fasta_file)
            new_genome.print_nucleotide(curr_nucl_fasta_file)
            curr_nucl_fasta_file.close()
            with_CDS = False
            for g in new_genome.genes:
                if g.gene_type == "CDS":
                    with_CDS = True
                    break
            log_file.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (with_CDS, new_genome.record_id, new_genome.genome_length, new_genome.gene_number, new_genome.organism, new_genome.record_name))
            if with_CDS:
                hmmscan_script_bash.write("$hmmscan --tblout %s.table --domtblout %s.table.domains -o /dev/null $hmm_database %s.fasta\n" % (curr_source, curr_source, curr_source))
                hmmscan_script_bat.write("%hmmscan% --tblout {0}.table --domtblout {0}.table.domains -o to_remove %hmm_database% {0}.fasta\n".format(curr_source))
                curr_prot_fasta_file = open("%s.fasta" % curr_source, "w")
                new_genome.print_proteins(curr_prot_fasta_file, "Unk")
                curr_prot_fasta_file.close()
                #new_genome.print_pseudo(pseudo_output)
        except IOError:
            log_file.write("    [WARNING]: file '%s' raised an exception (IOError)!\n" % single_file)
        except TypeError:
            log_file.write("    [WARNING]: file '%s' raised an exception (TypeError)!\n" % single_file)    
    hmmscan_script_bash.close()
    hmmscan_script_bat.close()

split_dir = "%s_split" % myargs.output
empty_dir = "%s_empty" % myargs.output
log_file = open("%s.log" % myargs.output, "w")
if (not os.path.isdir(split_dir)) and (not os.path.isdir(empty_dir)):
    os.mkdir(split_dir)
    os.mkdir(empty_dir)
    print ("Reading and splitting file '%s'..." % myargs.input_gbk)
    (c, f, name_prefix, org) = udav_read_GenBank.split_gbk(myargs.input_gbk, split_dir, empty_dir, remove_empty = False, list_under_source = True)
    print ("DONE. Total %i files with CDS and %i files without them were created!" % (c, f))
    log_file.write("# Read '%s' gbk file\n" % myargs.input_gbk)
    log_file.write("# Files with CDS data: %i\n" % f)
    log_file.write("# Files without CDS data: %i\n" % c)
log_file.write("#with_CDS\trecord_id\tlength_bp\tgene_number\torganism\trecord_name\n") # FIX: version 1.5.9 - name of the record is printed

# 1) READING FILES WITH CDS
read_data_from_gbk_dir(split_dir, myargs.output, log_file)
log_file.close()
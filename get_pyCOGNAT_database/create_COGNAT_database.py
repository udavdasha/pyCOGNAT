#!/usr/bin/env python
import sys, os, argparse
import re
import time
import udav_fasta, udav_soft

#========================================================================================
curr_version = 1.2
parser = argparse.ArgumentParser(description = 
"This script will create a COGNAT database based on the results of HMM search. \
Current version is %s" % curr_version 
)
parser.add_argument("-i", help = "Name of folder with fasta-formatted sequences ('URef' format, produced by <read_gbff_genomes.py> under -p option or <read_any_gbk.py>)", required = True, dest = "fasta_dir")
parser.add_argument("-c", help = "Folder with results of --domtblout option of <hmmscan> for COGs (should be obtained with the *_hmmscan.bash script)", required = True, dest = "COG_domtblout_dir")
parser.add_argument("-n", help = "Name of nucleotide sequences folder (produced by <read_gbff_genomes.py> under -n option or <read_any_gbk.py>)", required = True, dest = "nucl_fasta_dir")
parser.add_argument("-o", help = "Output directory for COGNAT database", required = True, dest = "output_dir")
parser.add_argument("-l", help = "Name of the log-file", required = True, dest = "log_file")
parser.add_argument("-m", help = "Which ID to use as primary: 'ID' for protein_id (DEFAULT) or 'locus' for locus_tag", required = False, default = "ID", dest = "mode")
parser.add_argument("-p", help = "(OPTIONAL) Folder with results of --domtblout option of <hmmscan> for Pfam (should be obtained with the *_hmmscan.bash script)", required = False, dest = "Pfam_domtblout_dir")
parser.add_argument("-t", help = "(OPTIONAL) Name of the TMHMM result in short format (one line per protein)", required = False, dest = "tmhmm_result")
parser.add_argument("-e", help = "E-value threshold for HMMer searches (i-Evalue) (DEFAULT = 0.1)", required = False, default = 0.1, dest = "e_value_threshold")
myargs = parser.parse_args()
myargs.e_value_threshold = float(myargs.e_value_threshold)

def get_primary_and_secondary_id(protein_object, mode):
    primary = None
    secondary = None
    if mode == "ID":
        primary = protein_object.protein_id
        secondary = protein_object.locus
    elif mode == "locus":
        primary = protein_object.locus
        secondary = protein_object.protein_id
    else:
        print ("FATAL ERROR: use either 'ID' or 'locus' as mode!")
        sys.exit()
    return (primary, secondary)

def check_input_dir_or_file(path_to_object, check_method, error_message, should_exist):
    if should_exist != check_method(path_to_object):
        print (error_message % path_to_object)    
        sys.exit()

def print_domain_hits_for_protein(protein_id, id_to_domains, protein_id_to_print, output_file):    
    if protein_id in id_to_domains:
        for domain_id in id_to_domains[protein_id].keys():
            domain_hits = id_to_domains[protein_id][domain_id].strip().split(" ")
            for hit in domain_hits:
                #hit: "%s..%s..%s..%s..%s..%s..%s" % (begin, end, evalue, score, hmm_begin, hmm_end, hmm_covered)                        
                hit_fields = hit.split("..")
                if len(hit_fields) != 7:
                    print ("FATAL ERROR: corrupt domain file, error for protein '%s'" % protein_id)
                    print (hit_fields)
                    sys.exit()
                #TM7_0005	COG2827	2.1e-08	32.0	55.0	8	67
                output_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein_id_to_print, domain_id, hit_fields[2], hit_fields[3], 
                                                                            hit_fields[6], hit_fields[0], hit_fields[1], hit_fields[4], hit_fields[5]))

check_input_dir_or_file(myargs.fasta_dir, os.path.isdir,"FATAL ERROR: Input directory '%s' with protein FASTA-files is missing", True) 
check_input_dir_or_file(myargs.COG_domtblout_dir, os.path.isdir, "FATAL ERROR: Input directory '%s' with domain tables for COG search is missing", True)
#check_input_dir_or_file(myargs.Pfam_domtblout_dir, os.path.isdir, "FATAL ERROR: Input directory '%s' with domain tables for Pfam search is missing", True)
check_input_dir_or_file(myargs.nucl_fasta_dir, os.path.isdir, "FATAL ERROR: Input directory '%s' with nucleotide FASTA-files is missing", True)
#check_input_dir_or_file(myargs.output_dir, os.path.isdir, "FATAL ERROR: Output directory '%s' was already created, please provide another name", False)
if not os.path.isdir(myargs.output_dir):
    os.mkdir(myargs.output_dir)

log_file = open(myargs.log_file, "w")
current_time = time.strftime("%a, %d %b %Y %X", time.localtime())
log_file.write("# Log file produced by the <create_COGNAT_database.py> script at: %s\n" % current_time)
log_file.write("#\n")
log_file.write("#assembly\tnum_of_gbk\tgbk_without_proteins\tprotein number\ttime spent\ttime global\n")
t_global_start = time.time()
#
tmhmm_features = dict()
if myargs.tmhmm_result != None:
    print ("Reading global TMHMM output '%s' into a dictionary and preparing it..." % myargs.tmhmm_result)
    tmhmm_features_unprepared = udav_soft.read_TMHMM_output(myargs.tmhmm_result)
    for protein_id in tmhmm_features_unprepared.keys():
        #value = [TMHMM] 13..32,52..74,87..109,119..138,145..167,172..189
        num_of_helices = len(tmhmm_features_unprepared[protein_id].split(" ")[1].split(","))
        tmhmm_features[protein_id] = num_of_helices
        #print ("%s\t%s\t%s" % (protein_id, tmhmm_features_unprepared[protein_id], num_of_helices))
    print ("DONE!")
else:
    print ("TMHMM file was not provided, no data on transmembrane helices will be put in the database!")

fasta_filenames = os.listdir(myargs.fasta_dir)
gbff_fasta_count = 0
for fasta_filename in fasta_filenames:
    prefix = re.match("(.+)\.fasta$", fasta_filename).group(1)
    if prefix == None:
        print ("WARNING: file with unproper name detected, skipping: '%s'" % fasta_filename)
        continue
    t_start = time.time()
    t_global = (t_start - t_global_start)/60
    gbff_fasta_count += 1
    print ("Reading assembly #%i (out of %i): '%s' (time running: %.1f min)" % (gbff_fasta_count, len(fasta_filenames), prefix, t_global))

    curr_assembly_dir_path = os.path.join(myargs.output_dir, prefix)
    if not os.path.isdir(curr_assembly_dir_path):
        os.mkdir(curr_assembly_dir_path)
    else:
        log_file.write("%s\tN/A\tN/A\tN/A\tN/A\tN/A\n" % prefix)
        continue    

    curr_COG_path = os.path.join(myargs.COG_domtblout_dir, "%s.table.domains" % prefix)
    check_input_dir_or_file(curr_COG_path, os.path.isfile, "FATAL ERROR: Domain table '%s' is not found", True) 
    (id_to_COGs, COGs) = udav_soft.read_Pfam_output(curr_COG_path, myargs.e_value_threshold, False, -100, do_not_get_features=True) # FIX: v.1.2 This was done for each nucleotide record

    if myargs.Pfam_domtblout_dir != None:
        curr_Pfam_path = os.path.join(myargs.Pfam_domtblout_dir, "%s.table.domains" % prefix)
        check_input_dir_or_file(curr_Pfam_path, os.path.isfile, "FATAL ERROR: Domain table '%s' is not found", True)
        (id_to_Pfam, Pfam) = udav_soft.read_Pfam_output(curr_Pfam_path, myargs.e_value_threshold, False, -100, do_not_get_features=True)

    curr_nucl_path = os.path.join(myargs.nucl_fasta_dir, "%s.nucl.fasta" % prefix) 
    check_input_dir_or_file(curr_nucl_path, os.path.isfile, "FATAL ERROR: Nucleotide sequence file '%s' is not found", True) 

    curr_fasta_path = os.path.join(myargs.fasta_dir, fasta_filename)
    (proteins, found) = udav_fasta.read_fasta(curr_fasta_path, None, None, dict(), 100500, False, None, None, None, None, "URef")
    print ("    Red %i proteins..." % len(proteins))
    nucl_records_to_proteins = dict() # Dict with nucleotide record accessions as keys and protein list as values
    for p in proteins:
        curr_nucl_record = p.source_record.split(" ")[0]
        if not curr_nucl_record in nucl_records_to_proteins:
            nucl_records_to_proteins[curr_nucl_record] = list()
        nucl_records_to_proteins[curr_nucl_record].append(p)

    (nucl_records, found) = udav_fasta.read_fasta (curr_nucl_path, None, None, dict(), 100500, False, None, None, None, None, "COGNAT_nucl")    
    print ("    Red %i nucleic acid records..." % len(nucl_records))
    protein_number = 0
    empty_nucl_number = 0
    for n in nucl_records:
        if not n.gi in nucl_records_to_proteins: # This nucleotide record does not contain proteins and thus is not required in the database
            empty_nucl_number += 1
            continue
        curr_gbk_dir_path = os.path.join(curr_assembly_dir_path, n.gi)
        if os.path.isdir(curr_gbk_dir_path): #FIX: version 1.1 (analysis can now be stopped)
            continue
        os.mkdir(curr_gbk_dir_path)              
        g_fasta_file = open(os.path.join(curr_gbk_dir_path, "g_%s.fasta" % n.gi), "w") # gbk sequence (in 'COGNAT_nucl' format)
        n.print_fasta(g_fasta_file, 60)
        g_fasta_file.close()

        a_gis_file = open(os.path.join(curr_gbk_dir_path, "a_%s.gis" % n.gi), "w") # list of protein GIs
        a_ids_file = open(os.path.join(curr_gbk_dir_path, "a_%s.ids" % n.gi), "w") # list of protein IDs
        c_txt_file = open(os.path.join(curr_gbk_dir_path, "c_%s.txt" % n.gi), "w") # COG database HMMer results
        d_txt_file = open(os.path.join(curr_gbk_dir_path, "d_%s.txt" % n.gi), "w") # Pfam database HMMer results
        m_txt_file = open(os.path.join(curr_gbk_dir_path, "m_%s.txt" % n.gi), "w") # TMHMM search
        p_fasta_file = open(os.path.join(curr_gbk_dir_path, "p_%s.fasta" % n.gi), "w") # protein sequences (in 'COGNAT' format)
        protein_list = nucl_records_to_proteins[n.gi]
        protein_number += len(protein_list)
        #id_to_domains[protein_id][domain_name] - dict of data
        for p in protein_list:       
            (primary_id, secondary_id) = get_primary_and_secondary_id(p, myargs.mode)                 
            a_gis_file.write("%s\n" % primary_id)
            a_ids_file.write("%s\n" % primary_id)
            full_protein_id = "gi|%s|ref|%s" % (p.gi, p.protein_id)
            #full_protein_id = p.protein_id
            #print_domain_hits_for_protein(full_protein_id, id_to_COGs, p.protein_id, c_txt_file)
            #print_domain_hits_for_protein(full_protein_id, id_to_COGs, p.locus, c_txt_file)
            print_domain_hits_for_protein(full_protein_id, id_to_COGs, primary_id, c_txt_file) 
            if myargs.Pfam_domtblout_dir != None: 
                #print_domain_hits_for_protein(full_protein_id, id_to_Pfam, p.protein_id, d_txt_file)
                print_domain_hits_for_protein(full_protein_id, id_to_Pfam, primary_id, d_txt_file)
            if full_protein_id in tmhmm_features:
               #m_txt_file.write("%s\t%i\n" % (p.protein_id, tmhmm_features[full_protein_id]))
               m_txt_file.write("%s\t%i\n" % (primary_id, tmhmm_features[full_protein_id]))
            else:
               #m_txt_file.write("%s\t%i\n" % (p.protein_id, 0))
               m_txt_file.write("%s\t%i\n" % (primary_id, 0))
            #p_fasta_file.write("%s" % p.get_name_in_format("ID", "COGNAT"))
            #p_fasta_file.write("%s" % p.get_name_in_format("locus", "COGNAT"))
            p_fasta_file.write("%s" % p.get_name_in_format(myargs.mode, "COGNAT"))
            p_fasta_file.write("%s\n\n" % p.sequence)
        a_gis_file.close()
        a_ids_file.close()
        c_txt_file.close()
        d_txt_file.close()
        m_txt_file.close()
        p_fasta_file.close()
    t_end = time.time()
    t_global = (t_end - t_global_start)/60
    log_file.write("%s\t%i\t%i\t%i\t%.1f sec\t%.1f min\n" % (prefix, len(nucl_records), empty_nucl_number, protein_number, t_end - t_start, t_global))
log_file.close()
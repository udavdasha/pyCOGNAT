#!/usr/bin/env python
"""
Module for basis COGNAT functionality. 
Used in the COGNAT.py and pyCOGNAT.py scripts.
"""
curr_version = 1.1

import os
import copy
import udav_soft, udav_fasta

def remove_none(any_list):
    i = 0
    while i < len(any_list):
          if any_list[i] == None:
              any_list.pop(i)
              i -= 1
          i += 1
    return    

def map_proteins_to_gbk(protein_ids, assembly_folder):
    """
    Reads an <assembly_folder> (a folder for a certain genome assembly, e.g. 'GCA_002839445.1_ASM283944v1')
    and returns a dictionary with gbk names as keys and lists of respective 'protein_ids' as values
    """
    gbk_to_ids = dict()

    gbk_folders = os.listdir(assembly_folder)
    for gbk in gbk_folders:
        if len(protein_ids) == 0: # All ids were mapped
            break
        ids_filename = os.path.join(assembly_folder, gbk, "a_%s.ids" % gbk)
        if not os.path.isfile(ids_filename):
            print ("[FATAL ERROR]: File with ids '%s' is missing" % ids_filename)
            raise ReadDatabaseError
        curr_ids = udav_soft.read_plain_list(ids_filename)
        m = 0
        for i in range(len(protein_ids)):
            if protein_ids[i] in curr_ids:
                if not gbk in gbk_to_ids: 
                    gbk_to_ids[gbk] = list()                           
                gbk_to_ids[gbk].append(protein_ids[i])
                m += 1
                protein_ids[i] = None
        remove_none(protein_ids)
        #print ("Total %i protein IDs found in a gbk folder '%s'" % (len(curr_ids), gbk))
        #print (curr_ids)
        #print ("    %i proteins mapped to this gbk (%i remained)" % (m, len(protein_ids)))
        #input()
    return gbk_to_ids

def map_proteins_to_gbff_and_gbk(protein_ids, COGNAT_database_dir, pbar = None, tkwindow = None):
    """
    Reads main COGNAT database folder <COGNAT_database_dir> and returns dictionary with requred 
    <protein_ids> as keys and a tuple of gbff and corresponding gbk as values.
    If no information was found, an id would be linked with None value.
    If <pbar> and <tkwindow> are given, the information about the progress of current operation
    will be written to the <pbar> instead of the console.
    """
    print ("Mapping proteins to gbff and gbk...")
    ids_to_gbff_and_gbk = dict()    
    gbff_folders = os.listdir(COGNAT_database_dir)
    n = 0
    for gbff in gbff_folders:
        n += 1
        if (pbar == None) or (tkwindow == None):
            print ("    Working with %i gbff directory (out of %i)" % (n, len(gbff_folders)))
        else:
            pbar["value"] = int(100 * n / len(gbff_folders))
            tkwindow.update()
        path_to_gbff = os.path.join(COGNAT_database_dir, gbff)
        gbk_to_ids = map_proteins_to_gbk(protein_ids, path_to_gbff)
        for gbk in gbk_to_ids.keys():
            for protein_id in gbk_to_ids[gbk]:
                ids_to_gbff_and_gbk[protein_id] = (gbff, gbk)

    if not (len(protein_ids) == 0): # Some protein ids were not assigned
        print ("[WARNING] Genome data was not found for %i proteins:" % len(protein_ids))
        for p in protein_ids:
            print ("%s" % p)
    else:
        print ("Genome data found for all %i proteins" % len(ids_to_gbff_and_gbk.keys()))
    return ids_to_gbff_and_gbk 

def read_alina_output(alina_output_filename):
    """
    This method reads output of the 'genomic_environment.py' module written by Alina Erina and returns:
    1) A dictionary of dictionaries 'gbff' -> 'gbk' -> list() of <udav_fasta.Annotated_sequence> objects
    corresponding to the prediction data;
    2) A dictionary with the assignment of proteins to COGs
    """
    add_proteins = dict()
    add_id_to_domains = dict()
    alina_output_file = open(alina_output_filename)
    i = 0
    for string in alina_output_file:
        #GCF000196135.1	NC_005090	479612..481042	COG1034	GCF_000196135.1_ASM19613v1/NC_005090
        #  [0]           [1]         [2]          [3]                 [4]          [5]          [6]                 [7]
        #addprot1	COG1034	GCF000196135.1	NC_005090	479612..481042	9.80e-110	342	GCF_000196135.1_ASM19613v1/NC_005090
        string = string.strip()
        if len(string) == 0:
            continue
        if string[0] == "#":
            continue
        i += 1
        fields = string.split("\t")
        prot_id = fields[0]
        #prot_id = "addprot%i" % i
        curr_gbff = fields[7].split("/")[0]
        curr_gbk = fields[7].split("/")[1]
        if not curr_gbff in add_proteins:
            add_proteins[curr_gbff] = dict()
        if not curr_gbk in add_proteins[curr_gbff]:
            add_proteins[curr_gbff][curr_gbk] = list()

        coord1 = int(fields[4].split("..")[0])
        coord2 = int(fields[4].split("..")[1])
        begin = min(coord1, coord2)
        end = max(coord1, coord2)
        direction = None
        if coord2 > coord1: # Direct
            direction = 1
        else:
            direction = -1

        curr_protein = udav_fasta.Annotated_sequence(prot_id, prot_id, prot_id, "Unknown product", "Unknown organism", begin, end, direction,
                                                     "Unknown taxonomy", "%s @ %s" % (curr_gbk, curr_gbff), "Unknown FASTA")
        add_proteins[curr_gbff][curr_gbk].append(curr_protein)

        curr_cog_mark = fields[1]
        evalue = fields[5]
        score = fields[6]
        add_id_to_domains[prot_id] = dict()
        add_id_to_domains[prot_id][curr_cog_mark] = "%s..%s..%s..%s..%s..%s..%s " % (1, end - begin + 1, float(evalue), float(score), "Unk", "Unk", "Unk")
    alina_output_file.close()
    print ("Obtained %i proteins from the additional prediction file '%s'" % (i, alina_output_filename))
    return (add_proteins, add_id_to_domains)

def add_predicted_genes(curr_gene_neighborhoods, predicted_proteins, missing):
    """
    This method takes a list of <curr_gene_neighborhoods> (each of which is a list of protein objects
    of the <udav_fasta.Annotated_sequence> class) and a list of additional <predicted_proteins>
    which should be inserted in the right places inside of them. Protein IDs of these proteins should
    be marked with the ^ signs (e.g. '^YP_00000001.1^') so that they could be properly shown.
    """
    placed_proteins = dict()
    for neighborhood in curr_gene_neighborhoods:
        neighborhood.sort(key = lambda k: k.gene_begin, reverse = False)
        first_begin = neighborhood[0].gene_begin
        first_gene_pid = neighborhood[0].protein_id
        last_end = neighborhood[len(neighborhood) - 1].gene_end
        last_gene_pid = neighborhood[len(neighborhood) - 1].protein_id
        for add_prot in predicted_proteins:
            if add_prot.gene_begin in range(first_begin, last_end): # This protein falls into this neighborhood
                add_prot.protein_id = "^%s^" % add_prot.protein_id
                neighborhood.append(copy.deepcopy(add_prot))
                placed_proteins[add_prot.protein_id] = True
                #for i in range(len(neighborhood) - 1):
                #    if add_prot.gene_begin in range (neighborhood[i].gene_end, neighborhood[i + 1].gene_begin): # This additional protein should be inserted after [i]
                #        add_prot.protein_id = "^%s^" % add_prot.protein_id
                #        neighborhood.insert(i + 1, copy.deepcopy(add_prot))
                #        placed_proteins[add_prot.protein_id] = True
                #        break
            else:
                print ("Protein '%s' (%i..%i) was not assigned to a neighborhood!" % (add_prot.protein_id, add_prot.gene_begin, add_prot.gene_end))
                print ("Begin of the neighborhood: %i (start of gene %s)" % (first_begin, first_gene_pid))
                print ("End of the neighborhood: %i (end of gene %s)" % (last_end, last_gene_pid))
        last_end = neighborhood[len(neighborhood) - 1].gene_end

        neighborhood.sort(key = lambda k: k.gene_begin, reverse = False)

    if len(predicted_proteins) != len(placed_proteins.keys()):
        for p in predicted_proteins:
            if not p.protein_id in placed_proteins:
                missing[p.protein_id] = True

def proceed_gbk_data(assembly_dir, gbk, req_protein_ids, neighbor_number, other_req_protein_ids = dict(), separate = False):
    """
    This method will return gene neighborhoods of genes under given <protein_ids> from the
    respective <gbk> in <assembly_dir>. Genes found in the <req_protein_ids> will have their 'protein_ids'
    marked with asterisks (like *PKL40181.1*)

    Size of the gene neighborhood should not be fixed to the <neighbor_number> to the left and to the right of the FIRST
    occuring gene in <req_protein_ids>, thus is it relative to the last occuring gene.

    Also if a dictionary of <other_req_protein_ids> is given, these proteins are also considered interesting and
    thus the neighborship is extended based on them.

    If <separate> is set to true, neighborhoods will be given for each protein in <req_protein_ids>, they will not be merged if
    they overlap (this should be used for COGNAT-like behavior) 
    """
    neighborhoods = list()
    fasta_filename = os.path.join(assembly_dir, gbk, "p_%s.fasta" % gbk)
    (protein_list, found) = udav_fasta.read_fasta(fasta_filename, None, None, dict(), 100500, False, None, None, None, None, "COGNAT", correct_gene_begin = True)
    protein_list = sorted(protein_list, key = lambda p: p.gene_begin)
    #print ("---- Proceeding gbk-file %s!" % gbk)    
    for i in range(len(protein_list)):
        prot = protein_list[i] 
        if not prot.protein_id in req_protein_ids:
            continue
        #print ("-- Neighborhood for the protein %s:" % prot.protein_id)
        curr_neighborhood = list()
        # 1) <---- Going left
        #print (req_protein_ids)
        #print (other_req_protein_ids)
        last_left_index = max(0, i - neighbor_number)    
        #print ("Going left: to %i" % last_left_index)
        l = i       
        while l >= last_left_index:
            if not separate: #FIX: COGNATE emulation should produce separate surroundings
                if (protein_list[l].protein_id in req_protein_ids) or (protein_list[l].protein_id in other_req_protein_ids): # This arrow will be colored
                    last_left_index = max(0, last_left_index - 1)                 
                    #print ("    last_left_index changed to: %i" % last_left_index)
                if protein_list[l].protein_id in req_protein_ids:                
                    req_protein_ids.remove(protein_list[l].protein_id)                
                    protein_list[l].protein_id = "*%s*" % protein_list[l].protein_id
            curr_protein_object = copy.deepcopy(protein_list[l])
            if (l == i) and separate: # This is target protein, but marking is not done for all proteins in <req_protein_ids> because <separate> is True
                curr_protein_object.protein_id = "*%s*" % curr_protein_object.protein_id
            curr_neighborhood.append(curr_protein_object)
            #print ("i = %i, l = %i, protein: %s" % (i, l, protein_id_to_add))
            l -= 1
        curr_neighborhood.reverse()
        # 2) Going right ---->
        last_right_index = min(len(protein_list) - 1, i + neighbor_number)
        #print ("Going right: to %i" % last_right_index)
        r = i + 1
        while r <= last_right_index:
            if not separate: #FIX: COGNATE emulation should produce separate surroundings
                if (protein_list[r].protein_id in req_protein_ids) or (protein_list[r].protein_id in other_req_protein_ids): # This arrow will be colored                                
                    last_right_index = min(len(protein_list) - 1, last_right_index + 1)                                 
                    #print ("    last_right_index changed to: %i" % last_right_index)
                if protein_list[r].protein_id in req_protein_ids:                
                    req_protein_ids.remove(protein_list[r].protein_id)                
                    protein_list[r].protein_id = "*%s*" % protein_list[r].protein_id
            curr_neighborhood.append(protein_list[r])
            #print ("i = %i, r = %i, protein: %s" % (i, r, protein_list[r].protein_id))
            r += 1
        neighborhoods.append(curr_neighborhood)
        #input()
    return neighborhoods

def make_deep_copy_of_each_element(array_of_objects):
    for i in range(len(array_of_objects)):
        copy_of_element = copy.deepcopy(array_of_objects[i])
        array_of_objects[i] = copy_of_element

def fix_neighborhood_direction(neighborhood, COGNAT_database_dir, gbff, gbk):
    """
    This method will turn all neighborhoods direction as if we look on the opposite DNA chain
    """
    fixed_neighborhood = list()    
    gbk_seq_records = read_nucl_data(COGNAT_database_dir, gbff, gbk)
    gbk_L = 0
    try:
        if gbk_seq_records[0].gene_direction == 0: # Normal COGNAT database
            gbk_L = len(gbk_seq_records[0].sequence)
        else:
            gbk_L = gbk_seq_records[0].gene_direction # Artificial place for the length
    except:
        print ("[WARNING] Cannot read file with gbk sequence: '%s'" % gbk_seq_filename)
    for p in neighborhood:
        fixed_neighborhood.append(p)
    fixed_neighborhood.reverse()
    for p in fixed_neighborhood:
        new_gene_begin = gbk_L - p.gene_end + 1
        new_gene_end = gbk_L - p.gene_begin + 1
        p.gene_begin = new_gene_begin
        p.gene_end = new_gene_end
        p.gene_direction = 0 - p.gene_direction
    return fixed_neighborhood

def get_sorted_neighborhoods(protein_ids, ids_to_gbff_and_gbk, COGNAT_database_dir, neighbor_number, fix_direction, add_proteins, pbar = None, tkwindow = None):
    """
    Arguments are as follows:
    * <protein_ids> is a list of protein_ids as they were presented by a user
    * <ids_to_gbff_and_gbk> is a dictionary created by the <map_proteins_to_gbff_and_gbk()> method
    * <COGNAT_database_dir> is a path to COGNAT database folder
    * <neighbor_number> is a maximal number of protein genes to show from both sides from the target gene
    * <fix_direction> - if True, all directions of neighborhoods will be the same
    * <add_proteins> is a dictionary following the same logic as the <gbff_to_gbk_and_ids> dictionary.
    Returns an ordered list of lists, each element of the latter is an <udav_fasta.Annotated_sequence>.
    IDs of proteins which are target will be flanked with the '*' sign while IDs of proteins which are
    added from the <add_proteins> are flanked with the '^' sign.

    If <pbar> and <tkwindow> are given, the information about the progress of current operation
    will be written to the <pbar> instead of the console.
    """
    print("Obtaining neighborhoods for proteins of interest...")
    missing = dict()
    # 1) Reversing <ids_to_gbff_and_gbk> dict
    gbff_to_gbk_and_ids = dict()   
    for protein_id in ids_to_gbff_and_gbk.keys():
        gbff = ids_to_gbff_and_gbk[protein_id][0]
        gbk = ids_to_gbff_and_gbk[protein_id][1]
        if not gbff in gbff_to_gbk_and_ids:
            gbff_to_gbk_and_ids[gbff] = dict()
        if not gbk in gbff_to_gbk_and_ids[gbff]:
            gbff_to_gbk_and_ids[gbff][gbk] = list()
        gbff_to_gbk_and_ids[gbff][gbk].append(protein_id)
    # 2) Getting neighborhoods
    protein_ids_in_neighborhoods = dict() # All proteins occuring in neighborhoods FIX: version 1.1 (values are True for proteins with reversed direction and False for other)
    target_id_to_neighborhood = dict() # Dictionary of protein ids which were targets (their IDs are marked) to respective neighborhood
    n = 0
    for gbff in gbff_to_gbk_and_ids.keys():
        n += 1
        if pbar != None:
            pbar["value"] = int(100 * n / len(gbff_to_gbk_and_ids.keys()))
            tkwindow.update()
        else:
            print ("    Working with %i gbff directory (out of %i)" % (n, len(gbff_to_gbk_and_ids.keys())))
        path_to_gbff = os.path.join(COGNAT_database_dir, gbff)
        for gbk in gbff_to_gbk_and_ids[gbff]:
            curr_gene_neighborhoods = proceed_gbk_data(path_to_gbff, gbk, gbff_to_gbk_and_ids[gbff][gbk], neighbor_number, dict(), separate = True) #FIX: option <separate> added not to merge neighborhoods
            make_deep_copy_of_each_element(curr_gene_neighborhoods) #FIX: if neighborhoods were overlapping, erroneous fixing occured
            if gbff in add_proteins: # ALINA'S PSEUDO PROTEINS SHOULD BE ADDED TO NEIGHBORHOODS
                if gbk in add_proteins[gbff]:
                    add_predicted_genes(curr_gene_neighborhoods, add_proteins[gbff][gbk], missing)
            for single_neighborhood in curr_gene_neighborhoods:
                #print ("New neighborhood: %s" % single_neighborhood)
                was_reversed = False #FIX: v. 1.1
                for p in single_neighborhood: 
                    real_protein_id = p.protein_id.strip("*")
                    if real_protein_id != p.protein_id: # This is the target gene
                        #print ("----Target gene found: '%s'" % real_protein_id)
                        if fix_direction and (p.gene_direction == -1): # Reverse direction, correction required
                            single_neighborhood = fix_neighborhood_direction(single_neighborhood, COGNAT_database_dir, gbff, gbk)
                            was_reversed = True
                            #print ("----Changed direction: '%s'" % single_neighborhood)
                        target_id_to_neighborhood[real_protein_id] = single_neighborhood
                for p in single_neighborhood: 
                    protein_ids_in_neighborhoods[p.protein_id.strip("*")] = was_reversed
    if len(missing.keys()) != 0:
        print ("WARNING: not all proteins were added to a genetic neighborhood, the following %i are missing:" % len(missing.keys()))
        ids_to_print = ""
        for p in missing.keys():
            ids_to_print += "%s, " % p
        print(ids_to_print.strip(" ,"))

    # 3) Sorting neighborhoods
    neighborhoods = list() 
    for protein_id in protein_ids:
        if protein_id in target_id_to_neighborhood:
            neighborhoods.append(target_id_to_neighborhood[protein_id])
        else:
            print("[WARNING] Genomic neighborhood was not found for protein '%s'" % protein_id)
            pseudo_protein = udav_fasta.Annotated_sequence(protein_id, protein_id, protein_id, "Unk", "Unk", 0, 1000, 1, "Unk", "Unk", "Unk")
            neighborhoods.append([pseudo_protein])
    return (neighborhoods, protein_ids_in_neighborhoods)

def read_nucl_data(COGNAT_database_dir, gbff, gbk):
    """
    Method will read all available nucleotide data for given parameters
    """
    path_to_nucl_seq = os.path.join(COGNAT_database_dir, gbff, gbk, "g_%s.fasta" % gbk)
    (nucl_records, found) = udav_fasta.read_fasta(path_to_nucl_seq, None, None, dict(), 100500, False, None, None, None, None, "COGNAT_nucl")
    return nucl_records

def assign_proteins_to_domains(id_to_domains, COGNAT_database_dir, gbff, gbk, req_protein_ids = None, max_e_value = 0.1, domain_type = "COG"):
    """
    Method will assign proteins to COGs based on the file from the COGNAT database 'c_<gbk>.txt'
    (or respective Pfam file 'd_<gbk>.txt', if <domain_type> is 'Pfam') and fill given <id_to_domains> 
    parameter. 
    If a dictionary of <req_protein_ids> is not None, only proteins from it are taken.
    Mimics functionality as the <udav_soft.read_Pfam_output> method (except for <unite_same> etc)
    Hits with i-evalue higher then given <max_e_value> are not shown.
    ! Check which e-value threshold was used for the COGNAT database creation. It should not be
      lower than the one given here
    """
    base_filename = "c_%s.txt" % gbk
    if domain_type == "Pfam":
        base_filename = "d_%s.txt" % gbk
    domain_hits_filename = os.path.join(COGNAT_database_dir, gbff, gbk, base_filename)
    domain_hits_file = open(domain_hits_filename)
    for string in domain_hits_file:
        #protein domain evalue  score profile % begin   end
        #MK0001	COG0430	6e-140	464.5	99	2	349
        string = string.strip()
        if len(string) == 0:
            continue
        fields = string.split("\t")
        protein_id = fields[0]
        if req_protein_ids != None:
            if not protein_id in req_protein_ids:
                continue 
        domain_name = fields[1]
        evalue = fields[2]
        score = fields[3]
        hmm_covered = fields[4]
        begin = fields[5]
        end = fields[6]
        #---- 1) E-value check
        if float(evalue) > max_e_value:
            continue
        #---- 2) Adding to hash with protein id as keys
        if not protein_id in id_to_domains:
            id_to_domains[protein_id] = dict()
            id_to_domains[protein_id][domain_name] = ""
        else:
            if not domain_name in id_to_domains[protein_id]:
                id_to_domains[protein_id][domain_name] = ""
        new_data = "%s..%s..%s..%s..%s..%s..%s " % (begin, end, evalue, score, "Unk", "Unk", hmm_covered) 
        if not new_data in id_to_domains[protein_id][domain_name]: #FIX: if data about this domain was already added in previous neighborhood analysis
            id_to_domains[protein_id][domain_name] += new_data
    domain_hits_file.close()

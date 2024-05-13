"""
Module for the very base classes used in all scripts
------- Version: 1.7.0

Methods included in this module:
        1) dict read_feature_file (feature_file, feature_type)
           Reads <feature_file>, fills dictionary <feature_type> with found feature names
           and returns a dictionary with a link between sequence id and the list of features

        2) void print_pure_sequences (seqs, output_filename, id_only, remove_gaps)
           Prints sequences of <seqs> list purified from gaps to the output_filename.
           If <id_only> parameter is True, name of each sequence will be shorteded to its ID.
           If <remove_gaps> is True, gaps will be removed

        3) list read_alignment (input_filename)
           Reads alignment into a list of <Sequence> objects from file <input_filename>

        4) void correspond_back(seqs, correspond_filename)

        5) list proceed_params (parameters)
      
Classes included in this module:
	1) Sequence (-> Alignment_sequence, -> Annotated_sequence, -> Featured_sequence)
           ~ Variables: ~
           str  name          - string with initial description line
           str  ID            - part of the name before 1st space (as considered by fasta format)
           str  sequence      

           ~ Methods: ~
           void print_fasta(fasta_file)           - print fasta format sequence into given file        
           bool length_in_range(min_len, max_len) - checks if sequence length is into the given range

	2) Feature
           ~ Variables: ~
           str  name          - name of the feature (e.g., 'TMHMM', 'COG0001' or 'HELIX')
           list regions       - list of regions (like '1..5')

           ~ Methods: ~
           int get_begin(region)   - returns begin value of given region (e.g. for '1..5' returns 1)
           int get_end(region)     - returns end value of given region (e.g. for '1..5' returns 5)

	3) Featured_sequence (<- Sequence)
           ~ Variables: ~
        -> str  name          - string with initial description line
        -> str  ID            - part of the name before 1st space (as considered by fasta format)
        -> str  sequence
           list features      - list of objects of the <Feature> class

           ~ Methods: ~
        -> void print_fasta(fasta_file)           - print fasta format sequence into given file        
        -> bool length_in_range(min_len, max_len) - checks if sequence length is into the given range

	4) Alignment_sequence (<- Sequence, -> Seq_vertex)
           ~ Variables: ~
	-> str name
        -> str ID
        -> str sequence

           ~ Methods: ~
           void remove_limits (long_names)                      - removes nasty symbols added by JalView
           str correct_organism (correct_orgs, expanded_orgs): - corrects organism name or expands it
"""
import os, sys
import re

def print_hash(the_hash, hash_name):
    print ("Printing hash with name '%s'" % hash_name)
    for key in the_hash.keys():
         print ("'%s'\t'%s'" % (key, the_hash[key]))

def bubble_sort_keys(dictionary, direct, values = True): 
    keys = dictionary.keys()
    i = 0
    while i < len(keys):
        j = i
        while j < len(keys):
            change = False
            if values == True:  # Sort by value itself (considered int or float)
                if (direct == True) and (dictionary[keys[i]] < dictionary[keys[j]]):
                    change = True
                if (direct == False) and (dictionary[keys[i]] > dictionary[keys[j]]):
                    change = True
            if values == False:  # Sort by length of value (considered a list or tuple)
                if (direct == True) and ( len(dictionary[keys[i]]) < len(dictionary[keys[j]]) ):
                    change = True
                if (direct == False) and ( len(dictionary[keys[i]]) > len(dictionary[keys[j]]) ):
                    change = True

            if change:
               i_key = keys[i]
               keys[i] = keys[j]
               keys[j] = i_key
            j += 1
        i += 1
    return keys
              
def print_pure_sequences(seqs, output_filename, id_only, remove_gaps):
    output_file = open(output_filename, "w")
    for s in seqs:
        name = s.name
        sequence = s.sequence
        if remove_gaps == True:
            sequence = sequence.replace("-", "")
        if id_only == True:
            name = s.ID
        output_file.write(">%s\n%s\n\n" %(name, sequence))
    output_file.close()    

    if id_only == True:
        output_correspond = open(output_filename + ".correspond", "w")
        for s in seqs:
            output_correspond.write("%s\t%s\n" % (s.ID, s.name))
        output_correspond.close()

def correspond_back(seqs, correspond_filename):
    corr = dict()
    correspond_file = open(correspond_filename, "r")
    for string in correspond_file:
        string = string.strip()
        if len(string) != 0:
            fields = string.split("\t")
            corr[fields[0]] = fields[1]
    correspond_file.close()
    
    n = 0
    for i in range(len(seqs)):
        if seqs[i].ID in corr:
            seqs[i].name = corr[seqs[i].ID]
            n += 1
        else:
            print ("WARNING: sequence with ID %s was not found in table %s!" % (seqs[i].ID, correspond_filename))
    print ("Total %s sequence names fixed (from %i expected)" % (n, len(seqs)))

def proceed_params (parameters):
    parameters[0] = os.path.abspath(parameters[0]) # Working directory
    i = 1
    while i < len(parameters):
        if parameters[i] != None:
            parameters[i] = os.path.join(parameters[0], parameters[i]) 
        i += 1
    return parameters

class Sequence:
    def __init__(self, name, sequence, define_proper_id = True):
        self.name = name
        self.ID = name.split(" ", 1)[0]
        self.sequence = sequence
        if (self.ID == self.name) and define_proper_id:
            self.define_id()

    def print_fasta (self, fasta_file, size = None, part_begin = None, part_end = None):
        fasta_file.write(">" + self.name + "\n")
        sequence = self.sequence
        if (part_begin != None) and (part_end != None):
            sequence = self.sequence[part_begin - 1:part_end]
        if size == None:
            fasta_file.write(sequence + "\n\n")
        else:
            start = 0
            while start < len(sequence):
                end = start + size
                fasta_file.write(sequence[start:end] + "\n")
                start += size
            fasta_file.write("\n")

    def length_in_range(self, min_len, max_len):
        result = True
        if (min_len != None) and (max_len != None):
            if (len(self.sequence) < min_len) or (len(self.sequence) > max_len):
                result = False
        return result

    def define_id(self):
        self.ID = self.name.split("|", 1)[0]
        #id_parts = self.ID.split("_") # In case ID contains cut data (e.g YP_000001_5-45)
        #if len(id_parts) > 2:
        #    self.ID = id_parts[0] + id_parts[1]
        self.organism = self.name.split("|")[-1]

class Feature:
    def __init__(self, single_feature, not_exact = False):
        #[TMHMM] 1..5,9..11,17..25
        self.name = single_feature.split(" ")[0].strip("[]")
        if not_exact:
            self.name = self.name.split("/", 1)[0]
        self.regions = single_feature.split(" ")[1].split(",")

    def get_begin(self, region):
        fields = region.split("..")
        if fields[0].isdigit():
            return int(fields[0])
        else:
            print ("WARNING: feature region '%s'" % region)
            return -1
            
    def get_end(self, region):
        fields = region.split("..")
        if len(fields) == 1:
            return int(fields[0])
        else:
            return int(fields[1])

    def get_length(self):
        length = 0
        for r in self.regions:
            fields = r.split("..")
            length += int(fields[1]) - int(fields[0]) + 1
        return length
                                
class Featured_sequence(Sequence):
    def __init__(self, name, sequence, features):
        Sequence.__init__(self, name, sequence)
        self.features = features

    def get_features_order(self): #Returns string with all features except TMHMM placed by their start position
        feature_order = list()
        for f in self.features:
            if f.name == "TMHMM":
                continue
            for r in f.regions:
                curr_start = f.get_begin(r)
                feature_order.append((f.name, curr_start))
        i = 0
        while i < len(feature_order):
            j = i
            while j < len(feature_order):
                if feature_order[j][1] < feature_order[i][1]:
                   i_tuple = feature_order[i]
                   feature_order[i] = feature_order[j]
                   feature_order[j] = i_tuple
                j += 1
            i += 1
        feature_string = ""
        for f in feature_order:
            feature_string += f[0] + " "
        feature_string = feature_string.strip()
        return feature_string
                      
    def range_features_length(self): #All features except for the first (TMHMM) will be ranged by length
        i = 1
        while i < len(self.features):
            j = i
            while j < len(self.features):
                if self.features[i].get_length() > self.features[j].get_length():
                    i_feature = self.features[i]
                    self.features[i] = self.features[j]
                    self.features[j] = i_feature   
                j += 1
            i +=1

    def correct_features(self, gap_symbol): 
        for i in range(len(self.sequence)):
            curr_letter = self.sequence[i]
            if curr_letter == gap_symbol:
                for f in range(len(self.features)):
                    for r in range(len(self.features[f].regions)):
                        curr_region = self.features[f].regions[r]
                        start = self.features[f].get_begin(curr_region)
                        end = self.features[f].get_end(curr_region)
                        if start >= i + 1:
                            start = start + 1
                        if end >= i + 1:
                            end = end + 1
                        self.features[f].regions[r] = "%i..%i" % (start, end)
        #print ("Correcting features in sequence %s DONE" % self.ID)

    def get_feature_string(self):
         feature_string = "%s\t" % self.ID
         for f in self.features:
             feature_string += "[%s] " % f.name
             for r in f.regions:
                 start = f.get_begin(r)
                 if start == -1:
                     print ("FATAL ERROR in protein %s" % self.ID)
                     sys.exit()
                 end = f.get_end(r)
                 feature_string += "%s..%s," % (start, end)
             feature_string = feature_string.strip(",")
             feature_string += "\t"
         feature_string = feature_string.strip("\t")
         return feature_string

    def vary_coordinate(self, coordinate, coord_type, fixed, range_to_vary, gap_symbol):
        original_gap_num = 0
        inc = 0
        if coord_type == "START":
            inc = 1
            original_gap_num = self.sequence[coordinate - 1 : fixed].count(gap_symbol)
        if coord_type == "END":
            inc = -1
            original_gap_num = self.sequence[fixed : coordinate].count(gap_symbol)
        best_gap_num = original_gap_num
        best = coordinate

        for i in range_to_vary: 
            new = coordinate + i            
            seq_part = None
            if coord_type == "START":
                seq_part = self.sequence[new - 1 : fixed] # Part for new region (for start)
            if coord_type == "END":
                seq_part = self.sequence[fixed : new] # Part for new region (for end) 
                seq_part = seq_part[::-1]
            
            for s in seq_part:
                if s == gap_symbol:
                    new += inc
                else:
                    break

            if coord_type == "START":
                seq_part = self.sequence[new - 1 : fixed] # Part for new region (for start)
            if coord_type == "END":
                seq_part = self.sequence[fixed : new] # Part for new region (for end) 

            new_gap_num = seq_part.count(gap_symbol)
            if new_gap_num * i * 2 < best_gap_num: # Moving by 1 letter = at least 2 gaps out
                if new_gap_num < best_gap_num:
                    best_gap_num = new_gap_num
                    best = new
        return best

    def correct_TM(self, gap_symbol, vary_value): # Transmembrane feature should always be the first one!
        if len(self.features) == 0:
            return None
        if self.features[0].name == "TMHMM":
            n = 0
            for r in range(len(self.features[0].regions)):
                curr_region = self.features[0].regions[r]
                start = self.features[0].get_begin(curr_region)
                end = self.features[0].get_end(curr_region)
                
                best_start = self.vary_coordinate(start, "START", end, range(0, vary_value + 1), gap_symbol)
                best_end = self.vary_coordinate(end, "END", start, range(0 - vary_value,0), gap_symbol)

                if (start < best_start) or (end > best_end):
                    n += 1
                    self.features[0].regions[r] = "%i..%i" % (best_start, best_end)
            #print ("Total %i cases fixed for sequence %s" % (n, self.ID))

class Alignment_sequence(Sequence):
    def __init__(self, name, sequence, define_proper_id = True):
        Sequence.__init__(self, name, sequence, define_proper_id)
                
    def remove_limits(self, long_names, replace_chars = False):
        if replace_chars: # Fix: version 1.2
            self.name = self.name.replace("/", ".")
            self.name = self.name.replace("|", "_")
            self.ID = self.ID.replace("/", ".")
            self.ID = self.ID.replace("|", "_")
        original_name = self.name
        fields = self.name.split("/", 1)
        req_part = fields[0]
        if long_names == True:
            fields = self.name.split(" ", 1)  
            first_part = fields[0].split("/")
            id_part = ""
            for p in range(len(first_part) - 1):
                id_part += first_part[p] + "/"

            id_part = id_part.strip("/")
            req_part = id_part

            if len(fields) != 1: # This name was really long
                desc_part = fields[1]
                req_part = id_part + " " + desc_part
                print (req_part)
        self.name = req_part
        if self.ID == original_name: #Case when ID is identical to name
            self.define_id()

    def correct_organism (self, correct_orgs, expanded_orgs):
        new_name = self.name
        if expanded_orgs != None:
            fields = name.split("|")
            if not len(fields) < 3:        
                AC = fields[1]
                short_name = fields[2].split("_")[1]
                if short_name in expanded_orgs:
                    new_name = ">" + AC + "|" + expanded_orgs[short_name]
        if correct_orgs != None:
            for k in correct_orgs.keys():
                new_name = name.replace(k, correct_orgs[k])
        return new_name

    def prepare_organism (self):
        new_name = self.name
        bad_symbols = "[\[\]\'\,\=]" # FIX: version 1.3 (symbol = is added as bad)
        if re.search(bad_symbols, self.name) != None: # Uneligble characters found
            new_name = re.sub(bad_symbols, "", self.name)
            print ("Replacing: OLD = %s" % self.name)
            print ("           NEW = %s" % new_name)
            
            self.name = new_name

def read_feature_file(feature_filename, feature_type, color_scheme, scheme_only, not_exact = False):
    #YP_00001	[TMHMM] 1..5,9..11,17..25	[PF000001] 5..67
    id_to_features = dict()
    feature_file = open(feature_filename, "r")
    for string in feature_file:
        string = string.strip()
        if len(string) == 0:
            continue
        if string[0] == "#":
            continue
        fields = string.split("\t")
        curr_id = fields[0]
        curr_features = list()
        for i in range(1, len(fields)):
            curr_features.append(Feature(fields[i], not_exact))
            if (not curr_features[-1].name in color_scheme) and scheme_only: # This feature is ignored
                curr_features.pop()
            else:
                if not curr_features[-1].name in feature_type:
                    feature_type[curr_features[-1].name] = True
        id_to_features[curr_id] = curr_features
    feature_file.close()
    return id_to_features

def read_alignment(input_filename, define_proper_id = True):
    seq_list = list()

    input_file = open (input_filename, "r")      #------------------------ Reading alignment
    for string in input_file:
        string = string.strip()
        if len(string) > 0:
            if string[0] == ">":
                string = string.strip(">")
                seq_list.append(Alignment_sequence(string, "", define_proper_id))
            else:
                seq_list[-1].sequence += string
    input_file.close()

    alignment_length = len(seq_list[0].sequence) #------------------------ Checking alignment
    for p in seq_list:
        if len(p.sequence) != alignment_length:
            error = "FATAL ERROR at sequence %s: wrong length %i, expected: %i" % (p.ID, len(p.sequence), alignment_length)
            print (error)
            return error
    return seq_list

def get_featured (sequence_filename, correspond_filename, feature_filename, feature_type, color_scheme, scheme_only):
    alignment = read_alignment(sequence_filename)
    vary_value = 5
    if correspond_filename != None:
        correspond_back(alignment, correspond_filename)
    sequence_features = read_feature_file(feature_filename, feature_type, color_scheme, scheme_only)
    featured_alignment = list()
    n = 0
    print ("Obtaining features for %i sequences..." % len(alignment))
    #debug = open("_AAAAA.txt", "w")
    for s in alignment:
        s.remove_limits(False)
        n += 1
        if s.ID in sequence_features:
            curr_seq = Featured_sequence(s.name, s.sequence, sequence_features[s.ID])
            before = curr_seq.get_feature_string()
            curr_seq.correct_features("-")
            after_correct = curr_seq.get_feature_string()
            curr_seq.correct_TM("-", vary_value)
            after_correct_TM = curr_seq.get_feature_string()
            curr_seq.range_features_length()
            #debug.write("BEFORE        : %s\n" % before)
            #debug.write("AFTER CORRECT : %s\n" % after_correct)
            #debug.write("AFTER TM (%i) : %s\n" % (vary_value, after_correct_TM))
        else:
            print ("WARNING: features not found for the protein '%s'!" % s.ID)
            curr_seq = Featured_sequence(s.name, s.sequence, list())
        featured_alignment.append(curr_seq)
        #print ("Features obtained for %i / %i sequences" % (n, len(alignment)))
    #debug.close()
    print ("DONE obtaining features!")
    return featured_alignment

def sort_by_features(featured_alignment):
    print ("Alignment is being sorted by its features...")
    sorted_alignment = list()
    feature_types = dict() # Here strings with features in them will be saved
    for s in featured_alignment:
        curr_order = s.get_features_order()
        if not curr_order in feature_types:
            feature_types[curr_order] = list()
        feature_types[curr_order].append(s)
    
    keys = bubble_sort_keys(feature_types, False, False) # Ascending sort by length
    for key in keys:
        sorted_alignment.extend(feature_types[key])
    print ("\tDone!")
                    
    return sorted_alignment    

def read_gi_to_tax(filename, strings_directly = False):
    """
    Method reads assignment of GIs to two (or one) taxonomy units and returns a dictionary
    of GI -> last taxon. This file can be produced for example by the <get_gi_to_taxa.py>.
    If <strings_directly> is True, not a filename but strings
    """
    gi_to_tax = dict()
    input_file = filename
    if not strings_directly:
        input_file = open(filename)
    for string in input_file:
        string = string.strip()
        if len(string) == 0:
            continue
        fields = string.split("\t")
        gi = fields[0]
        taxon = ""        
        if len(fields) == 1:
            print ("WARNING: no color is found for %s!" % gi)
            taxon = "Unknown"
        if len(fields) == 2:
            taxon = fields[1]
        if len(fields) > 2:
            taxon = fields[2]
        gi_to_tax[gi] = taxon
    if not strings_directly:
        input_file.close()
    return gi_to_tax

def read_unique_assignment(filename, what = "unique assignment"):
    print ("    Reading %s file..." % what)  
    key_to_value = dict()
    w = 0
    assign_file = open(filename, "r")
    for string in assign_file:
        string = string.strip()
        if len(string) == 0:
            continue
        if string[0] == "#":
            continue
        fields = string.split("\t", 1)
        if fields[0] in key_to_value:
            print ("    [WARNING]: %i is assigned to two different values: %s and %s!" % (fields[0], key_to_value[fields[0]], fields[1]))
            w += 1
        try:
            key_to_value[fields[0]] = fields[1]      
        except IndexError:
            print ("    [FATAL ERROR]: String '%s' does not contain tabulation marks!" % string)
    assign_file.close()
    print ("    Obtained %i key-value pairs; total %i warnings obtained" % (len(key_to_value.keys()), w))
    return key_to_value

def report_info(info, log_file = None):
    """
    Simple method will print given <info> string to STDOUT and also
    will write it with line-break to a given <log_file>
    """
    print (info)
    if log_file != None:
        log_file.write("%s\n" % info)

"""
Module for reading GenBank files 
------- Version: 2.24
        2.0  * My_Ref format replaced for URef
        2.1  * <protein_id> for database construction is now truncated to normal
               RefSeq id. However if this is not possible - "bad" version is retained
        2.2  * Genome_file now have boolean feature "circular"
        2.3  * Genome_file now reads its sequence (following ORIGIN string)
        2.4  * Qualifier 'ribosomal_slippage' is now supported
        2.5  * Qualifier 'trans_splicing' is now supported; GI is extracted as a separate feature
        2.6  * Proper class & function documentation added
        2.7  * Fixing of protein_id is now general
        2.8  * <get_sequence> method added to the <Genome_file> class
        2.9  * <header_strings> attribute added to the <Genome_file> class
             * <print_gbk()> method added to the <Genome_file> class
        2.10 * <Gene.fix_protein_id()> method is fixed to work properly with GenBank too
        2.11 * ACCESSION is now red properly for multiple records
        2.12 * GenPept molecule type is now red properly; DBSOURCE is obtained again
        2.14 * Intergene distances are calculated
        2.15 * <Genome_file.print_nucleotide> method implemented
        2.16 * <merge_until_size> method now works with long strings
        2.17 * Assembly data is added to the <print_proteins> method of the <Genome_file> class
        2.18 * Method <get_gbk_sequence> fixed as it was not printing the last line of the sequence
        2.19 * <counted_features> is added to differ from <req_features>, because of an error with
               the gene number (self.gene_number counted 'CDS', 'tRNA', 'rRNA' IN ADDITION to 'gene'!)
        2.20 * Method <get_old_locus_data> added
        2.21 * Argument 'huge_genome' agged to the <Genome_file.obtain_data()> method
        2.23 * Now 'huge_genome' allows to create a fasta record
        2.24 * <list_under_source> option added to the <split_gbk> method
"""
import os, sys, re
import udav_bio

def get_global_string(key, data, file_size = 79):
    """
     <-------21---------->
     <-5-><----key--->
    '     misc_feature    complement(230015..230611)'
    """
    spaced_size = 21
    num_prefix_spaces = 5
    num_suffix_spaces = spaced_size - num_prefix_spaces - len(key)
    spaced_key = "%s%s%s" % (num_prefix_spaces * " ", key, num_suffix_spaces * " ")
    data_length = file_size - spaced_size
    data_string = ""
    if len(data) < data_length:
        data_string = data
    else:
        print ("WARNING: coordinates are out of range! Current data:")
        print (data)
        data_string = data
    result = "%s%s" % (spaced_key, data_string)
    return result
          
def merge_until_size(fields, separator, max_length):
    """
    Method merges <fields> elements with the <separator> until result length
    exceeds <max_length>.

    Returns: (result_string, remaining_fields)
    """
    result_string = ""
    i = 0
    while i in range(len(fields)):
        new_suffix = "%s%s" % (separator, fields[i])
        if (len(result_string) + len(new_suffix)) > max_length:
            break
        else:
            result_string += new_suffix
        i += 1
    result_string = "%s\n" % result_string.strip(separator)

    if len(result_string) == 1: #FIX: v. 2.16 (too long strings are considered separately)
        if len(fields[0]) < max_length:
            print ("[WARNING] Error in <merge_until_size> method!")
            print ("          The following field is smaller then given <max_length> = %i" % max_length)
            print ("          '%s'" % fields[0])
        part1 = fields[0][:max_length + 1]
        part2 = fields[0][max_length + 1:]
        result_string = "%s\n" % part1
        fields[0] = part2
    return (result_string, fields[i:])

def get_local_string(qualifier, data, data_type, sep_symbol = " ", file_size = 79):
    """
     <-------21---------->        
    '                     /note="COG1812 Archaeal S-adenosylmethionine synthetase"'
    '                     /codon_start=1                                          '
    """
    spaced_size = 21
    data_string = "%s/%s=" % (spaced_size * " ", qualifier)
    add_symbols = 2
    if data_type != "str":
        data_string += "%s" % data
    else:
        add_symbols += 2                
        if (spaced_size + add_symbols + len(qualifier) + len(data)) < file_size:
            data_string += '"%s"' % data
        else:
            data = data.replace("/", " / ")
            fst_data_size = file_size - spaced_size - add_symbols - len(qualifier) + 1
            if qualifier == "translation":
                data_string += '"%s\n' % data[0:fst_data_size]
                new_string = ""
                for char in data[fst_data_size:]:
                     new_string += char
                     if len(new_string) >= (file_size - spaced_size):
                         data_string += "%s%s\n" % (spaced_size * " ", new_string)
                         new_string = ""
                data_string += "%s%s\n" % (spaced_size * " ", new_string)
            else:
                data_parts = data.split(sep_symbol)            
                (fst_string, data_parts) = merge_until_size(data_parts, sep_symbol, fst_data_size)
                data_string += '"%s' % fst_string
                i = 0
                while len(data_parts) > 0:
                    (new_string, data_parts) = merge_until_size(data_parts, sep_symbol, file_size - spaced_size)
                    data_string += "%s%s" % (spaced_size * " ", new_string)
                    i += 1
                    if i > 100:
                        print ("Cannot work out this string (i = %i): '%s'" % (i, data))
                        print (new_string)
                        print (data_parts)
                        sys.exit()
            data_string = data_string.strip("\n")
            data_string += '"'  
    return data_string

def get_gbk_sequence(sequence, column_width = 10, column_num = 6):
    begin_num = 1
    i = 0
    n = 0
    result = list()
    curr_string = ""
    while i < len(sequence):
        curr_part = sequence[i:i+column_width]
        if n == 0: # First block in a string
            begin_num = i + 1
            prefix_spaces = 10 - len("%s" % begin_num) - 1
            curr_string = "%s%s" % (prefix_spaces * " ", begin_num)                        
        curr_string += " %s" % curr_part
        n += 1
        if n >= column_num:
            n = 0
            result.append(curr_string)           
        i += column_width
    if n != 0:
        result.append(curr_string) #FIX: 2.18 (last string was omitted)
    return "\n".join(result)

class Gene:
    """
    Class describes a single gene (either RNA or protein)
    Attributes:
    str  coordinates   - as given by the genome file (with all the regions, if more than 1)
    int  begin         - the first coordinate in the genome (real begin)
    int  end           - the last coordinate in the genome (real end)                                          
    int  direction     - 1 or -1
    str  gene_type     - CDS, tRNA, rRNA etc           
    dict feature       - hash of all gene properties. Hash keys are qualifiers of 
                         the features in gbk file (e.g., 'db_xref', 'product' etc)
    """
    def __init__(self, coordinates, gene_type, begin = 0, end = 0, direction = 0):
        self.coordinates = coordinates
        self.begin = begin
        self.end = end
        self.direction = direction    
        self.gene_type = gene_type
        self.feature = dict()

    def try_get_feature(self, string, last_qualifier):
        """
        Tries to read one of the features from given string.
        Reads all possible features and writes them to the dictionary <self.feature>
        If the string contains a qualifier, returns it. Otherwise return argument qualifier
        """
        qualifier = last_qualifier 

        if string[0] == "/" and string.count("=") > 0: # This string contains a new qualifier
            fields = string.split("=")
            qualifier = fields[0][1:]
            value = fields[1].strip('"')
            if qualifier in self.feature: # This qualifier already was observed (i.e. db_xref)
                # Bug fixed: symbol | was used as a separator in fasta format and thus was replaced for <> 
                self.feature[qualifier] += "<>" + value
            else:                         # This is a new qualifier
                self.feature[qualifier] = value
        else:
            if string == "/pseudo": # This is pseudogene qualifier
                qualifier = "pseudo"
                self.feature[qualifier] = "Yes"
            elif string == "/ribosomal_slippage": # FIX: version 2.4
                qualifier = "ribosomal_slippage"
                self.feature[qualifier] = "Yes"
            elif string == "/trans_splicing": # FIX: version 2.5
                qualifier = "trans_splicing"
                self.feature[qualifier] = "Yes"
            else:                   # This should be appended
                if last_qualifier in self.feature:
                    if last_qualifier == "translation":  # Protein sequences does not require spaces
                        self.feature[last_qualifier] += string.strip('"')
                    else:
                        self.feature[last_qualifier] += " " + string.strip('"')
                        self.feature[last_qualifier] = self.feature[last_qualifier].replace("  ", " ")

        return qualifier

    def extract_data(self, fields):
        """
        Returns a string with required data of the gene
        <fields> -- list of elements with the data qualifiers ("protein_id" etc)
        """
        self.extract_GI()
        values = []
        result = "%s\t%i\t%i\t%i" % (self.gene_type, self.begin, self.end, self.direction)
        for i in range(len(fields)):
            values.append("")
            if fields[i] in self.feature:
               values[i] = self.feature[fields[i]]
            result += "\t" + values[i]
        result += "\n"
        return result

    def extract_GI(self):
        """
        If <self.feature> contains 'db_xref' feature, this method tries to find its GI part
        and upon success the new 'gi' feature is created. Works only for the genes which
        do not already contain 'gi' feature
        """
        if ("gi" in self.feature) or (not "db_xref" in self.feature): # Gi was extracted before or cannot be extracted
            return
       
        all_ref = self.feature["db_xref"].split("<>")
        for r in all_ref:
            key = r.split(":", 1)[0]
            value = r.split(":", 1)[1]
            if key == "GI":
                self.feature["gi"] = value;

    def fix_protein_id(self):
        if "protein_id" in self.feature:
             # FIX: Avoiding things like "REF_PseudoCAP:PA2278<>NP_250968.1"
             # FIX 2.10: This implementation fails for GenBank IDs: BAU40205.1 is converted into AU40205.1!
             #proper_id = re.search("\w{2}_*\d+\.\d+", self.feature["protein_id"])
             proper_id = re.search("\w+\.\d+", self.feature["protein_id"])
             if proper_id != None:
                 self.feature["protein_id"] = proper_id.group(0)
             else:
                 print ("WARNING: Failed to found proper id, remaining id is '%s'" % self.feature["protein_id"])

    def get_fasta(self, record = "Unk_record_name", org = "Unk_org", taxonomy = "Unk_tax", record_id = "Unk_record"):
        """
        Returns string of fasta-formated information. 
        <record>   (OPTIONAL) -- name of the record to which this gene belongs
        <org>      (OPTIONAL) -- name of the organism
        <taxonomy> (OPTIONAL) -- taxonomy of the organism
        
        Current format: Uref
        Example: gi|187930584|ref|YP_001901071.1 Rpic_3519|hypothetical protein|Ralstonia pickettii 12J chromosome 1, complete sequence.|Ralstonia pickettii 12J|3689970|3690518|-1|Bacteria; Proteobacteria; Betaproteobacteria; Burkholderiales;Burkholderiaceae; Ralstonia.Ralstonia pickettii 12J
        """
        protein_id = "Unk"
        locus_tag = "Unk"
        product = "Unk"
        gi = "Unk"

        if "protein_id" in self.feature:
             protein_id = self.feature["protein_id"]
        if "locus_tag" in self.feature:
             locus_tag = self.feature["locus_tag"]                 
        if "product" in self.feature:
             product = self.feature["product"]  
        if "db_xref" in self.feature:
            self.extract_GI()
            if "gi" in self.feature:
                gi = self.feature["gi"]
        #Fix: if protein id is missing it will be replaced with coordinates
        if protein_id == "Unk":
            protein_id = str(self.begin) + ".." + str(self.end)
                    
        fasta_string = ">gi|%s|ref|%s %s|%s|%s|%s|%i|%i|%i|%s|%s" % (gi, protein_id, locus_tag, product, record, 
                                                                     org, self.begin, self.end, self.direction, 
                                                                     taxonomy, record_id)
        return fasta_string

    def extend_N_terminus(self, genome_seq, border, translation_table, min_extend):
        """
        Method tries to extend current gene sequence to the N-terminus, but without
        overlap with the neighbouring gene
        <genome_seq> -- full genome sequence
        <border>     -- coordinate in the genome which specifies end of the gene neighbouring
                        current gene from the N-terminus (which can be either previous
                        gene end, if current gene is direct, or next gene begin, if
                        current gene is complement)
                        
        Returns:
        <real_extension> -- amino acid sequence of the extension ('' if no extension required)
        <real_start>     -- Translation of the previous start codon ('' if no extension could be made)
        """
        ext_seq = ""
        start_seq = ""
        region_begin = 0
        region_end = 0
        if self.direction == 1:
            region_begin = border # border is end position starting from 1 (thus -1 ), we exclude it (thus +1)
            region_end = self.begin - 1
            region_begin += (region_end - region_begin) % 3 # To remain in the same reading frame

            start_seq = genome_seq[self.begin - 1 : self.begin + 2] # Sequence of the start codon (could not be 'ATG') 
        else:
            region_begin = self.end
            region_end = border - 1
            region_end -= (region_end - region_begin) % 3 # To remain in the same reading frame

            start_seq = genome_seq[self.end - 3 : self.end] # Sequence of the start codon (could not be 'ATG')
        ext_seq = genome_seq[region_begin: region_end] # Sequence between the given border and N-terminus of the current gene (in the same frame as the gene)
  
        poss_extension = udav_bio.Nucl_sequence(ext_seq, self.direction, self.feature["protein_id"])           
        poss_extension.translate(translation_table)
        real_extension = ""
        real_start = ""
        for i in reversed(range(len(poss_extension.translation))):
            if poss_extension.translation[i] == "*": # Stop codon triggered
                break
            if poss_extension.translation[i] == "M": # Possible new start!
                real_extension = poss_extension.translation[i : ]
        if len(real_extension) >= min_extend:
            new_start_seq = udav_bio.Nucl_sequence(start_seq, self.direction, self.feature["protein_id"])
            new_start_seq.translate(translation_table)
            real_start = new_start_seq.translation
            #print "Range is [%i, %i)" % (region_begin, region_end)
            #print ext_seq
            #print "EXTENSION: dir = %i, external_border = %i, gene_begin = %i, gene_end = %i" % (self.direction, border, self.begin, self.end)
            #print "           POSSIBLE_NUCL: '%s'" % poss_extension.seq
            #print "           POSSIBLE:      '%s'" % poss_extension.translation
            #print "           REAL:          '%s'" % real_extension
            #print "           START CODON is now: '%s'" % real_start
            
        return (real_extension, real_start)

class Genome_file:
    """
    Class describes a single .gbk file (NCBI database format).
    Attributes:
    str  filename      
    str  record_id     - e.g. NC_000913
    str  record_name   - e.g. Escherichia coli str. K-12 substr. MG1655, complete genome.
    int  genome_length - length of the genome
    int  gene_number   - total number of genes
    str  organism      - full name of the organism according to the file
    list taxonomy      - list with taxonomy (from widest to smallest divisions)
    list genes         - type 'Gene'
    dict features      - all features found in the record
    int  CDS_num       - number of CDS-type records
    int  pseudo_num    - number of pseudegenes
    bool circular      - 'True' if the record was circular (like plasmid)
    str  sequence      -  If not 'None', sequence of the record (small letters)     
    """
    def __init__(self, filename):
        if os.path.isfile(filename):
            self.filename = filename
        else:
            print ("File %s cannot be opened. Object will be corrupted!" % filename)
            self.filename = None
        self.record_id = None
        self.record_name = None
        self.genome_length = None
        self.gene_number = None
        self.organism = None
        self.taxonomy = None
        self.genes = None
        self.features = None
        self.CDS_num = 0
        self.pseudo_num = 0
        self.circular = None
        self.sequence = None
        self.header_strings = None
        self.database_source = None

    def obtain_data(self, huge_genome, nucl_fasta = None):
        """
        This method fills the list of genes <self.genes> with the data from the
        file <self.filename>. Check variable <req_features> into the script body:
        it determines which features are considered genes.
        - Return values -
        True  -- if worked properly
        False -- if failed to open file
        """
        if self.filename == None:
            return False # Error
        self.genes = list()
        self.header_strings = list()
        self.gene_number = 0
        input_file = open (self.filename)
        global_key = ""
        feature_key = ""
        curr_qualifier = ""
        req_features = {"CDS":1, "tRNA":1, "rRNA":1, "gene":1}
        counted_features = {"CDS":1, "tRNA":1, "rRNA":1}
        header_finished = False
        self.features = dict()
        fst_fasta_string_written = False
        
        for s in input_file:
            s = s.strip("\n")
            if header_finished == False:
            #------------------------------- 1) Header part --------------------------
            #  [0:12]
            #<---------->
            #LOCUS       NC_019815             821813 bp    DNA     circular CON 11-JUN-2013
            #DEFINITION  Candidatus Kinetoplastibacterium crithidii (ex Angomonas deanei
            #            ATCC 30255), complete genome.
            #ACCESSION   NC_019815
            #VERSION     NC_019815.1  GI:429462518
            #DBLINK      Project: 183630
            #            BioProject: PRJNA183630
            #KEYWORDS    .
            #SOURCE      Candidatus Kinetoplastibacterium crithidii (ex Angomonas deanei
            #            ATCC 30255)
            #  ORGANISM  Candidatus Kinetoplastibacterium crithidii (ex Angomonas deanei
            #            ATCC 30255)
            #            Bacteria; Proteobacteria; Betaproteobacteria;
            #            Kinetoplastibacterium.
                self.header_strings.append(s)
                key = s[0:12]
                data = s[12:]
             
                if key.strip(" ") != "":
                    global_key = s[0:12].strip(" ")                    
 
                if global_key == "LOCUS":
                    parts = data.split()
                    self.genome_length = int(parts[1])
                    if len(parts) == 6: # # FIX: 2.12 GenPept type of record 
                        molecule_type = parts[3]
                    if len(parts) == 7: # GenBank type of recored
                        molecule_type = parts[4]

                    if molecule_type == "circular":
                        self.circular = True
                    elif molecule_type == "linear":
                        self.circular = False
                    else:
                        print ("WARNING: unknown type of molecule: '%s'" % data)
                        
                if global_key == "DEFINITION":
                    if self.record_name == None:
                        self.record_name = data
                    else:
                        self.record_name += " " + data
                if global_key == "ACCESSION":
                    if self.record_id == None: # FIX: 2.11
                        self.record_id = data
                    else:
                        self.record_id += " " + data
                #if global_key == "SOURCE":                    
                #    if self.organism == None:
                #        self.organism = data
                #    else:
                #        self.organism += " " + data

                if global_key == "DBSOURCE": # FIX: 2.12
                    #DBSOURCE    accession CP000964.1
                    #DBSOURCE    embl accession CR936257.1 # FIX: 2.13
                    data = data.replace("embl ", "")
                    self.database_source = data.split(" ")[1]

                if global_key == "ORGANISM":         
                    if (data.count(";") != 0 or data[-1] == ".") and self.organism != None:                        
                        #This is taxonomy string
                        if self.taxonomy == None:                
                            self.taxonomy = data
                        else:
                            self.taxonomy += " " + data 
                    else:
                        #This is organism string
                        if self.organism == None:
                            self.organism = data
                        else:
                            self.organism += " " + data
                    
                if global_key == "FEATURES":
                    header_finished = True                    
                    if self.taxonomy == None:
                        print ("WARNING: Taxonomy is empty. Organism: %s, record ID: %s" % (self.organism, self.record_id))
                        self.taxonomy = list()
                    else:
                        taxa_list = self.taxonomy.strip(" .").split (";")                    
                        for i in range(len(taxa_list)):
                            taxa_list[i] = taxa_list[i].strip(" ")
                        taxa_list.append(self.organism)
                        self.taxonomy = taxa_list
                continue                  
            #------------------------------- 2) Features part ------------------------
            key = s[:21].strip(" ")
            data = s[21:]
            #if key[0:6] == "CONTIG": # Record red, no genome will be here
            #    break
            if key[0:6] == "ORIGIN": # Genome sequence began
                self.sequence = list()
                continue
            if s[0:10].strip().isdigit(): #This only happens if the ORIGIN string reached
                #ORIGIN      
                #        1 acgacacaga ccgaggtagg acgtggcaca ggacgaagag ctctcccggg tctggggtca
                #       61 cgtggtgacc acgctcgagg agagcccgga catcacgcag cgtcagctcg cgttcgtccg
                # <...>
                #  3526441 a  
                #//
                data = "".join(s[10:].strip().split(" "))
                if not huge_genome:                    
                    self.sequence.append(data)
                else:
                    if nucl_fasta != None: # FIX: v.2.23
                        if fst_fasta_string_written == False:
                            nucl_fasta.write("%s\n" % self.get_nucleotide_fasta_name())
                            fst_fasta_string_written = True
                        nucl_fasta.write("%s\n" % data) 
                continue

            if key != "":
                if key in self.features:
                    self.features[key] += 1
                else:
                    self.features[key] = 1 
                feature_key = key          
                if feature_key in req_features:                             
                    new_gene = Gene(data.strip(), feature_key)
                    curr_qualifier = ""
                    self.genes.append(new_gene)
                    if feature_key in counted_features: #FIX: version 2.19, correct gene number counting
                        self.gene_number += 1
            if (key == "") and (feature_key in req_features):
                curr_qualifier = self.genes[-1].try_get_feature(data, curr_qualifier)
                if curr_qualifier == "":                                     
                    self.genes[-1].coordinates += data.strip()          
        input_file.close() 

        #---------------------------------- Pseudogene search----------------------
        locuses = dict()
        for i in range(len(self.genes)):
            self.read_coordinates(self.genes[i])
            if "product" in self.genes[i].feature: # Fix of idiotic usage of | symbol
                if self.genes[i].feature["product"].count("|") > 0:
                    print ("Fixing buggy usage of symbol '|' inside the product! Gene begin = %i" % self.genes[i].begin)
                    self.genes[i].feature["product"] = self.genes[i].feature["product"].replace("|", "_")
         
            if self.genes[i].gene_type == "CDS":
                self.CDS_num += 1
                if not "translation" in self.genes[i].feature:
                    self.pseudo_num += 1
                    self.genes[i].gene_type = "pseudo" 
            else:
                if "pseudo" in self.genes[i].feature:
                    self.pseudo_num += 1
                    self.genes[i].gene_type = "pseudo" 
                           
        length = len(self.genes)

        known = dict()
        for i in range(len(self.genes)):
            curr_type = self.genes[i].gene_type
            if curr_type == "CDS" or curr_type == "tRNA" or curr_type == "rRNA":
                known[self.genes[i].coordinates] = 1

        for i in range(len(self.genes)):
            if self.genes[i].gene_type == "gene":
               if not self.genes[i].coordinates in known:
                   self.genes[i].gene_type = "other"
        
        i = 0
        while i < length:                    
            if self.genes[i].gene_type == "gene":
                self.genes.pop(i)
                length -= 1
            else:
                i += 1          

        #---------------------------------- Sequence concatenation ----------------------
        if self.sequence != None:
            self.sequence = "".join(self.sequence)

        #---------------------------------- Fixing of protein_id ------------------------
        for i in range(len(self.genes)):
            self.genes[i].fix_protein_id()

        if (nucl_fasta != None) and huge_genome:
            nucl_fasta.write("\n")

        return True
    
    def get_intergene_dist(self, dist_direct, dist_reverse): #FIX: 2.14
        """
        Method calculates intergene distance between all genes in current gbk file
        and fills <dist_direct> and <dist_reverse> dictionaries
        <dist> is a dictionary of dictionaries:
         * dist["PROT-RNA"] : Gene.gene_type 'CDS' vs 'tRNA', 'rRNA'
         * dist["RNA-RNA"]  : Gene.gene_type 'tRNA', 'rRNA' vs 'tRNA', 'rRNA'
         * dist["PROT-PROT"]: Gene.gene_type 'CDS' vs 'CDS'
        """
        i = 0
        while i < len(self.genes) - 1:
            curr = self.genes[i] 
            next = self.genes[i + 1]
            distance = max(curr.begin, next.begin) - min(curr.end, next.end) - 1
            curr_prot_id = "Unk"
            next_prot_id = "Unk"
            feature_to_report = "locus_tag" # FIX: 2.22 - now not only protein_id can be reported
            if feature_to_report in curr.feature:
                curr_prot_id = curr.feature[feature_to_report]
            if feature_to_report in next.feature:
                next_prot_id = next.feature[feature_to_report]
            
            info_str = "%i\t%s\t%s\t%s\t%s..%s\t%s" % (distance, self.record_id, curr_prot_id, next_prot_id, curr.begin, curr.end, self.sequence[next.begin-1:curr.end])
            subdict = None

            # 1) Getting subdictionary information
            if (curr.gene_type == "CDS") and (next.gene_type == "CDS"):
                subdict = "PROT-PROT"
            elif (curr.gene_type == "CDS") and (next.gene_type in ["tRNA", "rRNA"]):
                subdict = "PROT-RNA"
            elif (curr.gene_type in ["tRNA", "rRNA"]) and (next.gene_type == "CDS"):
                subdict = "PROT-RNA"
            elif (curr.gene_type in ["tRNA", "rRNA"]) and (next.gene_type in ["tRNA", "rRNA"]):
                subdict = "RNA-RNA"
            else:
                subdict = "OTHER"
            #print "DIST = %i, sub = %s" % (distance, subdict)
            #print "CURR: begin = %i, end = %i" % (curr.begin, curr.end)
            #print "NEXT: begin = %i, end = %i" % (next.begin, next.end)
            #raw_input()

            # 2) Getting direction information
            if curr.direction == next.direction:
                if not distance in dist_direct[subdict]:
                    dist_direct[subdict][distance] = 0
                dist_direct[subdict][distance] += 1
                if distance == -1:
                    dist_direct[-1].append(info_str)
                if distance == -4:
                    dist_direct[-4].append(info_str)
            else:
                if not distance in dist_reverse[subdict]:
                    dist_reverse[subdict][distance] = 0
                dist_reverse[subdict][distance] += 1
                if distance == -4:
                    dist_reverse[-4].append(info_str)
            i += 1
        
    def print_table(self, output_filename, fields):
        """
        This method prints table with information
        <output_filename> -- name of the file with the table which should be created
        <fields>          -- list of elements with the data qualifiers ("protein_id" etc)       
        """
        output_file = open(output_filename, "w")
        output_file.write("#Chromosome table for the file '%s'\n" % self.filename)
        output_file.write("#Record ID:\t%s\n" % self.record_id)
        output_file.write("#Record name:\t%s\n" % self.record_name)
        output_file.write("#Organism:\t%s\n" % self.organism)
        output_file.write("#Taxonomy:\t")
        for t in self.taxonomy:
            output_file.write("%s; " % t)
        output_file.write("\n")  
        output_file.write("#Features in this file:\n")
        output_file.write("#Feature\tNumber\n")
        for key in self.features.keys():
            output_file.write("#%s\t%s\n" % (key, self.features[key]))
        output_file.write("#type\tbegin\tend\torientation")
        for f in fields:                                                           
           output_file.write("\t" + f)                                                   
        output_file.write("\n")                                                    
        
        for g in self.genes:
            output_file.write(g.extract_data(fields))
        output_file.close()

    def read_coordinates(self, feature_obj):
        """
        Read coordinates of the gene (non-proper format will result in a return of [0, 0, 0])
        Proper formats are:
        1) 7548..8180
        2) complement(406..786)
        3) complement(join(1390189..1390483,1390483..1391310,1391310..1393681))
        4) order instead of join, as: order(147423..148106,148108..149580)
        Supported formats are:
        1) 2516438 (a single begin coordinate, no ending); also inside join or other operators

        In case of multiple coordinate pairs, first value will be the beginning and last value 
        will be the end of the given <feature_obj> (Gene or sister objects)
        """
        string = feature_obj.coordinates
        begin = 0
        end = 0
        direction = 1
        string_ready = re.sub("[\>\<\(\)]", "", string)
        #string_ready = string.translate(None, "><()") # Removes symbols (Python 2.7)
        
        if string_ready.count("complement") > 0:
            string_ready = string_ready.replace("complement", "")
            direction = -1

        if string_ready.count("join") > 0:
            string_ready = string_ready.replace("join", "")
        if string_ready.count("order") > 0:
            string_ready = string_ready.replace("order", "")

        coordinates = string_ready.split(",")
        begin = coordinates[0].split("..")[0]
        for c in coordinates:
            second_split = c.split("..")
            if len(second_split) == 2:
                end = second_split[1]
            else:
                end = begin
                #print "WARNING: strange coordinate format. File: %s, organism %s" % (self.filename, self.organism)
                #print "Original string: " + string
                #print "Coordinates array: "
                #print coordinates  
        try: # FIX: 2.24 - found case like join(OM066871.1:29965..31731,15583..20330) !!11
            feature_obj.begin = int(begin)
        except ValueError:
            feature_obj.begin = 0
        try:
            feature_obj.end = int(end)
        except ValueError:
            feature_obj.end = 0
        feature_obj.direction = direction            

    def extend_protein_genes(self, translation_table, extended_file, min_extend):
        n = 0
        curr_taxonomy = ""
        for i in range(len(self.taxonomy) - 1):
            curr_taxonomy += self.taxonomy[i] + "; "
        curr_taxonomy = curr_taxonomy.strip("; ")

        for i in range(1, len(self.genes) - 1):
            if self.genes[i].gene_type == "CDS":
                border = 0
                if self.genes[i].direction == 1: # Direct chain           
                    border = self.genes[i - 1].end
                else: # Complement chain
                    border = self.genes[i + 1].begin 
                (result, new_start) = self.genes[i].extend_N_terminus(self.sequence, border, translation_table, min_extend)
                if result != "":
                    n += 1

                    extended_file.write(self.genes[i].get_fasta(self.record_name, self.record_id, self.organism, curr_taxonomy) + "\n")
                    extended_file.write("%s%s%s\n" % (result, new_start, self.genes[i].feature["translation"][1:]))
                    print ("%s\t%i..%i\t%s" % (self.record_id, self.genes[i].begin, self.genes[i].end, self.organism))
        return n

    def get_sequence(self, chars_in_string):
        seq_strings = list()
        seq_strings.append("")
        if self.sequence != None:     
            c = 0  
            for char in self.sequence:
                if c >= chars_in_string:
                    seq_strings.append("")
                    c = 0
                c += 1
                seq_strings[-1] += char                    
        return seq_strings
 
    def print_gbk(self, filename, feature_order, mode = "w"):
        new_file = open(filename, mode)
        for string in self.header_strings: #----------- 1) Header strings
            new_file.write("%s\n" % string)

        for gene in self.genes: #---------------------- 2) Genes and features                               
            new_file.write("%s\n" % get_global_string("gene", gene.coordinates))
            new_file.write("%s\n" % get_local_string("locus_tag", gene.feature["locus_tag"], "str"))
            new_file.write("%s\n" % get_global_string(gene.gene_type, gene.coordinates))
            for feature_name in feature_order:
                if feature_name in gene.feature:
                    feature_type = "int"
                    try:
                        int(gene.feature[feature_name])
                    except ValueError:
                        feature_type = "str"
                    new_file.write("%s\n" % get_local_string(feature_name, gene.feature[feature_name], feature_type))

        new_file.write("ORIGIN\n") #------------------- 3) Genome sequence
        new_file.write("%s\n" % get_gbk_sequence(self.sequence))
        new_file.write("//\n")
        new_file.close()  

    def get_nucleotide_fasta_name(self):
        fasta_name = ">ref|%s|%s|%s|%s" % (self.record_id.split(" ")[0], self.organism, "; ".join(self.taxonomy), self.record_name)
        return fasta_name

    def print_nucleotide(self, nucleotide_file):
        """
        Method will append information about nucleotide sequence of current
        genome to given <nucleotide_file>
        """
        if self.sequence == None:
            print ("WARNING: record from '%s' does not contain sequence data!" % self.filename)
        else:
            #nucleotide_file.write(">%s|%ibp|%s\n" % (self.record_id.replace(" ", "_"), len(self.sequence), self.organism.replace(" ", "_")))
            # Olesya's format for COGNAT database:
            #ref|AAXS01000001.1|candidate division TM7 genomosp. GTL1|Bacteria; Candidatus Saccharibacteria; unknown; unknown; unknown; TM7|Candidate division TM7 genomosp. GTL1 ctg1, whole genome shotgun sequence.
            nucleotide_file.write("%s\n" % self.get_nucleotide_fasta_name())
            seq_strings = self.get_sequence(60)
            for string in seq_strings:
                nucleotide_file.write("%s\n" % string)
            nucleotide_file.write("\n")

    def get_taxonomy(self):
        curr_taxonomy = ""
        curr_taxonomy_tab = ""
        for i in range(len(self.taxonomy) - 1):
            curr_taxonomy += self.taxonomy[i] + "; "
            curr_taxonomy_tab += self.taxonomy[i] + "\t"
        curr_taxonomy = curr_taxonomy.strip("; ")
        curr_taxonomy_tab = curr_taxonomy_tab.strip("\t")
        return (curr_taxonomy, curr_taxonomy_tab)

    def print_proteins(self, protein_file, curr_assembly = None, req_id = None, req_id_type = "protein_id"):
        """
        Method will append information about protein sequences annotated
        in the current genome with the 'CDS' qualifier to given <protein_file>.
 
        URef record data (last field) will contain given <curr_assembly>: FIX 2.17

        If <req_id> dict is given, only proteins from it will be printed (feature
        with the name <ref_id_type> is used).
        """
        (curr_taxonomy, curr_taxonomy_tab) = self.get_taxonomy()
        for g in self.genes:
            if g.gene_type == "CDS":
                take_gene = True
                fasta_data = g.get_fasta(self.record_name, self.organism, curr_taxonomy, self.record_id)
                if curr_assembly != None:
                    fasta_data += " @ %s" % curr_assembly
                sequence = ""
                if "translation" in g.feature:
                    sequence = g.feature["translation"]
                else:
                    print ("EXCEPTION: protein sequence is missing")
                    print ("Organism: %s, record: %s, gene begin: %i" % (self.organism, self.record_id, g.begin))

                if req_id != None:                    
                    if req_id_type in g.feature:
                        if not g.feature[req_id_type] in req_id:
                            take_gene = False
                    else:
                        take_gene = False
                if take_gene:
                    protein_file.write(fasta_data + "\n")
                    protein_file.write(sequence + "\n\n")

    def get_operons(self, required_genes, operon_output, table_output, id_type):
        operons = dict()
        o_num = 0
        for i in range(len(self.genes)):            
            if (self.genes[i].gene_type == "CDS") and (id_type in self.genes[i].feature):                              
                if self.genes[i].feature[id_type] in required_genes: # Current gene is required
                    o_num += 1
                    curr_dir = self.genes[i].direction
                    strings_to_print = dict()
                    j = i - 1

                    while j >= 0: #Going down
                        if (self.genes[j].direction == curr_dir) and (self.genes[j].gene_type == "CDS"):
                            if id_type in self.genes[j].feature:
                                operons[self.genes[j].feature[id_type]] = True                                 
                                try:
                                    strings_to_print[j] = "%s\t%s\t%i\t%s\t%s\t%s\t%i\t%i\t%i\n" % (self.record_id, self.organism, o_num,
                                                                         self.genes[j].feature["protein_id"], self.genes[j].feature["locus_tag"],
                                                                         self.genes[j].feature["product"],
                                                                         self.genes[j].begin, self.genes[j].end, self.genes[j].direction)
                                except:
                                    print ("Likely one of required features ('protein_id', 'locus_tag' and 'product') were missing!")
                                    print (self.genes[j].feature)
                        else:
                            break
                        j -= 1
                    j = i                    
                    while j < len(self.genes): #Going up
                        if (self.genes[j].direction == curr_dir) and (self.genes[j].gene_type == "CDS"):
                            if id_type in self.genes[j].feature:
                                operons[self.genes[j].feature[id_type]] = True
                                try:
                                    strings_to_print[j] = "%s\t%s\t%i\t%s\t%s\t%s\t%i\t%i\t%i\n" % (self.record_id, self.organism, o_num,
                                                                         self.genes[j].feature["protein_id"], self.genes[j].feature["locus_tag"],
                                                                         self.genes[j].feature["product"],
                                                                         self.genes[j].begin, self.genes[j].end, self.genes[j].direction)
                                except:
                                    print ("Likely one of required features ('protein_id', 'locus_tag' and 'product') were missing!")
                                    print (self.genes[j].feature)
                        else:
                            break
                        j += 1
                    reverse_sort = False
                    if curr_dir == -1:
                        reverse_sort = True
                    for key in sorted(strings_to_print.keys(), reverse=reverse_sort, key=lambda t: t):
                        table_output.write(strings_to_print[key])
                    table_output.write("#\n")
 
        if len (operons.keys()) != 0:
            operon_output.write("#* %s | %s | %s *\n" % (self.organism, self.record_name, self.record_id))
            self.print_proteins(operon_output, None, operons, id_type)            
            operon_output.write("\n")           

    def get_old_locus_data(self, old_locus_dict):
        for gene in self.genes:            
            if (gene.gene_type == "CDS") and ("locus_tag" in gene.feature):
                if not gene.feature["locus_tag"] in old_locus_dict:
                    if "old_locus_tag" in gene.feature:
                        old_locus_dict[gene.feature["locus_tag"]] = (gene.feature["old_locus_tag"], gene.feature["protein_id"], self.record_id)
                    else:
                        old_locus_dict[gene.feature["locus_tag"]] = (None, gene.feature["protein_id"], self.record_id)
                else:
                    print ("WARNING: gene locus '%s' found more than one time in CDS records!" % gene.feature["locus_tag"])

    def print_pseudo(self, fasta_file):
        for g in self.genes:
            if g.gene_type == "pseudo":                
                (curr_taxonomy, curr_taxonomy_tab) = self.get_taxonomy()
                fasta_data = g.get_fasta(self.record_name, self.organism, curr_taxonomy, self.record_id)
                coord = g.coordinates.replace("complement","")
                curr_begin = int(coord.split("..")[0].strip("<>()"))
                curr_end = int(coord.split("..")[-1].strip("<>()"))
                sequence = self.sequence[curr_begin - 1:curr_end]
                fasta_file.write("%s\n%s\n\n" % (fasta_data, sequence))

def read_Prodigal_output(filename):
    """
    Method reads Prodigal .sco format
    """
    seqnum_to_genes = dict()
    curr_seqnum = None
    sco_file = open(filename, "r")
    for string in sco_file:
        ## Sequence Data: seqnum=1;seqlen=70790;seqhdr="Aenigmarchaeota archaeon JGI 0000106-F11"
        ## Model Data: version=Prodigal.v2.6.2;run_type=Single;model="Ab initio";gc_cont=37.42;transl_table=11;uses_sd=1
        #>1_3_479_-
        string = string.strip()
        if string[0:15] == "# Sequence Data":
            data = string[17:]
            fields = data.split(";")
            curr_seqnum = int(fields[0].split("=", 1)[1])
            seqnum_to_genes[curr_seqnum] = list()
        if string[0] == ">":
            fields = string.strip(">").split("_")
            if len(fields) != 4:
                print ("FATAL ERROR: format of the .sco file is not recognized! String:")
                print ("seqnum=%i" % curr_seqnum)
                print (string)
                print (fields)
                sys.exit()
            temp_gene_id = "%s_%i_%s" % (prefix, curr_seqnum, fields[0])
            begin = int(fields[1])
            end = int(fields[2])
            coordinates = "%i_%i" % (begin, end)
            gene_type = "CDS"
            direction = None
            if fields[3] == "-":
                direction = -1
            else:
                direction = 1
            new_gene = Gene(coordinates, gene_type, begin, end, direction)
            new_gene.feature["protein_id"] = temp_gene_id
            
            seqnum_to_genes[curr_seqnum].append(new_gene)    
    sco_file.close()
    return seqnum_to_genes

def empty_gbk_file(filename):
    """
    Method checks if this gbk file does not contain any mapped contigs
    """
    do_not_store_nucleotide_sequence = True
    new_genome = Genome_file(filename)       
    new_genome.obtain_data(do_not_store_nucleotide_sequence)
    result = False
    if len(new_genome.genes) == 0: # This file is a simple contig with no genes
        result = True
    return result

def split_gbk(curr_filename, output_dir, empty_dir, remove_empty = False, list_under_source = False):
    """
    Method splits this gbk into separate files.
    If <list_under_source> is False, files will be named based on the <curr_filename> prefix (before extension).
    If <list_under_source> is True, they will be named based on the SOURCE data 
    """
    i = 1
    base_name = os.path.basename(curr_filename)
    fields = base_name.split(".")
    name_prefix = ".".join(fields[0:len(fields) - 1])
    new_filename = os.path.join(empty_dir, base_name)

    # ---- FIRST RUN FIX: 1.5
    org = "Unk"
    if remove_empty == True:
        curr_file = open (curr_filename, "r")
        translation_found = False
        for string in curr_file:                   
            if org == "Unk":
                res = re.search("^SOURCE +([^\n]+)", string)
                if res != None:
                    org = res.group(1)
                    print (org)

            if re.search("\/translation\=", string) != None:
                translation_found = True
                break
        curr_file.close()
        if translation_found == False:
            print ("File %s does not contain any protein records. It is to be moved away!" % base_name)
            if os.path.isfile(new_filename):
                os.remove(new_filename)
            os.rename(curr_filename, new_filename)
            return (100500, 0, name_prefix, org)

    # ---- SECOND RUN
    curr_file = open (curr_filename, "r")
    c = 0 # Number of 'bad' contig files
    f = 0 # Number of 'good' contig files 
    new_file = None
    new_name = None

    for string in curr_file:
        if org == "Unk": #FIX: 1.6
            res = re.search("^SOURCE +([^\n]+)", string)
            if res != None:
                org = res.group(1)

        if string[0:5] == "LOCUS":
           if new_file != None:
               new_file.close()

               if empty_gbk_file(new_name):
                   c += 1                   
                   if remove_empty:
                       os.remove(new_name)
                       i -= 1
               else:
                   f += 1
           new_name = None
           if list_under_source == False:
               new_name = os.path.join(output_dir, name_prefix + "_" + str(i) + ".gbk")
           else: #FIX: 2.24 version (for reading any gbk file)
               splitted_string = re.split("\s+", string)
               new_name = os.path.join(output_dir, splitted_string[1] + ".gbk")
               
           new_file = open (new_name, "w")
           i += 1
        if new_file != None:
            new_file.write(string)
    new_file.close()
    if empty_gbk_file(new_name): #FIX: 1.3 (last empty file was not removed)
        c += 1                   
        if remove_empty:
            os.remove(new_name)
            i -= 1
    else:
        f += 1                   
    curr_file.close()
    return (c, f, name_prefix, org)

def print_dist_file(filename, dist, begin = None, end = None):
    """
    Prints information from the <dist_direct> or <dist_reverse> dict to the
    file with the <filename>.
    If <begin> and <end> are given, prints 'flat' statistics from one to another only
    """
    output_file = open (filename, "w")
    dist_types = ["PROT-PROT", "PROT-RNA", "RNA-RNA", "OTHER"]
    #dist_types = ["PROT-PROT"]
    for dist_type in dist_types:
        output_file.write("# Distance type: %s\n" % dist_type)
        sorted_keys = sorted(dist[dist_type].keys())
        if len(sorted_keys) == 0:
            continue
        if begin == None:
            begin = sorted_keys[0]
        if end == None:
            end = sorted_keys[-1] + 1

        for key in sorted_keys:
            if key < begin:
                output_file.write("%s\t%s\n" % (key, dist[dist_type][key]))

        for i in range(begin, end):
            value = 0
            if i in dist[dist_type]:
                value = dist[dist_type][i]
            output_file.write("%s\t%s\n" % (i, value))

        for key in sorted_keys:
            if key >= end:
                output_file.write("%s\t%s\n" % (key, dist[dist_type][key]))

    distances = [-4]
    for d in distances:
        if d in dist:
            output_file.write("# Information about genes with this distance: %i\n" % d)
            for element in dist[d]:
                 output_file.write("%s\n" % element)
    output_file.close()

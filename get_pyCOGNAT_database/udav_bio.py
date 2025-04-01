"""
Module for manipulation with nucleotide sequences
------- Version: 1.1
        1.1  * Proper class & function documentation added
"""

import sys, os

class Nucl_sequence:
    """
    Class describing nucleotide sequence.
    Attributes:
    str seq         - nucleotide sequence
    str ID          - some kind of sequence identifier
    str translation - protein sequence coded by this nucleotide sequence
    """
    def __init__(self, seq, direction, ID):
        self.seq = seq.upper()
        self.ID = ID
        if direction == -1: # Sequence must be reversed
            self.reverse()
        self.translation = ""
         
    def translate (self, translation_table):
        """
        Builds translation and fills corresponding attribute using given <translation_table>
        in the form of a dictionary with codons as keys and amino acids as values
        """
        self.translation = list()
        for i in range (0, len(self.seq) - 3 + 1, 3):
            curr_codon = self.seq[i:i+3]
            if len(curr_codon) != 3:
                print ("FATAL ERROR: codon %s was detected at %i symbol of protein %s!" % (curr_codon, i, self.ID))
            curr_aa = "X"
            if curr_codon in translation_table:
                curr_aa = translation_table[curr_codon]
            else:
                print ("WARNING: Codon '%s' not found in translation table (protein %s)!" % (curr_codon, self.ID))
            self.translation.append(curr_aa)
        self.translation = "".join(self.translation)

    def reverse (self):
        """
        Method reverses nucleic acid sequence by getting its complement and making direct reverse
        """
        complement = {"A":"T", "T":"A", "G":"C", "C":"G"}
        rev = ""
        for s in self.seq:
            if s in complement:
                rev += complement[s]
            else:
                print ("FATAL ERROR: unknown nucleotide %s was detected in protein %s!" % (s, self.ID))
                sys.exit()
        rev = rev[::-1] #Extended slice syntax: reversing order
        self.seq = rev
       
def read_translation_table (filename):
    """
    Method reads translation table in NCBI format from the <filename>
    """
    table_hash = dict()    
    aa = ""
    base = list()
    
    table_file = open (filename, "r")
    for string in table_file:
        string = string.strip()
        name = string.split(" = ")[0].strip()
        data = string.split(" = ")[1].strip()      
        if name == "AAs":
            aa = data
        if name == "Base1" or name == "Base2" or name == "Base3":
            base.append(data)
    table_file.close()

    if len(base) != 3:
        print ("ERROR: translation table could not be written correctly!")
    else:
        for i in range(len(aa)):
            curr_codon = base[0][i] + base[1][i] + base[2][i]
            table_hash[curr_codon] = aa[i]
            table_hash[curr_codon.lower()] = aa[i]
    return table_hash    

class Chromosome_table:
    """
    This class describes a chromosome table created by the <read_Genome_release.py> script
    Main variable:
    genes - list of genes stored as hashes with keys corresponding 
            to values in chromosome table
    """
    def __init__(self, filename):
        self.genes = list()
        input_table = open(filename, "r")        
        for string in input_table:
            string = string.strip()
            if len(string) == 0:
                continue
            if string[0] == "#":
                continue
            fields = string.split("\t")
            if len(fields) < 4:
                print ("FATAL ERROR: Unexpected chromosome table format. Check that the following")
                print ("             fields are included: 'type', 'begin', 'end', 'orientation'.")
                print ("             Also these should be there: 'locus_tag', 'protein_id', 'gene', 'product'")
                print ("Number of fields: %i" % len(fields))
                print ("String = %s" % string)
                print (fields)
                sys.exit()
            new_gene = dict()
            new_gene["type"] = fields[0]
            new_gene["begin"] = int(fields[1])
            new_gene["end"] = int(fields[2])
            new_gene["orientation"] = int (fields[3])
            if len(fields) > 4:
                new_gene["locus_tag"] = fields[4]
            else:
                new_gene["locus_tag"] = "undef"
            if len(fields) > 5:
                if fields[5] == "":
                    new_gene["protein_id"] = "%i..%i" % (new_gene["begin"], new_gene["end"])
                else:
                    new_gene["protein_id"] = fields[5]
            else:
                new_gene["protein_id"] = "%i..%i" % (new_gene["begin"], new_gene["end"])
            if len(fields) > 6:
                new_gene["gene"] = fields[6]
            else:
                new_gene["gene"] = "undef"
            if len(fields) > 7:
                new_gene["product"] = fields[7]
            else:
                new_gene["product"] = "undef"
            self.genes.append(new_gene)
        input_table.close()

    def get_dist(self, left_gene, right_gene):
        return right_gene["begin"] - left_gene["end"] - 1

    def return_operon(self, target, operon_dist):
        operon = list()                         
        direction = 0
        #print "Searching for the operon\n\n"
        for i in range(len(self.genes)):                       
            if self.genes[i]["protein_id"] == target:                
                leftmost_found = False
                rightmost_found = False
                direction = self.genes[i]["orientation"]
                l = i
                r = i
                while (not leftmost_found) or (not rightmost_found):
                    if not leftmost_found: # Searching to the left; l is decreasing
                        l -= 1
                        if l <= 0: # 5'-terminus of the record riched 
                            leftmost_found = True
                        else:      
                            curr_dist = self.get_dist(self.genes[l], self.genes[l + 1])
                            #print "i = %i, l = %i. Left distance = %s, operon_dist = %i" % (i, l, curr_dist, operon_dist)
                            if (curr_dist > operon_dist) or (self.genes[l]["orientation"] != self.genes[i]["orientation"]):
                                leftmost_found = True
                                #print "Left found!"
                    if not rightmost_found: # Searching to the right; r is increasing
                        r += 1
                        if r >= len(self.genes): # 3'-terminus of the record riched 
                            rightmost_found = True
                        else:      
                            curr_dist = self.get_dist(self.genes[r - 1], self.genes[r])
                            #print "i = %i, r = %i. Right distance = %s, operon_dist = %i" % (i, r, curr_dist, operon_dist)
                            if (curr_dist > operon_dist) or (self.genes[r]["orientation"] != self.genes[i]["orientation"]):
                                rightmost_found = True
                                #print "Right found!"                
                for j in range (l + 1, r):
                    operon.append(self.genes[j])        

        if (direction == -1) and (len(operon) > 1):
            operon.reverse()
        return operon
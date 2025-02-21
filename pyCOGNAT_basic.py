# -*- coding: utf-8 -*-
import os, sys, re, codecs, platform, time
import copy
import tkinter
import tkinter.font

def check_filename(filename):
    if re.search("[^\w]", filename) != None:
        return False
    else:
        return True

class TextFrameWithLabel(tkinter.Frame):
    def __init__(self, parent, padding, background, text_color, label_text, label_font, back_and_thick = None):
        tkinter.Frame.__init__(self, parent)
        if back_and_thick != None:
            self.configure(highlightbackground = back_and_thick[0], highlightthickness = back_and_thick[1])
        self.p = padding
        self.back = background
        self.text_color = text_color
        self.label = None
        self.text_widget = None
        self.check = None
        self.wrap_text = tkinter.BooleanVar()
        self.wrap_text.set(True)
        self.font_size = 8
        self.create_UI(label_text, label_font)

    def create_UI(self, label_text, label_font):
        self.grid_columnconfigure(0, weight = 1)
        self.grid_rowconfigure(1, weight = 1)

        panel = tkinter.Frame(self)
        panel.grid_columnconfigure(0, weight = 1)
        panel.grid_rowconfigure(0, weight = 1)
        self.label = tkinter.Label(panel, text = label_text, font = label_font)
        self.label.grid(row = 0, column = 0, sticky = "NSW", padx = self.p, pady = self.p)
        self.check = tkinter.Checkbutton(panel, text = "Wrap text", variable = self.wrap_text, command = self.wrap_configure)
        self.check.grid(row = 0, column = 1, sticky ="NSE", padx = self.p, pady = self.p)
        panel.grid(row = 0, column = 0, columnspan = 2, sticky = "NSEW")

        text_scr_y = tkinter.Scrollbar(self, orient = tkinter.VERTICAL)
        text_scr_y.grid(row = 1, column = 1, sticky = "NSEW")
        text_scr_x = tkinter.Scrollbar(self, orient = tkinter.HORIZONTAL)
        text_scr_x.grid(row = 2, column = 0, sticky = "NSEW")
        self.normal_font = tkinter.font.Font(self.text_widget, ("Courier New", self.font_size))
        self.text_widget = tkinter.Text(self, state = tkinter.NORMAL, font = self.normal_font,
                                        yscrollcommand = text_scr_y.set, xscrollcommand = text_scr_x.set, wrap = tkinter.WORD)
        self.text_widget.grid(row = 1, column = 0, sticky = "NSEW")
        text_scr_y.configure(command = self.text_widget.yview)
        text_scr_x.configure(command = self.text_widget.xview)

        if platform.system() == "Linux":
            self.text_widget.bind("<Control-Button-4>", self.change_font)
            self.text_widget.bind("<Control-Button-5>", self.change_font)
        else:
            self.text_widget.bind("<Control-MouseWheel>", self.change_font)

    def wrap_configure(self):
        if self.wrap_text.get() == False:
            self.text_widget.configure(wrap = tkinter.NONE)
        else:
            self.text_widget.configure(wrap = tkinter.WORD)

    def change_font(self, event):
        change = 0
        if platform.system() == "Linux":
            if event.num == 4: # Scroll up
                change = 1
            if event.num == 5: # Scroll down
                change = -1
        else:
            if event.delta > 0: # Scroll up
                change = 1
            else:
                change = -1
        if (change < 0) and (self.font_size == 1):
            return

        self.font_size += change
        self.normal_font.configure(size = self.font_size)

    def write_into_file(self, filename, use_codecs = False):
        """
        Method returns True if the were was written correctly and False if some errors occured
        """
        curr_text = self.text_widget.get(1.0, tkinter.END)
        if len(curr_text.strip()) == 0: # Text is empty
            return False
        if filename == "": # No filename provided, save was canceled
            return False
        strings = curr_text.split("\n")
        if use_codecs:
            output_file = codecs.open(filename, "w", encoding = "utf_8")
        else:
            output_file = open(filename, "w")
        for string in strings:
            output_file.write("%s\n" % string)
        output_file.close()
        return True

    def read_from_file(self, filename, use_codecs = False):
        if not os.path.isfile(filename):
            return
        self.text_widget.configure(state = tkinter.NORMAL)
        self.text_widget.delete(1.0, tkinter.END)
        if use_codecs:
            input_file = codecs.open(filename, "r", encoding = "utf_8")
        else:
            input_file = open(filename, "r")
        for string in input_file:
            self.text_widget.insert(tkinter.END, string)
        input_file.close()

        stripped_text = self.text_widget.get(1.0, tkinter.END).strip()
        self.text_widget.delete(1.0, tkinter.END)
        self.text_widget.insert(tkinter.END, stripped_text)

    def get_content_as_list(self):
        """
        Returns content of the <self.text_widget> as list of strings, each stripped.
        Empty strings and strings with '#' in the beginning are omitted
        """
        strings = self.text_widget.get(1.0, tkinter.END).split("\n")
        result_strings = list()
        for s in strings:
            s = s.strip()
            if len(s) == 0:
                continue
            if s[0] == "#":
                continue
            result_strings.append(s)
        return result_strings

    def get_content_as_dict(self, separator):
        """
        Returns content of the <self.text_widget> as a dictionary of the first 'column'
        to all other part of the string (stripped)
        Empty strings and strings with '#' in the beginning are omitted
        """
        strings = self.text_widget.get(1.0, tkinter.END).split("\n")
        result = dict()
        for s in strings:
            s = s.strip()
            if len(s) == 0:
                continue
            if s[0] == "#":
                continue
            parts = s.split(separator, 1)
            result[parts[0]] = parts[1]
        return result

    def add_string(self, string):
        self.text_widget.insert(tkinter.END, string)

    def add_tagged_string(self, string, tag_name):
        self.add_string("%s\n" % string)
        start = self.text_widget.search(string, 1.0, stopindex = tkinter.END)
        end = '%s+%dc' % (start, len(string))
        self.text_widget.tag_add(tag_name, start, end)

    def clear(self, disable = False):
        self.text_widget.configure(state = tkinter.NORMAL)
        self.text_widget.delete(1.0, tkinter.END)
        if disable:
            self.text_widget.configure(state = tkinter.DISABLED)
            self.check.configure(state = tkinter.DISABLED)

    def enable_text(self):
        self.text_widget.configure(state = tkinter.NORMAL)
        self.check.configure(state = tkinter.NORMAL)

    def disable_text(self):
        self.text_widget.configure(state = tkinter.DISABLED)
        self.check.configure(state = tkinter.DISABLED)

def clear_and_disable_entry_widget(entry_widget, do_not_disable = False):
    entry_widget.configure(state = tkinter.NORMAL)
    entry_widget.delete(0, tkinter.END)
    if not do_not_disable:
        entry_widget.configure(state = tkinter.DISABLED)

def read_plain_list(filename):
    result = list()
    plain_file = open(filename)
    for string in plain_file:
        string = string.strip()
        if len(string) == 0:
            continue
        result.append(string)
    plain_file.close()
    return result

def write_plain_list(filename, elements):
    plain_file = open(filename, "w")
    for element in elements:
        plain_file.write("%s\n" % element)
    plain_file.close()

def filter_table(initial_filename, new_filename, good_first_column):
    """
    Reads tab-separated <initial_filename>, checks which of its first column values
    occurs in <good_first_column> list and prints corresponding strings to <new_filename>.
    Comment lines (starting with '#') are taken 'as is'
    """
    req = dict()
    for element in good_first_column:
        req[element] = True
    initial_file = open(initial_filename)
    new_file = open(new_filename, "w")
    for string in initial_file:
        string = string.strip()
        if len(string) == 0:
            continue
        fields = string.split("\t", 1)
        if (fields[0] in req) or (string[0] == "#"):
            new_file.write("%s\n" % string)
    initial_file.close()
    new_file.close()

def merge_overlapping_ranges(ranges):
    """
    Given <list_of_ranges> should be a list of lists of two values - <begin> and <end>.
    This method will merge overlapping ranges
    """
    ranges.sort(key = lambda k: k[0], reverse = False)
    i = 0
    while i < len(ranges) - 1:
        overlap = min(ranges[i][1], ranges[i + 1][1]) - max(ranges[i][0], ranges[i + 1][0]) + 1
        if overlap > 0:
            new_begin = min(ranges[i][0], ranges[i + 1][0])
            new_end = max(ranges[i][1], ranges[i + 1][1])
            print ("Found overlap between ranges: [%i, %i] and [%i, %i]. New range: [%i, %i]" % (ranges[i][0], ranges[i][1], ranges[i + 1][0], ranges[i + 1][1], new_begin, new_end))
            ranges[i] = [new_begin, new_end]
            ranges.pop(i + 1)
        else:
            i += 1

def obtain_database_slice(database_path, output_path, script_path, protein_ids, max_gene_number, ibar, rbar, lbar, root):
    """
    Method reads COGNAT database and creates its copy, but only for the proteins
    specified in the <protein_ids> list and their surroundings of <max_gene_number>
    to the left and right. Progress will be reported to the <ibar>, <rbar> and <lbar>
    """
    ids_to_gbff_and_gbk = None
    try:
        sys.path.append(script_path)
        import COGNAT_basic
        ids_to_gbff_and_gbk = COGNAT_basic.map_proteins_to_gbff_and_gbk(protein_ids, database_path, ibar, root)
        del COGNAT_basic
    except ImportError:
        print ("No module 'COGNAT_basic' and/or its dependencies was found in the script directory '%s'" % script_path)
        return
    try:
        sys.path.append(script_path)
        import udav_fasta
    except ImportError:
        print ("No module 'udav_fasta' and/or its dependencies was found in the script directory '%s'" % script_path)
        return

    # 1) Obtaining information about proteins and nucleotide regions which should be taken to database slice
    gbff_to_gbk_to_req_proteins = dict() # Dict of dict of dicts
    gbff_to_gbk_to_ranges = dict()
    all_proteins = dict()
    req_and_found = dict()
    n = 0
    for protein_id in ids_to_gbff_and_gbk.keys():
        n += 1
        rbar["value"] = int(100 * n / len(ids_to_gbff_and_gbk.keys()))
        root.update()
        gbff_prefix = ids_to_gbff_and_gbk[protein_id][0]
        if not gbff_prefix in gbff_to_gbk_to_req_proteins:
            gbff_to_gbk_to_req_proteins[gbff_prefix] = dict()
            gbff_to_gbk_to_ranges[gbff_prefix] = dict()
        gbk_prefix = ids_to_gbff_and_gbk[protein_id][1]
        if not gbk_prefix in gbff_to_gbk_to_req_proteins[gbff_prefix]:
            gbff_to_gbk_to_req_proteins[gbff_prefix][gbk_prefix] = dict() # dictionary of proteins which should be taken
            gbff_to_gbk_to_ranges[gbff_prefix][gbk_prefix] = list()

        curr_prot_path = os.path.join(database_path, gbff_prefix, gbk_prefix, "p_%s.fasta" % gbk_prefix)
        (prot_records, found) = udav_fasta.read_fasta(curr_prot_path, None, None, dict(), 100500, False, None, None, None, None, "COGNAT")
        prot_records.sort(key = lambda k: k.gene_begin, reverse = False)
        for i in range(len(prot_records)):
            curr_id = prot_records[i].protein_id
            if (curr_id in ids_to_gbff_and_gbk) and not (curr_id in req_and_found): # This is one of required proteins
                min_index = max(i - max_gene_number, 0)
                max_index = min(i + max_gene_number, len(prot_records) - 1)
                for j in range(min_index, max_index + 1):
                    gbff_to_gbk_to_req_proteins[gbff_prefix][gbk_prefix][prot_records[j].protein_id] = True
                    all_proteins[prot_records[j].protein_id] = prot_records[j]
                gbff_to_gbk_to_ranges[gbff_prefix][gbk_prefix].append([prot_records[min_index].gene_begin, prot_records[max_index].gene_end])
                req_and_found[curr_id] = True

    # 2) Printing information to the database_slice
    n = 0
    for gbff in gbff_to_gbk_to_req_proteins.keys():
        n += 1
        lbar["value"] = int(100 * n / len(gbff_to_gbk_to_req_proteins.keys()))
        root.update()
        gbff_path = os.path.join(output_path, gbff)
        if not os.path.isdir(gbff_path):
            os.mkdir(gbff_path)
        for gbk in gbff_to_gbk_to_req_proteins[gbff].keys():
            merge_overlapping_ranges(gbff_to_gbk_to_ranges[gbff][gbk])
            gbk_path = os.path.join(gbff_path, gbk)
            if not os.path.isdir(gbk_path):
                os.mkdir(gbk_path)
            # 2.1) Writing 'a_' and 'p_' files
            protein_ids = sorted(list(gbff_to_gbk_to_req_proteins[gbff][gbk].keys()), key = lambda k: all_proteins[k].gene_begin, reverse = False)
            write_plain_list(os.path.join(gbk_path, "a_%s.gis" % gbk), protein_ids)
            write_plain_list(os.path.join(gbk_path, "a_%s.ids" % gbk), protein_ids)
            p_fasta_file = open(os.path.join(gbk_path, "p_%s.fasta" % gbk), "w")
            for protein_id in protein_ids:
                p_fasta_file.write("%s" % all_proteins[protein_id].get_name_in_format("ID", "COGNAT"))
                p_fasta_file.write("%s\n\n" % all_proteins[protein_id].sequence)

            # 2.2) Writing 'c_', 'd_' and 'm_' files
            initial_gbk_path = os.path.join(database_path, gbff, gbk)
            filter_table(os.path.join(initial_gbk_path, "c_%s.txt" % gbk), os.path.join(gbk_path, "c_%s.txt" % gbk), protein_ids)
            filter_table(os.path.join(initial_gbk_path, "d_%s.txt" % gbk), os.path.join(gbk_path, "d_%s.txt" % gbk), protein_ids)
            filter_table(os.path.join(initial_gbk_path, "m_%s.txt" % gbk), os.path.join(gbk_path, "m_%s.txt" % gbk), protein_ids)

            # 2.3) Writing 'g_' file
            (nucl_records, found) = udav_fasta.read_fasta(os.path.join(initial_gbk_path, "g_%s.fasta" % gbk), None, None, dict(), 100500, False, None, None, None, None, "COGNAT_nucl")
            g_fasta_file = open(os.path.join(gbk_path, "g_%s.fasta" % gbk), "w") # gbk sequence (in 'COGNAT_nucl' format)
            for [begin, end] in gbff_to_gbk_to_ranges[gbff][gbk]:
                original_name = nucl_records[0].name
                nucl_records[0].name += "|%i|%i" % (begin, len(nucl_records[0].sequence))
                nucl_records[0].print_fasta(g_fasta_file, 60, begin, end)
                nucl_records[0].name = original_name
            g_fasta_file.close()
            #print ("gbff = '%s', gbk = '%s'" % (gbff, gbk))
            #print ("Selected %i proteins in %i sequence ranges" % (len(protein_ids), len(gbff_to_gbk_to_ranges[gbff][gbk])))
    del udav_fasta
    return len(gbff_to_gbk_to_req_proteins.keys())

class DrawableObject:
    def __init__(self, object_type, data, coordinates, fill, outline, non_clickable = False):
        if not object_type in ["line", "arrow", "domain"]:
            raise ValueError("Given <object_type> for a <DrawableObject> is wrong: '%s'" % object_type)
        self.object_type = object_type
        self.data = data
        self.coordinates = coordinates # List of coordinates
        self.fill = fill
        self.outline = outline

    def shift_by_x(self, x):
        """
        Method adds given <x> value to each x coordinate (even in the list of coordinates)
        """
        for i in range(len(self.coordinates)):
            if (i % 2) == 0: # even value, so it is an x
                self.coordinates[i] += x

    def get_best_value(self, coordinate_type, compare_method, compare_with):
        """
        Method returns the 'best' value out of given value <compare_with> and a coordinate of
        specified <coordinate_type> among all coordinates of this object. Comparison is done with
        <compare_method> (e.g. min or max)
        """
        coordinate_type = coordinate_type.lower()
        best_value = compare_with
        if not coordinate_type in ["x", "y"]:
            print ("ERROR: unknow coordinate type '%s' given to <DrawableObject.get_better_value()>, should be either 'x' or 'y'" % coordinate_type)
            return None
        try:
            for i in range(len(self.coordinates)):
                if (coordinate_type == "x") and ((i % 2) == 0): # even value, so it is an x
                    best_value = compare_method(best_value, self.coordinates[i])
                if (coordinate_type == "y") and ((i % 2) != 0): # odd value, so it is an y
                    best_value = compare_method(best_value, self.coordinates[i])
        except:
            print ("ERROR: something went wrong with given comparison_method '%s' given to <DrawableObject.get_better_value()>" % compare_method)
        return best_value

def get_sorted_domains(domains):
    sorted_regions = list()
    try:
        for domain_name in domains:
            domain_string = domains[domain_name].strip()
            regions = domain_string.split(" ")
            for r in regions:
                fields = r.split("..")
                hmm_begin = 0
                hmm_end = 0
                try:
                    hmm_begin = int(fields[4])
                    hmm_end = int(fields[5])
                except ValueError:
                    pass
                values = (int(fields[0]), int(fields[1]), float(fields[2]), float(fields[3]), hmm_begin, hmm_end, domain_name)
                sorted_regions.append(values)
    except IndexError:
        print ("IndexError occured in <get_sorted_domains>")
        print (domains)
        raise
    sorted_regions.sort(key = lambda k: k[0], reverse = False)
    return sorted_regions

def get_domain_rectangles(additional_data, x0, y0, direction, domain_height, gene_length, arrow_tip_start, arrow_tip_size, domains, domain_to_color, base_domain_color):
    #            [0]   [1]    [2]    [3]      [4]      [5]       [6]
    #COG0001 : "begin..end..evalue..score..hmm_begin..hmm_end..hmm_cover <...>"
    domain_objects = list()
    colored_domain_occured = False
    # 1) Getting sorted regions
    sorted_regions = get_sorted_domains(domains)
    # 2) Creating domain rectangles
    for region in sorted_regions:
        domain = region[6]
        data = "%i..%i..%s..%s..%i..%i..%s" % tuple(region)
        data += "..%s" % additional_data
        color = base_domain_color
        if domain in domain_to_color:
            color = domain_to_color[domain]
            colored_domain_occured = True
        domain_start = region[0] * arrow_tip_start / gene_length
        domain_end = region[1] * arrow_tip_start / gene_length
        if direction == 1:
            domain_objects.append(DrawableObject("domain", data, [x0 + domain_start, y0 + 1, x0 + domain_end, y0 + domain_height - 1], color, "", non_clickable = False))
        else:
            domain_objects.append(DrawableObject("domain", data, [x0 + arrow_tip_start + arrow_tip_size - domain_end, y0 + 1, x0 + arrow_tip_start + arrow_tip_size - domain_start, y0 + domain_height - 1], color, "", non_clickable = False))
    return (domain_objects, colored_domain_occured)

def get_domain_architecture(domains):
    domain_architecture = ""
    sorted_regions = get_sorted_domains(domains)
    for region in sorted_regions:
        domain_architecture += "%s|" % region[6]
    domain_architecture = domain_architecture.strip("|")
    return domain_architecture

def get_global_arrow_color(base_color, domains, domain_to_color, base_domain_color):
    """
    Returns <base_color> if no domains are found, <base_domain_color> if no color specified for current
    domain architecture and a color from <domain_to_color> dict (domain names should be connected with spaces)
    """
    arrow_color = base_color
    if domains != None:
        domain_architecture = get_domain_architecture(domains)
        if domain_architecture in domain_to_color:
            arrow_color = domain_to_color[domain_architecture]
        else:
            arrow_color = base_domain_color
    return arrow_color

def get_gene_arrow_and_domains(additional_data, x0, y0, gene_length, protein_length, nucl_per_pixel, arrow_height, arrow_direction, domains, domain_to_color, base_domain_color, color_by_architecture):
    """
    Method will return a list of <DrawableObject> with an arrow representing a gene and
    rectangles representing domains in it.
    <domains> is a dictionary of domain objects.
    If <color_by_architecture> is True, color will be assigned to the whole gene (arrow)
    """
    gold = 0.618
    base_color = "#e6e6e6"
    gene_not_colored = True # Method returns now if this object was colored
    #print ("protein_id = '%s', x0 = %s, y0 = %s, gene_length = %s, nucl_per_pixel = %s, arrow_height = %s, arrow_dir = %s" % (protein_id, x0, y0, gene_length, nucl_per_pixel, arrow_height, arrow_direction))
    objects = list()

    arrow_length = gene_length / nucl_per_pixel
    arrow_tip_size = min((1 - gold) * arrow_length, arrow_height / 2)
    arrow_tip_start = arrow_length - arrow_tip_size
    polygon_coord = list()
    if arrow_direction == 1:
        polygon_coord.append(x0)
        polygon_coord.append(y0 + (arrow_height / 4))
        polygon_coord.append(x0 + arrow_tip_start)
        polygon_coord.append(y0 + (arrow_height / 4))
        polygon_coord.append(x0 + arrow_tip_start)
        polygon_coord.append(y0)
        polygon_coord.append(x0 + arrow_length)
        polygon_coord.append(y0 + (arrow_height / 2))
        polygon_coord.append(x0 + arrow_tip_start)
        polygon_coord.append(y0 + arrow_height)
        polygon_coord.append(x0 + arrow_tip_start)
        polygon_coord.append(y0 + (3 * arrow_height / 4))
        polygon_coord.append(x0)
        polygon_coord.append(y0 + (3 * arrow_height / 4))
    else:
        polygon_coord.append(x0 + arrow_length)
        polygon_coord.append(y0 + (arrow_height / 4))
        polygon_coord.append(x0 + arrow_length - arrow_tip_start)
        polygon_coord.append(y0 + (arrow_height / 4))
        polygon_coord.append(x0 + arrow_length - arrow_tip_start)
        polygon_coord.append(y0)
        polygon_coord.append(x0)
        polygon_coord.append(y0 + (arrow_height / 2))
        polygon_coord.append(x0 + arrow_length - arrow_tip_start)
        polygon_coord.append(y0 + arrow_height)
        polygon_coord.append(x0 + arrow_length - arrow_tip_start)
        polygon_coord.append(y0 + (3 * arrow_height / 4))
        polygon_coord.append(x0 + arrow_length)
        polygon_coord.append(y0 + (3 * arrow_height / 4))
    if color_by_architecture == True:
        new_color = get_global_arrow_color(base_color, domains, domain_to_color, base_domain_color)
        if (new_color != base_color) and (new_color != base_domain_color): # This gene is colored somehow
            gene_not_colored = False
    arrow_object = DrawableObject("arrow", additional_data, polygon_coord, base_color, "#000000", non_clickable = False)
    objects.append(arrow_object)
    if (domains != None) and (color_by_architecture == False):
        (domain_objects, colored_domain_occured) = get_domain_rectangles(additional_data, x0, y0 + (arrow_height / 4), arrow_direction, arrow_height / 2, protein_length, arrow_tip_start, arrow_tip_size, domains, domain_to_color, base_domain_color)
        if colored_domain_occured:
            gene_not_colored = False
        objects.extend(domain_objects)
    return (objects, gene_not_colored)

def report_nucl_seq_to_widget(widget, data):
    #           [0]       [1]         [2]      [3]       [4]         [5]       [6]
    # data: "sequence..real_begin..real_end..record..was_reversed..gene1_id..gene2_id"
    [sequence, real_begin, real_end, record, was_reversed, gene1_id, gene2_id] = data.split("..")
    if was_reversed == "True":
        widget.add_string(">%s_complement(%i-%i)\n" % (record, int(real_begin), int(real_end)))
    else:
        widget.add_string(">%s_%i-%i\n" % (record, int(real_begin), int(real_end)))
    widget.add_string("%s\n\n" % sequence)

def report_protein_seq_to_widget(widget, p, curr_org):
    widget.add_string(">%s %s [%s]\n" % (p.protein_id.strip("*"), p.product, curr_org))
    widget.add_string("%s\n\n" % p.sequence)

def report_gene_to_widget(widget, p, real_begin, real_end, real_direction, nucl_seq, curr_org, curr_gbff, curr_taxonomy, domains, domain_to_name):
    widget.text_widget.tag_configure("header", background = "#BBBBBB")
    widget.add_tagged_string("Protein gene was selected:", "header")
    widget.add_string("organism   : %s\n" % curr_org)
    widget.add_string("taxonomy   : %s\n" % curr_taxonomy)
    widget.add_string("assembly   : %s\n" % curr_gbff)
    widget.add_string("record     : %s\n\n" % p.source_record)

    widget.add_string("protein_id : %s\n" % p.protein_id.strip("*"))
    widget.add_string("product    : %s\n" % p.product)
    widget.add_string("length (aa): %s\n" % len(p.sequence))
    widget.add_string("gene_begin : %s\n" % real_begin)
    widget.add_string("gene_end   : %s\n" % real_end)
    widget.add_string("direction  : %s\n" % real_direction)
    widget.add_string("coordinates: %s..%s\n\n" % (real_begin, real_end))

    widget.add_tagged_string("The following domains were detected:", "header")
    if domains == None:
        widget.add_string("No domains detected under current parameters\n")
    else:
        sorted_regions = list()
        try:
            for domain_name in domains:
                domain_string = domains[domain_name].strip()
                regions = domain_string.split(" ")
                for r in regions:
                    fields = r.split("..")
                    domain_full_name = "Unknown domain name"
                    if domain_name in domain_to_name:
                        domain_full_name = domain_to_name[domain_name]
                    values = (int(fields[0]), int(fields[1]), float(fields[2]), float(fields[3]), fields[4], fields[5], domain_name, domain_full_name)
                    sorted_regions.append(values)
        except IndexError:
            print ("IndexError occured in <report_gene_to_widget>")
            raise
        sorted_regions.sort(key = lambda k: int(k[0]), reverse = False)
        widget.add_string("Domain\tbegin\tend\tscore\tevalue\thmm_b\thmm_e\n")
        for region in sorted_regions:
            # [0]   [1]    [2]    [3]     [4]        [5]     [6]         [7]
            #begin..end..evalue..score..hmm_begin..hmm_end..domain..domain_full_name
            widget.add_string("%s\t%i\t%i\t%.1f\t%s\t%s\t%s\t%s\n" % (region[6], region[0], region[1], region[3], region[2], region[4], region[5], region[7]))
    widget.add_string("\n")
    widget.add_tagged_string("Protein sequence (FASTA):", "header")
    report_protein_seq_to_widget(widget, p, curr_org)

    widget.add_tagged_string("Gene sequence (FASTA, 5'->3'):", "header")
    if real_direction == 1:
        widget.add_string(">%s:%i-%i\n" % (p.source_record, real_begin, real_end))
    else:
        widget.add_string(">%s:complement(%i-%i)\n" % (p.source_record, real_begin, real_end))
    widget.add_string("%s\n" % nucl_seq)

def report_domain_to_widget(widget, p, data, domain_to_name):
    # data: "begin..end..evalue..score..domain..protein_id"
    widget.text_widget.tag_configure("header", background = "#BBBBBB")
    widget.add_tagged_string("Protein domain was selected for the following protein:", "header")
    widget.add_string("protein_id : %s\n" % p.protein_id.strip("*"))
    widget.add_string("product    : %s\n" % p.product)
    widget.add_string("length (aa): %s\n\n" % len(p.sequence))

    widget.add_tagged_string("Current domain data:", "header")
    widget.add_string("Domain\tbegin\tend\tscore\tevalue\thmm_b\thmm_e\tname\n")
    fields = data.split("..")
    begin = int(fields[0])
    end = int(fields[1])
    domain_name = fields[6]
    domain_full_name = "Unknown domain name"
    if domain_name in domain_to_name:
        domain_full_name = domain_to_name[domain_name]
    widget.add_string("%s\t%i\t%i\t%.1f\t%s\t%s\t%s\t%s\n\n" % (domain_name, begin, end, float(fields[3]), float(fields[2]), fields[4], fields[5], domain_full_name))
    widget.add_tagged_string("Current region sequence:", "header")
    widget.add_string(">%s_%i-%i\n" % (p.protein_id.strip("*"), begin, end))
    widget.add_string("%s\n" % p.sequence[begin - 1:end + 1])

def report_intergene_to_widget(widget, data):
    #           [0]       [1]         [2]      [3]       [4]         [5]       [6]
    # data: "sequence..real_begin..real_end..record..was_reversed..gene1_id..gene2_id"
    widget.text_widget.tag_configure("header", background = "#BBBBBB")
    widget.add_tagged_string("Intergene region was selected:", "header")
    [sequence, real_begin, real_end, record, was_reversed, gene1_id, gene2_id] = data.split("..")

    widget.add_string("record     : %s\n" % record)
    widget.add_string("begin      : %s\n" % real_begin)
    widget.add_string("end        : %s\n" % real_end)
    widget.add_string("coordinates: %s..%s\n" % (real_begin, real_end))
    widget.add_string("reversed   : %s\n\n" % was_reversed)
    widget.add_string("gene1_pid  : %s\n" % gene1_id)
    widget.add_string("gene2_pid  : %s\n\n" % gene2_id)

    widget.add_tagged_string("Current region sequence:", "header")
    
    if was_reversed == "True":
        widget.add_string(">%s_complement(%i-%i)\n" % (record, int(real_begin), int(real_end)))
    else:
        widget.add_string(">%s_%i-%i\n" % (record, int(real_begin), int(real_end)))
    widget.add_string("%s\n" % sequence)

def get_complement(sequence):
    """
    Method reverses the nucleic acid <sequence> by getting its complement and making direct reverse
    """
    complement = {"A":"T", "a":"t", "T":"A", "t":"a", "G":"C", "g":"c", "C":"G", "c":"g", "n":"n", "N":"N", "x":"x", "X":"X"}
    rev = ""
    for s in sequence:
        if s in complement:
            rev += complement[s]
        else:
            print ("FATAL ERROR: unknown nucleotide '%s' was detected in <get_complement> method!" % s)
            return None
    rev = rev[::-1] #Extended slice syntax: reversing order
    return rev

def get_real_coordinates_and_sequence(begin, end, was_reversed, nucl_data_list):
    """
    Method will return a tuple of begin and end in real sequence, as well as the sequence of this
    part and id of nucleic record where is was identified
    """
    begin_in_seq = begin # These values are used for the sequence taking
    end_in_seq = end
    sequence = ""
    record_id = ""
    for nucl_data in nucl_data_list:
        if nucl_data.gene_begin == -1: # This is complete record
            if was_reversed:
                L = len(nucl_data.sequence)
                begin_in_seq = L - begin_in_seq + 1
                end_in_seq = L - end_in_seq + 1
                begin_in_seq, end_in_seq = end_in_seq, begin_in_seq
            sequence = nucl_data.sequence[begin_in_seq : end_in_seq + 1]
            record_id = nucl_data.gi
            if len(nucl_data_list) != 1:
                print ("WARNING: more than one nucleic record found, but no information about its start known: '%s'" % nucl_data.gi)
                break
        else:
            nucl_begin = nucl_data.gene_begin # an artificial information about real begin of this region is stored in the <gene_begin>
            nucl_end = nucl_begin + len(nucl_data.sequence) - 1
            if was_reversed:
                L = nucl_data.gene_direction # an artificial information about current record length
                begin_in_seq = L - begin_in_seq + 1
                end_in_seq = L - end_in_seq + 1
                begin_in_seq, end_in_seq = end_in_seq, begin_in_seq
            if (begin_in_seq >= nucl_begin) and (end_in_seq <= nucl_end):
                sequence = nucl_data.sequence[begin_in_seq - nucl_begin : end_in_seq + 1 - nucl_begin]
                record_id = nucl_data.gi
                break
    if sequence == "": # No sequence found
        sequence = "XXX"
    return (begin_in_seq, end_in_seq, sequence, record_id)

def read_taxonomy_colors(input_filename, hex_mode = False):
    """
    Method reads a tab-separated assignment between gbff identifier, color (in hex-code mode, if
    <hex_mode> is True) and additional taxonomy details.
    """
    gbff_colors = dict()
    gbff_order = list()
    input_file = open(input_filename)
    for string in input_file:
        #GCA_000246735.1_ASM24673v1	153	0	0	Archaea	Candidatus Thermoplasmatota
        string = string.strip()
        if len(string) == 0:
            continue
        if string[0] == "#":
            continue
        fields = string.split("\t")
        gbff = fields.pop(0)
        color = None
        if hex_mode:
            color = fields.pop(0)
        else:
            red = str( hex(int(fields.pop(0))) ).replace("0x", "")
            if len(red) == 1:
               red = "0" + red
            green = str( hex(int(fields.pop(0))) ).replace("0x", "")
            if len(green) == 1:
               green = "0" + green
            blue = str( hex(int(fields.pop(0))) ).replace("0x", "")
            if len(blue) == 1:
               blue = "0" + blue
            color = "#" + red + green + blue
        gbff_colors[gbff] = color
        gbff_order.append([gbff, "|".join(fields)])
    input_file.close()
    return (gbff_colors, gbff_order)

def print_simple_legend(name_to_color, order_list, svg_filename):
    tax_types = dict()
    names_and_colors = list()
    for i in range(len(order_list)):
        curr_gbff = order_list[i][0]
        curr_name = order_list[i][1]
        if not curr_name in tax_types: # New type of taxonomy
            curr_color = name_to_color[curr_gbff]
            names_and_colors.append([curr_name, curr_color])
            tax_types[curr_name] = True

    import drawsvg
    letter_w = 7.195 * 2 #7.2
    letter_h = 9.6 * 2
    max_name = 0
    for gbff, name in order_list:
        if len(name) > max_name:
            max_name = len(name)
    legend_size_x = letter_w * 5
    legend_size_y = letter_h
    spacer = 20
    field = 5
    max_x_size = legend_size_x + spacer + (max_name * letter_w) + (field * 2)
    max_y_size = (legend_size_y * len(names_and_colors)) + (field * 2)

    drawing = drawsvg.Drawing(max_x_size, max_y_size, origin = (0, 0))
    for i in range(len(names_and_colors)):
        curr_name = names_and_colors[i][0]
        curr_color = names_and_colors[i][1]
        drawing.append(drawsvg.Rectangle(field, (i * letter_h) + field, legend_size_x, legend_size_y, fill = curr_color))
        drawing.append(drawsvg.Text(curr_name, letter_h, x = field + legend_size_x + spacer, y = (i * letter_h) + field + 2, fill = curr_color, text_anchor = "start", dominant_baseline = "hanging", font_family = "Arial"))
    drawing.save_svg(svg_filename)
    del drawsvg

def proceed_gene_not_colored_list(gene_not_colored_list, to_remove):
    for pair in gene_not_colored_list:
        # pair = (protein_id, True/False)
        protein_id = pair[0]
        not_colored = pair[1]
        if not_colored:
            to_remove[protein_id] = True
        else:
            break

def hide_uncolored(objects_to_draw, gene_not_colored_list):
    """
    Method finds genes which flank the neighborhood and are not colored and
    filters <objects_to_draw> in order to remove them.
    """
    to_remove = dict()
    proceed_gene_not_colored_list(gene_not_colored_list, to_remove)
    gene_not_colored_list.reverse()
    proceed_gene_not_colored_list(gene_not_colored_list, to_remove)
    i = 0
    while i < len(objects_to_draw):
        data_parts = objects_to_draw[i].data.split("..")
        if objects_to_draw[i].object_type in ["arrow", "domain"]:
            protein_id = data_parts[-1]
            if protein_id in to_remove:
                objects_to_draw.pop(i)
                i -= 1
        elif objects_to_draw[i].object_type == "line":
            this_protein_id = data_parts[-1]
            next_protein_id = data_parts[-2]
            if (this_protein_id in to_remove) or (next_protein_id in to_remove):
                objects_to_draw.pop(i)
                i -= 1
        i += 1
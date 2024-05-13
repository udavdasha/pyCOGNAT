# -*- coding: utf-8 -*-
import os, re, time
import tkinter
import tkinter.ttk
import pyCOGNAT_basic

class ProjectData(tkinter.Frame):
    def __init__(self, parent, host):
        tkinter.Frame.__init__(self, parent)
        self.host = host
        self.p = self.host.p
        self.project_name_widget = None    # Entry for the project name      |
        self.max_gene_number = None        # Entry for num. of genes to save |
        self.id_filename_widget = None     # Entry for the id filename       | project data frame
        self.database_widget = None        # Entry for the database filename |
        self.initital_progress = None      # Progressbar                     |
        self.reading_progress = None       # Progressbar                     |
        self.loading_progress = None       # Progressbar                     |
        self.add_clade_button = None       # Button to add new clade         |

        self.evalue_widget = None          # Entry for the e-value threshold       |
        self.overlap_widget = None         # Entry for the domain overlap widget   |
        self.unite_check = None            # Checkbutton for the unite domains     |
        self.unite_domains = None          # tkinter boolean variable              | domain search options
        self.max_distance_widget = None    # Entry for the max. distance widget    |
        self.max_hmm_overlap_widget = None # Entry for the max. HMM overlap widget |

        self.domain_colors = None

        self.protein_ids = None

        self.man_log = None
        self.console = None

        self.create_UI()

    def create_UI(self):
        self.grid_rowconfigure(0, weight = 1)
        self.grid_columnconfigure(0, weight = 1)

        central_panel = tkinter.PanedWindow(self, orient = tkinter.HORIZONTAL, sashwidth = self.p * 2, sashrelief = tkinter.RIDGE, background = self.host.back)
        central_panel.grid(row = 0, column = 0, sticky = "NSEW")

        left_part = tkinter.Frame(central_panel)
        left_part.grid_rowconfigure(2, weight = 1)
        left_part.grid_columnconfigure(0, weight = 1)

        data_frame = tkinter.Frame(left_part)
        data_frame.grid_columnconfigure(1, weight = 1)
        tkinter.Label(data_frame, text = "Project data:", font = ("Arial", 10, "bold"), justify = tkinter.LEFT).grid(row = 0, column = 0, sticky = "NSW", padx = self.p, pady = self.p)
        tkinter.Label(data_frame, text = "Project name:", justify = tkinter.LEFT).grid(row = 1, column = 0, sticky = "NSW", padx = self.p, pady = self.p)
        self.project_name_widget = tkinter.Entry(data_frame, width = 20, font = ("Courier New", 12, "bold"), background = self.host.header, foreground = "#FFFFFF")
        self.project_name_widget.grid(row = 1, column = 1, sticky = "NSEW", padx = self.p, pady = self.p)
        tkinter.Label(data_frame, text = "Max gene num.:", justify = tkinter.LEFT).grid(row = 2, column = 0, sticky = "NSW", padx = self.p, pady = self.p)
        self.max_gene_number = tkinter.Entry(data_frame, width = 20, font = ("Courier New", 10))
        self.max_gene_number.insert(tkinter.END, "30")
        self.max_gene_number.grid(row = 2, column = 1, sticky = "NSEW", padx = self.p, pady = self.p)
        tkinter.Label(data_frame, text = "ID filename:", justify = tkinter.LEFT).grid(row = 3, column = 0, sticky = "NSW", padx = self.p, pady = self.p)
        self.id_filename_widget = tkinter.Entry(data_frame, width = 20, font = ("Courier New", 10), state = tkinter.DISABLED)
        self.id_filename_widget.grid(row = 3, column = 1, sticky = "NSEW", padx = self.p, pady = self.p)
        tkinter.Label(data_frame, text = "Database:", justify = tkinter.LEFT).grid(row = 4, column = 0, sticky = "NSW", padx = self.p, pady = self.p)
        self.database_widget = tkinter.Entry(data_frame, width = 20, font = ("Courier New", 10), state = tkinter.DISABLED)
        self.database_widget.grid(row = 4, column = 1, sticky = "NSEW", padx = self.p, pady = self.p)
        tkinter.Label(data_frame, text = "Initialising:", justify = tkinter.LEFT).grid(row = 5, column = 0, sticky = "NSW", padx = self.p, pady = self.p)
        self.initital_progress = tkinter.ttk.Progressbar(data_frame, value = 0, style = "success.Striped.Horizontal.TProgressbar")
        self.initital_progress.grid(row = 5, column = 1, sticky = "NSEW", padx = self.p, pady = self.p)
        tkinter.Label(data_frame, text = "Reading:", justify = tkinter.LEFT).grid(row = 6, column = 0, sticky = "NSW", padx = self.p, pady = self.p)
        self.reading_progress = tkinter.ttk.Progressbar(data_frame, value = 0, style = "success.Striped.Horizontal.TProgressbar")
        self.reading_progress.grid(row = 6, column = 1, sticky = "NSEW", padx = self.p, pady = self.p)
        tkinter.Label(data_frame, text = "Loading:", justify = tkinter.LEFT).grid(row = 7, column = 0, sticky = "NSW", padx = self.p, pady = self.p)
        self.loading_progress = tkinter.ttk.Progressbar(data_frame, value = 0, style = "success.Striped.Horizontal.TProgressbar")
        self.loading_progress.grid(row = 7, column = 1, sticky = "NSEW", padx = self.p, pady = self.p)
        self.add_clade_button = tkinter.Button(data_frame, text = "Add clade", font = ("Arial", 10, "bold"), command = self.host.add_clade_tab, state = tkinter.DISABLED)
        self.add_clade_button.grid(row = 8, column = 0, sticky = "NSEW", padx = self.p, pady = self.p)
        data_frame.grid(row = 0, column = 0, sticky = "NSEW")

        domain_frame = tkinter.Frame(left_part, highlightbackground = "#888888", highlightthickness = 1)
        domain_frame.grid_columnconfigure(1, weight = 1)
        tkinter.Label(domain_frame, text = "Domain search options:", font = ("Arial", 10, "bold"), justify = tkinter.LEFT).grid(row = 0, column = 0, sticky = "NSW", padx = self.p, pady = self.p)
        tkinter.Label(domain_frame, text = "E-value threshold:", justify = tkinter.LEFT).grid(row = 1, column = 0, sticky = "NSW", padx = self.p, pady = self.p)
        self.evalue_widget = tkinter.Entry(domain_frame, width = 10, font = ("Courier New", 10), state = tkinter.DISABLED)
        self.evalue_widget.grid(row = 1, column = 1, sticky = "NSEW", padx = self.p, pady = self.p)
        tkinter.Label(domain_frame, text = "Overlap threshold (%):", justify = tkinter.LEFT).grid(row = 2, column = 0, sticky = "NSW", padx = self.p, pady = self.p)
        self.overlap_widget = tkinter.Entry(domain_frame, width = 10, font = ("Courier New", 10), state = tkinter.DISABLED)
        self.overlap_widget.grid(row = 2, column = 1, sticky = "NSEW", padx = self.p, pady = self.p)
        self.unite_domains = tkinter.BooleanVar()
        self.unite_domains.set(False)
        self.unite_check = tkinter.Checkbutton(domain_frame, text = "Unite domains", variable = self.unite_domains, command = self.configure_unite, state = tkinter.DISABLED)
        self.unite_check.grid(row = 3, column = 0, sticky ="NSW", padx = self.p, pady = self.p)
        tkinter.Label(domain_frame, text = "Max. distance (aa):", justify = tkinter.LEFT).grid(row = 4, column = 0, sticky = "NSW", padx = self.p, pady = self.p)
        self.max_distance_widget = tkinter.Entry(domain_frame, width = 10, font = ("Courier New", 10), state = tkinter.DISABLED)
        self.max_distance_widget.grid(row = 4, column = 1, sticky = "NSEW", padx = self.p, pady = self.p)
        tkinter.Label(domain_frame, text = "Max. HMM overlap (%):", justify = tkinter.LEFT).grid(row = 5, column = 0, sticky = "NSW", padx = self.p, pady = self.p)
        self.max_hmm_overlap_widget = tkinter.Entry(domain_frame, width = 10, font = ("Courier New", 10), state = tkinter.DISABLED)
        self.max_hmm_overlap_widget.grid(row = 5, column = 1, sticky = "NSEW", padx = self.p, pady = self.p)
        domain_frame.grid(row = 1, column = 0, sticky = "NSEW")

        self.domain_colors = pyCOGNAT_basic.TextFrameWithLabel(left_part, self.p, self.host.back, "#000000", "Domain color code:", ("Arial", 10, "bold"), None)
        self.domain_colors.grid(row = 2, column = 0, sticky = "NSEW")
        self.domain_colors.disable_text()
        central_panel.add(left_part)

        self.protein_ids = pyCOGNAT_basic.TextFrameWithLabel(central_panel, self.p, self.host.back, "#000000", "List of protein IDs:", ("Arial", 10, "bold"))
        central_panel.add(self.protein_ids)

        right_part = tkinter.Frame(central_panel)
        right_part.grid_rowconfigure(0, weight = 1)
        right_part.grid_columnconfigure(0, weight = 1)
        self.man_log = pyCOGNAT_basic.TextFrameWithLabel(right_part, self.p, self.host.back, "#000000", "Manual log:", ("Arial", 10, "bold"))
        self.man_log.grid(row = 0, column = 0, sticky = "NSEW")
        self.console = pyCOGNAT_basic.TextFrameWithLabel(right_part, self.p, self.host.back, "#000000", "Console:", ("Arial", 10, "bold"))
        self.console.grid(row = 1, column = 0, sticky = "NSEW")
        central_panel.add(right_part)

        self.update_idletasks()
        central_panel.sash_place(0, 335, 1)
        central_panel.sash_place(1, 670, 1)

    def get_project_name(self):
        return self.project_name_widget.get().strip()

    def get_max_gene_number(self):
        return int(self.max_gene_number.get().strip())

    def get_evalue(self):
        return float(self.evalue_widget.get().strip())

    def get_overlap(self):
        return float(self.overlap_widget.get().strip())

    def get_max_distance(self):
        return int(self.max_distance_widget.get().strip())

    def get_max_hmm_overlap(self):
        return float(self.max_hmm_overlap_widget.get().strip())

    def get_domain_to_color(self):
        domain_to_color = dict()
        strings = self.domain_colors.get_content_as_list()
        for string in strings:
            fields = string.split("\t")
            domain_to_color[fields[0]] = fields[1]
        return domain_to_color

    def print_to_console(self, text):
        self.console.text_widget.insert(tkinter.END, "\n\n")
        t = time.localtime()
        time_string = "%02i.%02i.%i\t%02i:%02i:%02i" % (t.tm_mday, t.tm_mon, t.tm_year, t.tm_hour, t.tm_min, t.tm_sec)
        self.console.text_widget.insert(tkinter.END, "%s\n" % time_string)
        self.console.text_widget.insert(tkinter.END, "%s\n" % text)

    def configure_unite(self):
        if self.unite_domains.get() == False:
            self.max_distance_widget.configure(state = tkinter.DISABLED)
            self.max_hmm_overlap_widget.configure(state = tkinter.DISABLED)
        else:
            self.max_distance_widget.configure(state = tkinter.NORMAL)
            self.max_hmm_overlap_widget.configure(state = tkinter.NORMAL)

    def clear(self):
        print ("---- Project clearing ----")
        pyCOGNAT_basic.clear_and_disable_entry_widget(self.project_name_widget, do_not_disable = True)
        pyCOGNAT_basic.clear_and_disable_entry_widget(self.max_gene_number, do_not_disable = True)
        pyCOGNAT_basic.clear_and_disable_entry_widget(self.id_filename_widget)
        pyCOGNAT_basic.clear_and_disable_entry_widget(self.database_widget)

        pyCOGNAT_basic.clear_and_disable_entry_widget(self.evalue_widget)
        pyCOGNAT_basic.clear_and_disable_entry_widget(self.overlap_widget)
        pyCOGNAT_basic.clear_and_disable_entry_widget(self.max_distance_widget)
        pyCOGNAT_basic.clear_and_disable_entry_widget(self.max_hmm_overlap_widget)

        self.domain_colors.clear(disable = True)
        self.protein_ids.clear(disable = True)
        self.man_log.clear(disable = False)
        self.console.clear(disable = False)

    def enable_analysis(self, add_defaults = True):
        self.domain_colors.enable_text()
        self.evalue_widget.configure(state = tkinter.NORMAL)
        self.overlap_widget.configure(state = tkinter.NORMAL)
        self.unite_check.configure(state = tkinter.NORMAL)
        self.max_distance_widget.configure(state = tkinter.NORMAL)
        self.max_hmm_overlap_widget.configure(state = tkinter.NORMAL)
        self.add_clade_button.configure(state = tkinter.NORMAL, background = self.host.header, foreground = "#FFFFFF")
        if add_defaults:
            self.evalue_widget.insert(tkinter.END, "1e-5")
            self.overlap_widget.insert(tkinter.END, "5")
            self.max_distance_widget.insert(tkinter.END, "50")
            self.max_hmm_overlap_widget.insert(tkinter.END, "30")

    def fix_status(self, ids_path = None, database_path = None):
        self.project_name_widget.configure(state = tkinter.DISABLED)
        self.max_gene_number.configure(state = tkinter.DISABLED)

        if ids_path != None:
            self.id_filename_widget.configure(state = tkinter.NORMAL)
            self.id_filename_widget.insert(tkinter.END, ids_path)
        self.id_filename_widget.configure(state = tkinter.DISABLED)

        if database_path != None:
            self.database_widget.configure(state = tkinter.NORMAL)
            self.database_widget.insert(tkinter.END, database_path)
        self.database_widget.configure(state = tkinter.DISABLED)

        self.protein_ids.disable_text()

    def get_filename(self, file_type, external = None):
        prefix = self.get_project_name()
        if external != None:
            prefix = external
        type_to_filename = {"params"        : "%s_params.txt" % prefix,
                            "domain_colors" : "%s_domain_colors.txt" % prefix,
                            "protein_ids"   : "%s_protein_ids.txt" % prefix,
                            "man_log"       : "%s_man_log.txt" % prefix,
                            "console"       : "%s_console.txt" % prefix}
        if not file_type in type_to_filename.keys():
            print ("Unknown file type '%s', cannot get proper name" % file_type)
            return None
        result = type_to_filename[file_type]
        if external == None:
            work_dir = self.host.settings.work_dir
            result = os.path.join(work_dir, prefix, result)
        return result

    def get_widget_name_to_widget(self):
        widget_name_to_widget = {"project_name_widget"    : self.project_name_widget,
                                 "max_gene_number"        : self.max_gene_number,
                                 "id_filename_widget"     : self.id_filename_widget,
                                 "database_widget"        : self.database_widget,
                                 "evalue_widget"          : self.evalue_widget,
                                 "overlap_widget"         : self.overlap_widget,
                                 "max_distance_widget"    : self.max_distance_widget,
                                 "max_hmm_overlap_widget" : self.max_hmm_overlap_widget,
                                 "unite_domains"          : self.unite_domains}
        return widget_name_to_widget

    def save(self):
        params_filename = self.get_filename("params")
        widget_name_to_widget = self.get_widget_name_to_widget()
        params_file = open(params_filename, "w")
        for wname in widget_name_to_widget.keys():
            wvalue = widget_name_to_widget[wname].get()
            if wname != "unite_domains":
                wvalue = widget_name_to_widget[wname].get().strip()
            params_file.write("%s\t%s\n" % (wname, wvalue))
        params_file.close()
        self.domain_colors.write_into_file(self.get_filename("domain_colors"))
        self.protein_ids.write_into_file(self.get_filename("protein_ids"))
        self.man_log.write_into_file(self.get_filename("man_log"), use_codecs = True)
        self.console.write_into_file(self.get_filename("console"))

    def load(self, project_dir):
        print ("---- Project loading ----")
        external_project_name = os.path.basename(project_dir)
        params_filename = os.path.join(project_dir, self.get_filename("params", external = external_project_name))
        widget_name_to_widget = self.get_widget_name_to_widget()
        params = pyCOGNAT_basic.read_plain_list(params_filename)
        for element in params:
            fields = element.split("\t")
            if len(fields) == 1: # No data should be loaded to this widget:
                continue
            if fields[0] in widget_name_to_widget:
                wname = fields[0]
                wvalue = fields[1]
                if wname != "unite_domains":
                    widget_name_to_widget[wname].configure(state = tkinter.NORMAL)
                    widget_name_to_widget[wname].insert(tkinter.END, wvalue)
                else:
                    if wvalue == "True":
                        widget_name_to_widget[wname].set(True)
                    else:
                        widget_name_to_widget[wname].set(False)
            else:
                print ("Unknown parameter '%s' with value '%s' in parameters file, skipping..." % (fields[0], fields[1]))

        self.domain_colors.read_from_file(os.path.join(project_dir, self.get_filename("domain_colors", external = external_project_name)))
        self.protein_ids.read_from_file(os.path.join(project_dir,self.get_filename("protein_ids", external = external_project_name)))
        self.man_log.read_from_file(os.path.join(project_dir,self.get_filename("man_log", external = external_project_name)), use_codecs = True)
        self.console.read_from_file(os.path.join(project_dir,self.get_filename("console", external = external_project_name)))
        self.fix_status()
        self.initital_progress["value"] = 100
        self.loading_progress["value"] = 100
        self.reading_progress["value"] = 100
        self.enable_analysis(add_defaults = False)

# -*- coding: utf-8 -*-
import sys, os, re, time, platform, random
import tkinter
import tkinter.ttk
import tkinter.font
import tkinter.colorchooser
import pyCOGNAT_basic

def get_canvas(parent, width, height):
    h = tkinter.ttk.Scrollbar(parent, orient = tkinter.HORIZONTAL)
    v = tkinter.ttk.Scrollbar(parent, orient = tkinter.VERTICAL)
    canvas = tkinter.Canvas(parent, scrollregion=(0, 0, width, height), bg = "#ffffff", yscrollcommand = v.set, xscrollcommand = h.set)
    h["command"] = canvas.xview
    v["command"] = canvas.yview

    canvas.grid(column = 0, row = 0, sticky = "NWES")
    h.grid(column = 0, row = 1, sticky = "WE")
    v.grid(column = 1, row = 0, sticky = "NS")
    parent.grid_columnconfigure(0, weight = 1)
    parent.grid_rowconfigure(0, weight = 1)
    return canvas

class CladeFrame(tkinter.Frame):
    def __init__(self, parent, host):
        tkinter.Frame.__init__(self, parent)
        self.host = host
        self.p = self.host.p
        self.clade_name_widget = None
        self.clade_ids = None
        self.neighbor_number_widget = None
        self.neighbor_progress = None # progress bar for neighborhood search
        self.sort_progress = None # progress bar for domain obtaining
        self.add_colors = None
        self.report_widget = None
        self.main_canvas = None
        self.main_canvas_frame = None
        self.action_mode = None # Frame for buttons which specify click event reaction
        self.min_occur = None # Entry for minimal gene occur
        self.legend_canvas = None
        self.domain_type = None # tkinter string variable
        self.hide_uncolored = None # tkinter boolean variable
        self.color_type = None # tkinter string variable
        self.show_org = None # tkinter boolean variable
        self.export_legend = None # tkinter boolean variable
        self.x_mult_widget = None
        self.y_mult_widget = None
        self.offset_widget = None

        self.neighborhoods = None
        self.id_to_domains = None
        self.proteins = None # Current dictionary of proteins (<Annotated_sequence> objects)
        self.protein_ids_safe = None
        self.nucl_id_to_data = None # Current assignment of intergene regions to their data
        self.protein_id_to_nucl_data = None
        self.domain_to_color = None
        self.domain_occur = None # Occurence of domain (as in legend)
        self.drawable_objects = None # Current list of <pyCOGNAT_basic.DrawableObject> objects  |
        self.text_objects = None                                                              # | Main canvas objects,
        self.max_x = None                                                                     # | size and text data
        self.max_y = None                                                                     # |

        self.legend_objects = None                                                            # | Legend canvas objects
        self.legend_text_objects = None                                                       # |

        self.base_domain_color = "#888888"

        self.selected = None

        self.create_UI()
        self.main_canvas.bind("<Button-1>", self.canvas_clicked)
        self.main_canvas.bind("<Double-Button-1>", self.canvas_doubleclicked)
        self.main_canvas.bind("<MouseWheel>", self.on_mousewheel)
        self.legend_canvas.bind("<Double-Button-1>", self.legend_doubleclicked)

    def create_UI(self):
        self.grid_rowconfigure(0, weight = 1)
        self.grid_columnconfigure(0, weight = 1)

        central_panel = tkinter.PanedWindow(self, orient = tkinter.HORIZONTAL, sashwidth = self.p * 2, sashrelief = tkinter.RIDGE, background = self.host.back)
        central_panel.grid(row = 0, column = 0, sticky = "NSEW")

        left_part = tkinter.Frame(central_panel)
        left_part.grid_rowconfigure(2, weight = 1)
        left_part.grid_columnconfigure(0, weight = 1)
        tkinter.Label(left_part, text = "Clade name:", font = ("Arial", 10, "bold"), justify = tkinter.LEFT).grid(row = 0, column = 0, columnspan = 2, sticky = "NSW", padx = self.p, pady = self.p)
        self.clade_name_widget = tkinter.Entry(left_part, width = 20, font = ("Courier New", 10), background = self.host.header, foreground = "#FFFFFF")
        self.clade_name_widget.grid(row = 1, column = 0, columnspan = 2, sticky = "NSEW", padx = self.p, pady = self.p)
        self.clade_ids = pyCOGNAT_basic.TextFrameWithLabel(left_part, self.p, self.host.back, "#000000", "List of protein IDs:", ("Arial", 10, "bold"), None)
        self.clade_ids.grid(row = 2, column = 0, columnspan = 2, sticky = "NSEW")
        tkinter.Label(left_part, text = "Neighbors num.:", justify = tkinter.LEFT).grid(row = 3, column = 0, sticky = "NSW", padx = self.p, pady = self.p)
        self.neighbor_number_widget = tkinter.Entry(left_part, width = 6, font = ("Courier New", 10))
        self.neighbor_number_widget.insert(tkinter.END, "2")
        self.neighbor_number_widget.grid(row = 3, column = 1, sticky = "NSEW", padx = self.p, pady = self.p)
        self.domain_type = tkinter.StringVar()
        self.domain_type.set("COG")
        tkinter.Label(left_part, text = "Domain type:", justify = tkinter.LEFT).grid(row = 4, column = 0, sticky = "NSW", padx = self.p, pady = self.p)
        domain_type_choice = tkinter.ttk.Combobox(left_part, width = 7, state = "readonly", textvariable = self.domain_type)
        domain_type_choice["values"] = ("COG", "Pfam")
        domain_type_choice.grid(row = 4, column = 1, sticky = "NSW", padx = self.p, pady = self.p)
        domain_type_choice.current(0)
        self.hide_uncolored = tkinter.BooleanVar()
        self.hide_uncolored.set(False)
        h = tkinter.Checkbutton(left_part, text = "Hide non-colored edges", variable = self.hide_uncolored)
        h.grid(row = 5, column = 0, columnspan = 2, sticky ="NSW", padx = self.p, pady = self.p)
        #domain_type_choice.bind("<<ComboboxSelected>>", self.draw_data)
        tkinter.Button(left_part, text = "Draw", font = ("Arial", 10, "bold"), background = self.host.header, foreground = "#FFFFFF", command = self.draw_data).grid(row = 6, column = 0, columnspan = 2, sticky = "NS", padx = self.p, pady = self.p)
        self.neighbor_progress = tkinter.ttk.Progressbar(left_part, value = 0, style = "success.Striped.Horizontal.TProgressbar")
        self.neighbor_progress.grid(row = 7, column = 0, sticky = "NSEW", padx = self.p, pady = self.p)
        self.sort_progress = tkinter.ttk.Progressbar(left_part, value = 0, style = "success.Striped.Horizontal.TProgressbar")
        self.sort_progress.grid(row = 7, column = 1, sticky = "NSEW", padx = self.p, pady = self.p)
        self.domain_progress = tkinter.ttk.Progressbar(left_part, value = 0, style = "success.Striped.Horizontal.TProgressbar")
        self.domain_progress.grid(row = 8, column = 0, sticky = "NSEW", padx = self.p, pady = self.p)
        self.draw_progress = tkinter.ttk.Progressbar(left_part, value = 0, style = "success.Striped.Horizontal.TProgressbar")
        self.draw_progress.grid(row = 8, column = 1, sticky = "NSEW", padx = self.p, pady = self.p)
        self.add_colors = pyCOGNAT_basic.TextFrameWithLabel(left_part, self.p, self.host.back, "#000000", "Add. domain colors:", ("Arial", 10, "bold"), None)
        self.add_colors.grid(row = 9, column = 0, columnspan = 2, sticky = "NSEW")
        central_panel.add(left_part)

        central_part = tkinter.Frame(central_panel)
        central_part.grid_rowconfigure(1, weight = 1)
        central_part.grid_columnconfigure(0, weight = 1)
        top_central_frame = tkinter.Frame(central_part)
        top_central_frame.grid_rowconfigure(0, weight = 1)
        top_central_frame.grid_columnconfigure(9, weight = 1)
        self.color_type = tkinter.StringVar()
        self.color_type.set("Single domain")
        tkinter.Label(top_central_frame, text = "Color mode:", justify = tkinter.LEFT).grid(row = 0, column = 0, sticky = "NSW", padx = self.p, pady = self.p)
        color_type_choice = tkinter.ttk.Combobox(top_central_frame, width = 20, state = "readonly", textvariable = self.color_type)
        color_type_choice["values"] = ("Single domain", "Domain architecture")
        color_type_choice.grid(row = 0, column = 1, sticky = "NSW", padx = self.p, pady = self.p)
        color_type_choice.current(0)
        color_type_choice.bind("<<ComboboxSelected>>", self.refresh)
        self.show_org = tkinter.BooleanVar()
        self.show_org.set(True)
        b = tkinter.Checkbutton(top_central_frame, text = "Show organism", variable = self.show_org)
        b.grid(row = 0, column = 2, sticky ="NSW", padx = self.p, pady = self.p)
        tkinter.Label(top_central_frame, text = "X axis mult.:", justify = tkinter.LEFT).grid(row = 0, column = 3, sticky = "NSW", padx = self.p, pady = self.p)
        self.x_mult_widget = tkinter.Entry(top_central_frame, width = 4, font = ("Courier New", 10))
        self.x_mult_widget.insert(tkinter.END, "2.0")
        self.x_mult_widget.grid(row = 0, column = 4, sticky = "NSEW", padx = self.p, pady = self.p)
        tkinter.Label(top_central_frame, text = "Y axis mult.:", justify = tkinter.LEFT).grid(row = 0, column = 5, sticky = "NSW", padx = self.p, pady = self.p)
        self.y_mult_widget = tkinter.Entry(top_central_frame, width = 4, font = ("Courier New", 10))
        self.y_mult_widget.insert(tkinter.END, "0.5")
        self.y_mult_widget.grid(row = 0, column = 6, sticky = "NSEW", padx = self.p, pady = self.p)
        tkinter.Label(top_central_frame, text = "Label offset:", justify = tkinter.LEFT).grid(row = 0, column = 7, sticky = "NSW", padx = self.p, pady = self.p)
        self.offset_widget = tkinter.Entry(top_central_frame, width = 4, font = ("Courier New", 10))
        self.offset_widget.insert(tkinter.END, "0")
        self.offset_widget.grid(row = 0, column = 8, sticky = "NSW", padx = self.p, pady = self.p)
        self.redraw_button = tkinter.Button(top_central_frame, text = "Redraw", font = ("Arial", 10, "bold"), command = self.redraw, state = tkinter.DISABLED)
        self.redraw_button.grid(row = 0, column = 9, sticky = "NSW", padx = self.p, pady = self.p)
        self.export_legend = tkinter.BooleanVar()
        self.export_legend.set(True)
        c = tkinter.Checkbutton(top_central_frame, text = "Export legend", variable = self.export_legend)
        c.grid(row = 0, column = 10, sticky ="NSE", padx = self.p, pady = self.p)
        tkinter.Button(top_central_frame, text = "Export", command = self.export).grid(row = 0, column = 11, sticky = "NSE", padx = self.p, pady = self.p)
        tax_export_status = tkinter.DISABLED
        if self.host.gbff_colors != None:
            tax_export_status = tkinter.NORMAL
        tkinter.Button(top_central_frame, text = "Export tax. colors", command = self.host.export_tax_colors, state = tax_export_status).grid(row = 0, column = 12, sticky = "NSE", padx = self.p, pady = self.p)
        tkinter.Button(top_central_frame, text = "Close", command = self.close).grid(row = 0, column = 13, sticky = "NSE", padx = self.p, pady = self.p)
        top_central_frame.grid(row = 0, column = 0, sticky = "NSEW", padx = self.p, pady = self.p + 1)
        self.main_canvas_frame = tkinter.Frame(central_part)
        self.main_canvas = get_canvas(self.main_canvas_frame, 1080, 1080)
        self.main_canvas_frame.grid(row = 1, column = 0, sticky = "NSEW", padx = self.p, pady = self.p)
        central_panel.add(central_part)

        right_part = tkinter.PanedWindow(central_panel, orient = tkinter.VERTICAL, sashwidth = self.p * 2, sashrelief = tkinter.RIDGE, background = self.host.back)
        top_right_part = tkinter.Frame(central_panel)
        top_right_part.grid_rowconfigure(1, weight = 1)
        top_right_part.grid_columnconfigure(0, weight = 1)
        self.action_mode = ActionModeFrame(top_right_part, self)
        self.action_mode.grid(row = 0, column = 0, sticky = "NSEW", padx = self.p, pady = self.p)
        self.report_widget = pyCOGNAT_basic.TextFrameWithLabel(top_right_part, self.p, self.host.back, "#000000", "Information about selected entity:", ("Arial", 10, "bold"), None)
        self.report_widget.wrap_text.set(False)
        self.report_widget.wrap_configure()
        self.report_widget.grid(row = 1, column = 0, sticky = "NSEW")
        self.action_mode.report_widget = self.report_widget
        right_part.add(top_right_part)

        bottom_right_frame = tkinter.Frame(right_part)
        bottom_right_frame.grid_rowconfigure(1, weight = 1)
        bottom_right_frame.grid_columnconfigure(0, weight = 1)
        tkinter.Label(bottom_right_frame, text = "Min. occurence:", justify = tkinter.RIGHT).grid(row = 0, column = 0, sticky = "NSE", padx = self.p, pady = self.p)
        self.min_occur = tkinter.Entry(bottom_right_frame, width = 5, font = ("Courier New", 10))
        self.min_occur.insert(tkinter.END, "3")
        self.min_occur.grid(row = 0, column = 1, sticky = "NSE", padx = self.p, pady = self.p)
        tkinter.Button(bottom_right_frame, text = "Random colors", command = self.random_colors).grid(row = 0, column = 2, sticky = "NSE", padx = self.p, pady = self.p)
        tkinter.Button(bottom_right_frame, text = "Export colors", command = self.export_colors).grid(row = 0, column = 3, sticky = "NSE", padx = self.p, pady = self.p)
        legend_canvas_frame = tkinter.Frame(bottom_right_frame)
        self.legend_canvas = get_canvas(legend_canvas_frame, 200, 300)
        legend_canvas_frame.grid(row = 1, column = 0, columnspan = 4, sticky = "NSEW")
        right_part.add(bottom_right_frame)
        central_panel.add(right_part)

        self.update_idletasks()
        central_panel.sash_place(0, 225, 1)
        central_panel.sash_place(1, 1450, 1)
        right_part.sash_place(0, 1, 700)

    def export_colors(self):
        curr_min_occur = self.min_occur.get().strip()
        how_many = None
        try:
            how_many = int(curr_min_occur)
        except ValueError:
            pass
        export_strings = list()
        if self.domain_occur != None: # Domains were already obtained
            ordered_by_occur = sorted(list(self.domain_occur.keys()), key = lambda k: self.domain_occur[k], reverse = True)
            for domain in ordered_by_occur:
                if how_many != None:
                    curr_occur = self.domain_occur[domain]
                    if curr_occur < how_many:
                        break
                if domain in self.domain_to_color:
                    export_strings.append("%s\t%s" % (domain, self.domain_to_color[domain]))
            self.add_colors.add_string("\n")
            self.add_colors.add_string("\n".join(export_strings))
            if how_many == None:
                self.host.set_status("DONE", "Current colors for all domains were exported", "green")
            else:
                self.host.set_status("DONE", "Current colors for domains occuring at least %i times were exported" % how_many, "green")
        else:
            self.host.set_status("ERROR", "Domains were", "green")

    def random_colors(self):
        curr_min_occur = self.min_occur.get().strip()
        new_colors = dict()
        if self.domain_occur != None: # Domains were already obtained
            for domain in self.domain_occur.keys():
                if curr_min_occur != "":
                    curr_occur = self.domain_occur[domain]
                    if curr_occur < int(curr_min_occur):
                        continue
                RGB_range = range(0,255)
                rand_color = "#%02x%02x%02x" % (random.choice(RGB_range), random.choice(RGB_range), random.choice(RGB_range))
                new_colors[domain] = rand_color
            self.redraw(new_colors)
            self.host.set_status("DONE", "The canvas was remade with random colors", "green")
        else:
            self.host.set_status("ERROR", "Cannot make random colors for not an empty canvas", "red")

    def legend_doubleclicked(self, event):
        real_x = self.legend_canvas.canvasx(event.x) # This calculates canvas coordinates, not the screen
        real_y = self.legend_canvas.canvasy(event.y)
        item = self.legend_canvas.find_closest(real_x, real_y)
        tags = self.legend_canvas.itemcget(item, "tags")
        if tags != "":
            domain_id = tags.split(" ")[1]
            curr_color = self.base_domain_color
            if domain_id in self.domain_to_color:
                curr_color = self.domain_to_color[domain_id]
            new_color = tkinter.colorchooser.askcolor(curr_color, title = "Choose a color for the following domain: %s" % domain_id)
            new_color_hex = new_color[1]
            if new_color_hex != None:
                self.domain_to_color[domain_id] = new_color_hex
                self.legend_canvas.itemconfigure(item, fill = new_color_hex)
                self.redraw(self.domain_to_color)

    def on_mousewheel(self, event):
        scroll_value = None
        system = platform.system()
        if system == "Windows":
            scroll_value = int(-1*(event.delta/120))
        elif system == "Darwin":
            scroll_value = int(-1*event.delta)
        else:
            print ("WARNING: System defined by Python is '%s', cannot use <MouseWheel> event!")
            return
        self.main_canvas.yview_scroll(scroll_value, "units")

    def canvas_doubleclicked(self, event):
        real_x = self.main_canvas.canvasx(event.x)
        real_y = self.main_canvas.canvasy(event.y)
        item = self.main_canvas.find_closest(real_x, real_y)
        tags = self.main_canvas.itemcget(item, "tags")
        if tags != "":
            tags = tags.split(" ")
            if tags[0] == "domain":
                data = tags[1]
                curr_domain = data.split("..")[6]
                if not curr_domain in self.domain_to_color: # Domain was not colored
                    RGB_range = range(0,255)
                    rand_color = "#%02x%02x%02x" % (random.choice(RGB_range), random.choice(RGB_range), random.choice(RGB_range))
                    self.domain_to_color[curr_domain] = rand_color
                    self.redraw(self.domain_to_color)
                else:
                    self.domain_to_color.pop(curr_domain)
                    self.redraw(self.domain_to_color)
        self.canvas_clicked(event)

    def canvas_clicked(self, event):
        dash_pattern = (6, 4)
        if self.selected != None:
            if self.selected[1] == "gene":
                self.main_canvas.itemconfigure(self.selected[0], outline = "#000000", dash = "")
            if self.selected[1] == "intergene":
                self.main_canvas.itemconfigure(self.selected[0], fill = "#000000", dash = "")
            if self.selected[1] == "domain":
                self.main_canvas.itemconfigure(self.selected[0], outline = "", dash = "")
        self.report_widget.clear()
        real_x = self.main_canvas.canvasx(event.x)
        real_y = self.main_canvas.canvasy(event.y)
        item = self.main_canvas.find_closest(real_x, real_y)
        tags = self.main_canvas.itemcget(item, "tags")
        domain_to_name = None
        if self.domain_type.get() == "COG":
            domain_to_name = self.host.cog_to_name
        if self.domain_type.get() == "Pfam":
            domain_to_name = self.host.pfam_to_name

        if tags != "":
            tags = tags.split(" ")
            if tags[0] == "gene":
                self.main_canvas.itemconfigure(item, outline = "red", dash = dash_pattern)
                data = tags[1]
                parts = data.split("..")
                try:
                    real_begin = int(parts[0])
                    real_end = int(parts[1])
                    real_direction = int(parts[2])
                    nucl_seq = parts[3]
                    protein_id = parts[4].strip("*")
                    curr_gbk = self.proteins[protein_id].source_record
                    curr_org, curr_gbff, curr_taxonomy = "Unk", "Unk", "Unk"
                    if curr_gbk in self.nucl_id_to_data:
                        curr_org = self.nucl_id_to_data[curr_gbk][0]
                        curr_gbff = self.nucl_id_to_data[curr_gbk][1]
                        curr_taxonomy = self.nucl_id_to_data[curr_gbk][2]
                    domain_data = None
                    if protein_id in self.id_to_domains:
                        domain_data = self.id_to_domains[protein_id]
                    pyCOGNAT_basic.report_gene_to_widget(self.report_widget, self.proteins[protein_id], real_begin, real_end, real_direction, nucl_seq, curr_org, curr_gbff, curr_taxonomy, domain_data, domain_to_name)
                except KeyError as e:
                    warning = "WARNING: you clicked on a gene with data '%s', no data for it was saved!" % data
                    self.report_widget.add_string("%s\n" % warning)
                    #self.report_widget.add_string("%s\n" % self.id_to_domains)
                    print(warning)
                    raise e
                except IndexError as e:
                    warning = "WARNING: you clicked on a gene with data '%s', but is contain less than five data pieces!" % data
                    self.report_widget.add_string("%s\n" % warning)
                    print(warning)
                    raise e
                except TypeError as e:
                    warning = "WARNING: you clicked on a gene with data '%s', but it cannot be interpreted correctly!" % data
                    self.report_widget.add_string("%s\n" % warning)
                    print(warning)
                    raise e
            if tags[0] == "intergene":
                data = tags[1]
                self.main_canvas.itemconfigure(item, fill = "red", dash = dash_pattern)
                pyCOGNAT_basic.report_intergene_to_widget(self.report_widget, data)
            if tags[0] == "domain":
                self.main_canvas.itemconfigure(item, outline = "red", dash = dash_pattern)
                data = tags[1]
                protein_id = data.split("..")[-1].strip("*")
                try:
                    pyCOGNAT_basic.report_domain_to_widget(self.report_widget, self.proteins[protein_id], data, domain_to_name)
                except KeyError:
                    warning = "WARNING: you clicked on a domain '%s', no data for a corresponding protein was saved!" % data
                    self.report_widget.add_string("%s\n" % warning)
                    print(warning)
            self.selected = (item, tags[0])

    def draw_data(self):
        # 1) Renamning the tab
        self.rename_tab()
        # 2) Obtaining parameters
        protein_ids = self.clade_ids.get_content_as_list()
        self.protein_ids_safe = self.clade_ids.get_content_as_list()
        protein_num = len(protein_ids)
        neighbor_number = int(self.neighbor_number_widget.get().strip())
        max_neighbor_number = self.host.project_data_tab.get_max_gene_number()
        if max_neighbor_number < neighbor_number:
            print ("[..ERROR..] Too much neighbors requested (%i), only %i available" % (neighbor_number, max_neighbor_number))
            self.host.set_status("ERROR", "Too much neighbors requested, see console for details", "red")
            return
        self.host.set_status("STARTED", "[1 out of 2] Loading neighborhood information...", "#f08650")
        project_name = self.host.project_data_tab.get_project_name()
        script_path = self.host.settings.script_dir
        sliced_database_path = os.path.join(self.host.settings.work_dir, project_name, "COGNAT_database")
        ids_to_gbff_and_gbk = None
        self.protein_id_to_nucl_data = dict()
        self.neighborhoods = None
        self.pids_in_neighborhoods = None
        self.id_to_domains = dict() # All proteins from the database slice will be held here
        self.proteins = dict()
        self.domain_occur = None
        self.drawable_objects = None
        self.nucl_id_to_data = dict() # gbk id to a tuple of organism name, gbff and taxonomy
        self.neighbor_progress["value"] = 0
        self.sort_progress["value"] = 0
        self.domain_progress["value"] = 0
        self.draw_progress["value"] = 0
        # 3) Obtaining neighborhood obejcts
        try:
            sys.path.append(script_path)
            import COGNAT_basic
            ids_to_gbff_and_gbk = COGNAT_basic.map_proteins_to_gbff_and_gbk(protein_ids, sliced_database_path, self.neighbor_progress, self.host.parent)
            (self.neighborhoods, self.pids_in_neighborhoods) = COGNAT_basic.get_sorted_neighborhoods(self.protein_ids_safe, ids_to_gbff_and_gbk, sliced_database_path, neighbor_number, True, list(), self.sort_progress, self.host.parent)
            evalue = self.host.project_data_tab.get_evalue()
            i = 0
            for p_id in ids_to_gbff_and_gbk.keys():
                i += 1
                self.domain_progress["value"] = int(100 * i / len(ids_to_gbff_and_gbk.keys()))
                self.host.parent.update()
                gbff = ids_to_gbff_and_gbk[p_id][0]
                gbk = ids_to_gbff_and_gbk[p_id][1]
                COGNAT_basic.assign_proteins_to_domains(self.id_to_domains, sliced_database_path, gbff, gbk, None, evalue, self.domain_type.get())
                self.protein_id_to_nucl_data[p_id] = COGNAT_basic.read_nucl_data(sliced_database_path, gbff, gbk)
                fst_nucl_data = self.protein_id_to_nucl_data[p_id][0] # Filling <self.nucl_id_to_data>
                if not fst_nucl_data.gi in self.nucl_id_to_data:
                    self.nucl_id_to_data[fst_nucl_data.gi] = (fst_nucl_data.organism, gbff, fst_nucl_data.taxonomy)
            # Removing proteins which do not belong to current neighborhood number
            filtered_dict = dict()
            for p_id in self.id_to_domains.keys():
                if p_id in self.pids_in_neighborhoods:
                    filtered_dict[p_id] = self.id_to_domains[p_id]
            self.id_to_domains = filtered_dict
            import udav_soft
            if self.host.project_data_tab.unite_domains.get() == True: # Domain should be united under given options
                max_distance = self.host.project_data_tab.get_max_distance()
                max_hmm_overlap = self.host.project_data_tab.get_max_hmm_overlap()
                print ("Hits of the same domain will be united if a) they are separated by less than %i residues and b) regions in a profile HMM which cover its hits overlap by less than %.2f percent" % (max_distance, max_hmm_overlap))
                self.id_to_domains = udav_soft.unite_same_Pfam_hits(self.id_to_domains, max_distance, max_hmm_overlap)
            overlap = self.host.project_data_tab.get_overlap()
            self.id_to_domains = udav_soft.filter_Pfam_hits(self.id_to_domains, overlap)
            del COGNAT_basic
            del udav_soft
            self.host.set_status("STARTED", "[2 out of 2] Obtained %i neighborhoods for %i protein ids, drawing data..." % (len(self.neighborhoods), protein_num), "#f08650")
        except ImportError:
            print ("No module 'COGNAT_basic' and/or its dependencies was found in the script directory")
            return
        except ValueError:
            error_text = "Non-int value given as required neighbor number: '%s'" % neighbor_number
            print (error_text)
            self.host.set_status("ERROR", error_text, "red")
            raise
            return
        except IndexError:
            error_text = "It is likely that domain to color assignment was corrupted, it should be tab-separated data"
            print (error_text)
            self.host.set_status("ERROR", error_text, "red")
            return
        self.redraw()
        self.redraw_button.configure(state = tkinter.NORMAL, background = self.host.header, foreground = "#FFFFFF")

    def redraw(self, preset_colors = None):
        self.draw_progress["value"] = 0
        x_mult = float(self.x_mult_widget.get())
        y_mult = float(self.y_mult_widget.get())
        nucl_per_pixel = x_mult * 5
        arrow_height = y_mult * 50
        v_sep = 2
        top_left_y = v_sep
        n = 0
        lines_of_objects = list() # list of lists of pyCOGNAT_basic.DrawableObject
        self.domain_architectures = dict() # Domain architecture to a list of lengths of proteins which have this architecture
        self.domain_to_color = None
        if preset_colors == None: # Colors are not given directly
            self.domain_to_color = self.host.project_data_tab.get_domain_to_color()
            add_colors = self.add_colors.get_content_as_dict("\t") # These overrides global
            for c in add_colors.keys():
                self.domain_to_color[c] = add_colors[c]
        else:
            self.domain_to_color = preset_colors

        for neighborhood in self.neighborhoods:
             #print ("n = %i, protein_id = '%s', neighborhood: '%s'" % (n, self.protein_ids_safe[n], neighborhood))
             for p in neighborhood:
                 pid_stripped = p.protein_id.strip("*")
                 if not pid_stripped in self.proteins:
                     self.proteins[pid_stripped] = p
             curr_top_left_y = top_left_y + ((arrow_height + v_sep) * n)
             nucl_data_list = list()
             if self.protein_ids_safe[n] in self.protein_id_to_nucl_data: # This is a real protein found in the database
                 nucl_data_list = self.protein_id_to_nucl_data[self.protein_ids_safe[n]]
             new_line_of_objects = self.prepare_neighborhood(curr_top_left_y, arrow_height, nucl_per_pixel, neighborhood, nucl_data_list, self.id_to_domains, self.domain_to_color, self.pids_in_neighborhoods)
             lines_of_objects.append(new_line_of_objects)
             n += 1
             self.draw_progress["value"] = int(100 * n / len(self.neighborhoods))

        protein_ids_to_show = self.protein_ids_safe
        colors = None
        if self.host.gbff_colors != None:
            colors = list()
            for i in range(len(protein_ids_to_show)):
                color = None
                try:
                    curr_gbk = self.proteins[protein_ids_to_show[i]].source_record
                    if curr_gbk in self.nucl_id_to_data:
                        curr_gbff = self.nucl_id_to_data[curr_gbk][1]
                        if curr_gbff in self.host.gbff_colors:
                            color = self.host.gbff_colors[curr_gbff]
                except KeyError:
                    print ("Internat error in the redraw() method, cannot detect proper color for the organism")
                if color == None:
                    if "Unknown" in self.host.gbff_colors:
                        color = self.host.gbff_colors["Unknown"]
                    else:
                        color = "#888888"
                colors.append(color)
        if self.show_org.get() == True:
            protein_ids_to_show = list()
            for i in range(len(self.protein_ids_safe)):
                protein_ids_to_show.append(self.protein_ids_safe[i])
                curr_gbk = self.proteins[protein_ids_to_show[i]].source_record
                curr_org = "Unk"
                if curr_gbk in self.nucl_id_to_data: # Real protein
                    curr_org = self.nucl_id_to_data[curr_gbk][0]
                protein_ids_to_show[i] += " [%s]" % curr_org
        self.draw_lines_of_objects(lines_of_objects, arrow_height, protein_ids_to_show, colors)
        self.drawable_objects = lines_of_objects
        self.domain_occur = self.draw_domain_legend(2 * nucl_per_pixel / 3, arrow_height / 2)
        self.host.set_status("DONE", "Obtained %i lines of drawable objects" % len(self.drawable_objects), "green")

    def prepare_neighborhood(self, top_left_y, arrow_height, nucl_per_pixel, neighborhood, nucl_data_list, id_to_domains, domain_to_color, pids_in_neighborhoods):
        """
        This method will create a list of <pyCOGNAT_basic.DrawableObject> which will represent all objects
        for current <neighborhood>. Uses information from <pids_in_neighborhoods>: it is True if protein
        direction was changed and False otherwise
        """
        objects_to_draw = list()
        curr_x0 = 0
        #print ("    Working with new neighborhood:")
        fst_gene_start = neighborhood[0].gene_begin
        gene_not_colored_list = list()
        for i in range(len(neighborhood)):
            # 1) Creating arrow & domain rectangles for current gene
            curr_gene = neighborhood[i]
            curr_gene_length = curr_gene.gene_end - curr_gene.gene_begin + 1
            curr_prot_length = len(curr_gene.sequence)
            curr_direction = curr_gene.gene_direction
            curr_domains = None
            stripped_pid = curr_gene.protein_id.strip("*")
            was_reversed = False
            if stripped_pid in pids_in_neighborhoods:
                was_reversed = pids_in_neighborhoods[stripped_pid]
            #print ("CURR: %s" % curr_gene.__dict__)
            if stripped_pid in id_to_domains:
                curr_domains = id_to_domains[stripped_pid]
                domain_arch = pyCOGNAT_basic.get_domain_architecture(curr_domains)
                if not domain_arch in self.domain_architectures:
                    self.domain_architectures[domain_arch] = list()
                self.domain_architectures[domain_arch].append(curr_prot_length)
            (begin_seq, end_seq, gene_sequence, record_id) = pyCOGNAT_basic.get_real_coordinates_and_sequence(curr_gene.gene_begin, curr_gene.gene_end, was_reversed, nucl_data_list)
            real_direction = curr_direction
            if was_reversed:
                real_direction = 0 - real_direction
            if real_direction == -1:
                gene_sequence = pyCOGNAT_basic.get_complement(gene_sequence)
            additional_data = "%s..%s..%s..%s..%s" % (begin_seq, end_seq, real_direction, gene_sequence, curr_gene.protein_id)
            color_by_architecture = False
            if self.color_type.get() == "Domain architecture":
                color_by_architecture = True
            (new_objects, gene_not_colored) = pyCOGNAT_basic.get_gene_arrow_and_domains(additional_data, curr_x0, top_left_y, curr_gene_length, curr_prot_length, nucl_per_pixel, arrow_height, curr_direction, curr_domains, domain_to_color, self.base_domain_color, color_by_architecture)
            gene_not_colored_list.append((stripped_pid, gene_not_colored))
            objects_to_draw.extend(new_objects)
            # 2) Creating intergene line between current and the next gene
            if i != len(neighborhood) - 1: # This is not the last gene
                next_gene = neighborhood[i + 1]
                intergene_begin = curr_gene.gene_end + 1
                intergene_end = next_gene.gene_begin - 1
                if intergene_begin <= intergene_end: # Which means that this intergene region exists
                    (intergene_begin_seq, intergene_end_seq, sequence, record_id) = pyCOGNAT_basic.get_real_coordinates_and_sequence(intergene_begin, intergene_end, was_reversed, nucl_data_list)
                    if was_reversed:
                        sequence = pyCOGNAT_basic.get_complement(sequence)
                    sequence = "%s..%i..%i..%s..%s..%s..%s" % (sequence, intergene_begin_seq, intergene_end_seq, record_id, was_reversed, stripped_pid, next_gene.protein_id.strip("*"))
                    line_x0 = float(intergene_begin - fst_gene_start) / nucl_per_pixel
                    line_x1 = float(intergene_end - fst_gene_start) / nucl_per_pixel
                    line_y = top_left_y + (arrow_height / 2)
                    objects_to_draw.append(pyCOGNAT_basic.DrawableObject("line", sequence, [line_x0, line_y, line_x1, line_y], "#000000", ""))
                curr_x0 = (next_gene.gene_begin - fst_gene_start) / nucl_per_pixel
        if self.hide_uncolored.get(): # Genes which are not colored and flank a neighborhood should be filtered from the <objects_to_draw>
            pyCOGNAT_basic.hide_uncolored(objects_to_draw, gene_not_colored_list)
        return objects_to_draw

    def draw_lines_of_objects(self, lines_of_objects, arrow_height, strings = None, colors = None):
        """
        This method places lists of <DrawableObject> given in <lines_of_objects> on the <self.main_canvas>
        """
        self.main_canvas.delete("all")
        # 1) Finding start of the first required gene
        fst_req_x_coord = None
        for single_line in lines_of_objects:
            if len(single_line) == 1:
                continue
            if fst_req_x_coord != None:
                break
            for obj in single_line:
                if (obj.object_type == "arrow") and (obj.data.count("*") != 0):
                    fst_req_x_coord = obj.coordinates[0]
                    break
        # 2) Shifting neighborhoods to the required gene, finding <min_x>
        min_x = 0
        for single_line in lines_of_objects:
            diff = 0
            for obj in single_line:
                if (obj.object_type == "arrow") and (obj.data.count("*") != 0):
                    curr_x = obj.coordinates[0] # Required gene is always facing in one (direct) direction
                    diff = fst_req_x_coord - curr_x
                    break
            for obj in single_line:
                obj.shift_by_x(diff)
                min_x = obj.get_best_value("x", min, min_x)
                #if obj.coordinates[0] < min_x:
                #    min_x = obj.coordinates[0]
        # 3) Shifting neighborhoods so that all x coordinates are positive, getting <max_x> and <max_y>
        max_x = 0
        max_y = 0
        for single_line in lines_of_objects:
            for obj in single_line:
                obj.shift_by_x(self.p - min_x)
                max_x = obj.get_best_value("x", max, max_x)
                max_y = obj.get_best_value("y", max, max_y)
                #if obj.coordinates[0] > max_x:
                #    max_x = obj.coordinates[0]
                #if obj.coordinates[1] > max_y:
                #    max_y = obj.coordinates[1]
        # 4) Drawing objects
        for single_line in lines_of_objects:
            for obj in single_line:
                #print (obj.__dict__)
                if obj.object_type == "arrow":
                    self.main_canvas.create_polygon(obj.coordinates, tags = ("gene", obj.data), fill = obj.fill, outline = obj.outline, width = 2)
                if obj.object_type == "line":
                    self.main_canvas.create_line(obj.coordinates, tags = ("intergene", obj.data), fill = obj.fill, width = 2)
                if obj.object_type == "domain":
                    self.main_canvas.create_rectangle(obj.coordinates, tags = ("domain", obj.data), fill = obj.fill, outline = obj.outline)
        # 5) Adding text information
        max_text_width = 0
        text_font = tkinter.font.Font(self.main_canvas, ("Courier New", 0 - int(arrow_height)))
        n = 0
        offset = 0
        try:
            offset = int(self.offset_widget.get().strip())
        except ValueError:
            print ("Text offset '%s' is not an integer value, please enter a proper number. Offset is set to 0" % self.offset_widget.get().strip())
        self.text_objects = list()
        if strings != None:
            for string in strings:
                x0 = max_x + self.p - offset
                y0 = (n * (arrow_height + 2))
                curr_color = "#888888"
                if colors != None:
                    try:
                        curr_color = colors.pop(0)
                    except IndexError:
                        print ("Warning: no color remained in the color list, default color was used")
                text_object = self.main_canvas.create_text(x0, y0, tags = ("text"), fill = curr_color, font = text_font, anchor = tkinter.NW, text = string)
                bounds = self.main_canvas.bbox(text_object)  # returns a tuple like (x1, y1, x2, y2)
                max_text_width = max(bounds[2] - bounds[0], max_text_width)
                text_height = bounds[3] - bounds[1]
                self.text_objects.append([x0, y0, string, text_height])
                n += 1
        # 6) Canvas size update
        max_x += self.p + max_text_width
        max_y += self.p
        self.main_canvas.configure(width = max_x, height = max_y, scrollregion = (0, 0, max_x, max_y))
        self.max_x = max_x
        self.max_y = max_y

    def draw_domain_legend(self, aa_per_pixel, element_height):
        """
        The method draws a <self.legend_canvas> object
        """
        self.legend_canvas.delete("all")
        self.legend_objects = list()
        self.legend_text_objects = list()
        occured_in_proteins = dict()
        total_size_of_hits = dict()
        # 1) Calculating number of hits and total hit size
        if self.color_type.get() == "Single domain":
            for pid in self.id_to_domains:
                for domain in self.id_to_domains[pid].keys():
                    if not domain in occured_in_proteins:
                        occured_in_proteins[domain] = 0
                        total_size_of_hits[domain] = 0
                    regions = self.id_to_domains[pid][domain].strip().split(" ")
                    for r in regions:
                        parts = r.split("..")
                        r_length = int(parts[1]) - int(parts[0]) + 1
                        occured_in_proteins[domain] += 1
                        total_size_of_hits[domain] += r_length
        if self.color_type.get() == "Domain architecture":
            for domain_arch in self.domain_architectures.keys():
                if not domain_arch in occured_in_proteins:
                    occured_in_proteins[domain_arch] = 0
                    total_size_of_hits[domain_arch] = 0
                protein_lengths = self.domain_architectures[domain_arch]
                for plen in protein_lengths:
                    occured_in_proteins[domain_arch] += 1
                    total_size_of_hits[domain_arch] += plen
        sorted_by_occur = list(occured_in_proteins.keys())
        sorted_by_occur.sort(key = lambda k: occured_in_proteins[k], reverse = True)

        # 2) Printing legend
        text_font = tkinter.font.Font(self.legend_canvas, ("Courier New", 0 - int(element_height)))
        bold_font = tkinter.font.Font(self.legend_canvas, ("Courier New", 0 - int(element_height), "bold"))
        number_width = text_font.measure("12345")
        y0 = self.p
        x0 = self.p
        i = 0
        max_width = 0
        for domain in sorted_by_occur:
            y0 = i * ((element_height) + self.p)
            domain_size = (total_size_of_hits[domain] / occured_in_proteins[domain]) / aa_per_pixel
            self.legend_canvas.create_text(x0, y0, fill = "#888888", font = text_font, anchor = tkinter.NW, text = "%i" % occured_in_proteins[domain])
            self.legend_text_objects.append([x0, y0, "%i" % occured_in_proteins[domain], int(element_height), domain])
            fill_color = self.base_domain_color
            if domain in self.domain_to_color:
                fill_color = self.domain_to_color[domain]
            domain_rectangle_x_end = x0 + number_width + self.p + int(domain_size)
            self.legend_canvas.create_rectangle(x0 + number_width + self.p, y0, domain_rectangle_x_end, y0 + element_height, tags = ("domain", domain), fill = fill_color, outline = "")
            self.legend_objects.append([x0 + number_width + self.p, y0, domain_rectangle_x_end, y0 + element_height, domain])
            max_width = max(domain_rectangle_x_end, max_width)
            i += 1

        y0 = self.p
        x0 = self.p + max_width
        i = 0
        max_descr_text_width = 0
        for domain in sorted_by_occur:
            y0 = i * ((element_height) + self.p)
            curr_domain_name = None
            curr_font = text_font
            try:
                if self.domain_type.get() == "COG":
                    curr_domain_name = self.host.cog_to_name[domain]
                if self.domain_type.get() == "Pfam":
                    curr_domain_name = self.host.pfam_to_name[domain]
            except KeyError:
                if self.color_type.get() == "Single domain":
                    curr_domain_name = "Domain name not found"
                else:
                    curr_domain_name = "Domain architecture"
                    curr_font = bold_font
            domain_text = "%s [%s]" % (domain, curr_domain_name)
            descr_text = self.legend_canvas.create_text(x0, y0, fill = "#888888", font = curr_font, anchor = tkinter.NW, text = domain_text)
            self.legend_text_objects.append([x0, y0, domain_text, int(element_height), domain])
            bounds = self.legend_canvas.bbox(descr_text)  # returns a tuple like (x1, y1, x2, y2)
            max_descr_text_width = max(bounds[2] - bounds[0], max_descr_text_width)
            i += 1

        max_x = max_width + max_descr_text_width + self.p
        max_y = len(sorted_by_occur) * (element_height + self.p)
        self.legend_canvas.configure(width = max_x, height = max_y, scrollregion = (0, 0, max_x, max_y))
        self.legend_max_x = max_x
        self.legend_max_y = max_y
        return occured_in_proteins

    def export(self):
        try:
            import drawsvg
        except ModuleNotFoundError:
            self.host.set_status("ERROR", "Failed to import <drawsvg> module, export cannot be done", "red")
            return
        if self.drawable_objects == None:
            self.host.set_status("ERROR", "No drawing was yet obtained, export cannot be done", "red")
            return

        project_name = self.host.project_data_tab.get_project_name()
        clade_name = self.clade_name_widget.get().strip()
        default_filename = "%s.%s" % (project_name, clade_name)
        svg_filename = tkinter.filedialog.asksaveasfilename(filetypes = (("SVG files", "*.svg"), ("All", "*.*")), title = "Please enter a filename for the SVG scheme", initialfile = default_filename)
        if svg_filename == "": # Cancel
            return
        if re.search("\.svg$", svg_filename) == None:
            svg_filename += ".svg"
        export_w = self.max_x
        export_h = self.max_y
        if self.export_legend.get() == True:
            n_to_export = 0
            for obj in self.legend_objects:
                #obj = [x0, y0, x1, y1, domain_id]
                if obj[4] in self.domain_to_color:
                    n_to_export += 1
            export_w = max(self.max_x, self.legend_max_x)
            export_h += (self.legend_max_y * n_to_export / len(self.legend_objects)) + 25

        drawing = drawsvg.Drawing(export_w, export_h, origin = (0, 0))
        # 1) Drawing arrows, domains and lines
        for object_line in self.drawable_objects:
            line = drawsvg.Group()
            for obj in object_line:
                c = obj.coordinates
                if obj.object_type == "arrow":
                    arrow = drawsvg.Lines(c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12], c[13], close = True, fill = obj.fill, stroke = obj.outline)
                    line.append(arrow)
                if obj.object_type == "domain":
                    rh = c[3] - c[1]
                    rw = c[2] - c[0]
                    line.append(drawsvg.Rectangle(c[0], c[1], rw, rh, fill = obj.fill))
                if obj.object_type == "line":
                    line.append(drawsvg.Line(c[0], c[1], c[2], c[3], stroke = obj.fill))
            drawing.append(line)
        # 2) Drawing text
        text_data = drawsvg.Group()
        for text in self.text_objects:
            #text = [x0, y0, string, text_height]
            text_data.append(drawsvg.Text(text[2], text[3], x = text[0], y = text[1], fill = "#888888", text_anchor = "start", dominant_baseline = "hanging", font_family="Courier New"))
        drawing.append(text_data)

        # 3) Drawing legend
        if self.export_legend.get() == True:
            legend_data = drawsvg.Group()
            legend_padding = self.max_y + 20
            skipped_y = list()
            for obj in self.legend_objects:
                #obj = [x0, y0, x1, y1, domain_id]
                if obj[4] in self.domain_to_color:
                    curr_fill = self.domain_to_color[obj[4]]
                    rh = obj[3] - obj[1]
                    rw = obj[2] - obj[0]
                    curr_y = obj[1]
                    if len(skipped_y) != 0: # Skipped elements exist
                        curr_y = skipped_y.pop(0)
                        skipped_y.append(obj[1])
                    legend_data.append(drawsvg.Rectangle(obj[0], curr_y + legend_padding, rw, rh, fill = curr_fill))
                else:
                    skipped_y.append(obj[1])

            skipped_y = list()
            names_started = False
            for text in self.legend_text_objects:
                #text = [x0, y0, string, text_height, domain]
                if (text[2].count("[") != 0) and (not names_started): # This is not occurence count, but domain description started
                    skipped_y = list()
                    names_started = True
                if text[4] in self.domain_to_color:
                    curr_y = text[1]
                    if len(skipped_y) != 0: # Skipped elements exist
                        curr_y = skipped_y.pop(0)
                        skipped_y.append(text[1])
                    legend_data.append(drawsvg.Text(text[2], text[3], x = text[0], y = curr_y + legend_padding, fill = "#888888", text_anchor = "start", dominant_baseline = "hanging", font_family = "Courier New"))
                else:
                    skipped_y.append(text[1])
            drawing.append(legend_data)
        drawing.save_svg(svg_filename)
        del drawsvg
        self.host.set_status("DONE", "Saved information from clade %s to a file '%s'" % (clade_name, svg_filename), "green")

    def rename_tab(self):
        clade_name = self.clade_name_widget.get().strip()
        if clade_name == "":
            return
        if clade_name.count(" ") != 0:
            clade_name = clade_name.replace(" ", "_")
            self.host.set_status("WARNING", "Clade name contained spaces which were replaced with the underlines", "#f08650")
        curr_index = self.host.tabs.index("current")
        self.host.tabs.tab(curr_index, text = clade_name)

    def close(self):
        curr_tab_index = self.host.tabs.index("current")
        self.host.tabs.forget(curr_tab_index)
        self.host.clade_tabs.pop(curr_tab_index)

    def refresh(self, event):
        print ("Selected domain type is: '%s'" % self.domain_type.get())
        print ("Selected coloring type is: '%s'" % self.color_type.get())
        if self.neighborhoods != None:
            self.redraw()

    def save(self):
        project_name = self.host.project_data_tab.get_project_name()
        clade_name = self.clade_name_widget.get().strip()
        if clade_name == "":
            print ("[..WARNING..] a clade failed to save because no name was specified for it!")
            raise ValueError
        if pyCOGNAT_basic.check_filename(clade_name) == False:
            print ("[..WARNING..] a clade failed to save because its name is unproper!")
            raise ValueError
        clade_ids_filename = os.path.join(self.host.settings.work_dir, project_name, "%s.%s.ids" % (project_name, clade_name))
        self.clade_ids.write_into_file(clade_ids_filename)
        add_colors_filename = os.path.join(self.host.settings.work_dir, project_name, "%s.%s.colors" % (project_name, clade_name))
        self.add_colors.write_into_file(add_colors_filename)

    def load(self, clade_name, clade_ids_filename, add_colors_filename):
        self.clade_ids.read_from_file(clade_ids_filename)
        if add_colors_filename != None:
            self.add_colors.read_from_file(add_colors_filename)
        self.clade_name_widget.delete(0, tkinter.END)
        self.clade_name_widget.insert(tkinter.END, clade_name)

class ActionModeFrame(tkinter.Frame):
    def __init__(self, parent, host):
        tkinter.Frame.__init__(self, parent, host)
        self.host = host
        self.p = host.p
        self.report_widget = None
        self.mode = None # tkinter.StringVar
        self.images = dict()
        self.buttons = dict()
        self.create_UI()

    def create_UI(self):
        bborder = 5
        self.images["i_button"] = tkinter.PhotoImage(file = r"i_button.png")
        self.images["this_button"] = tkinter.PhotoImage(file = r"this_button.png")
        self.images["g_5button"] = tkinter.PhotoImage(file = r"5g_button.png")
        self.images["g_3button"] = tkinter.PhotoImage(file = r"3g_button.png")
        self.images["i_5button"] = tkinter.PhotoImage(file = r"5i_button.png")
        self.images["i_3button"] = tkinter.PhotoImage(file = r"3i_button.png")

        self.mode = tkinter.StringVar(self, "i")
        self.buttons["i_button"] = tkinter.Radiobutton(self, text = "i_button", indicator = 0, variable = self.mode, value = "i", image = self.images["i_button"], command = self.show_info)
        self.buttons["i_button"].grid(row = 0, column = 0, sticky = "NS", padx = 1, pady = self.p)
        self.buttons["this_button"] = tkinter.Radiobutton(self, text = "this_button", indicator = 0, variable = self.mode, value = "this", image = self.images["this_button"], command = self.show_info, state = tkinter.DISABLED)
        self.buttons["this_button"].grid(row = 0, column = 1, sticky = "NS", padx = 1, pady = self.p)
        self.buttons["g_5button"] = tkinter.Radiobutton(self, text = "g_5button", indicator = 0, variable = self.mode, value = "g_5'", image = self.images["g_5button"], command = self.show_info, state = tkinter.DISABLED)
        self.buttons["g_5button"].grid(row = 0, column = 2, sticky = "NS", padx = 1, pady = self.p)
        self.buttons["g_3button"] = tkinter.Radiobutton(self, text = "g_3button", indicator = 0, variable = self.mode, value = "g_3'", image = self.images["g_3button"], command = self.show_info, state = tkinter.DISABLED)
        self.buttons["g_3button"].grid(row = 0, column = 3, sticky = "NS", padx = 1, pady = self.p)
        self.buttons["i_5button"] = tkinter.Radiobutton(self, text = "i_5button", indicator = 0, variable = self.mode, value = "i_5'", image = self.images["i_5button"], command = self.show_info, state = tkinter.DISABLED)
        self.buttons["i_5button"].grid(row = 0, column = 4, sticky = "NS", padx = 1, pady = self.p)
        self.buttons["i_3button"] = tkinter.Radiobutton(self, text = "i_3button", indicator = 0, variable = self.mode, value = "i_3'", image = self.images["i_3button"], command = self.show_info, state = tkinter.DISABLED)
        self.buttons["i_3button"].grid(row = 0, column = 5, sticky = "NS", padx = 1, pady = self.p)

    def show_info(self):
        self.report_widget.clear()
        self.report_widget.add_string("Current mode is: '%s'" % self.mode.get())
        return

#!/usr/bin/env python
"""
This is a main script of the <pyCOGNAT> program
@ Daria Dibrova aka udavdasha
"""
curr_version = "0.0.1"
import tkinter
import tkinter.messagebox as tkMessageBox
import tkinter.filedialog
import tkinter.ttk as ttk
import sys, os, platform, re
import Settings
import ProjectData, CladeFrame
import pyCOGNAT_basic

LOGO_FILENAME = "pyCOGNAT_logo.gif"
ICON_FILENAME = "pyCOGNAT.ico"
INI_FILENAME = "settings.ini"
if platform.system() == "Linux":
    ICON_FILENAME = "@pyCOGNAT.xbm"
if len(sys.argv) > 1: #settings filename was given as an argument
    INI_FILENAME = sys.argv[1]

class pyCOGNAT(tkinter.Frame):
    """
    Main window of the pyCOGNAT program
    """
    def __init__(self, parent, settings_filename):
        tkinter.Frame.__init__(self, parent)
        self.parent = parent
        self.parent.protocol("WM_DELETE_WINDOW", self.close_app)
        self.p = 5 # Padding for internal frames
        self.back = "#d9d9d9"
        self.header = "#ff0066"
        self.status_label = None
        self.tabs = None
        self.clade_tabs = None
        self.project_data_tab = None
        self.cog_to_name = None
        self.pfam_to_name = None
        self.gbff_colors = None
        self.gbff_order = None

        req_settings = ["script_dir", "work_dir", "cognat_database", "cog_descr_file", "pfam_descr_file"]
        self.settings = Settings.read_settings_file(settings_filename, req_settings)
        valid_settings = self.check_settings()
        if valid_settings:
            self.cog_to_name = self.read_profile_descr_file(self.settings.cog_descr_file)
            self.pfam_to_name = self.read_profile_descr_file(self.settings.pfam_descr_file)
            try:
                (self.gbff_colors, self.gbff_order) = pyCOGNAT_basic.read_taxonomy_colors(self.settings.tax_colors)
            except AttributeError:
                print ("Taxonomy coloring filename was not provided, names will not be colored according to the taxonomy")
        self.create_UI()
        #self.load_project(os.path.join(self.settings.work_dir, "COG1883"))

    def create_UI(self):
        self.grid_rowconfigure(0, weight = 1)
        self.grid_columnconfigure(0, weight = 1)

        self.tabs = ttk.Notebook(self)
        self.tabs.grid(row = 0, column = 0, sticky = "NSEW", padx = self.p, pady = self.p)

        self.clade_tabs = list()

        self.project_data_tab = ProjectData.ProjectData(self.tabs, self)
        self.tabs.add(self.project_data_tab, text = "Project data")

        self.status_label = tkinter.Label(self, text = "", foreground = self.header, font = ("Courier New", 12, "bold"))
        self.set_status("READY", "Enter a new project name or load an existing project", self.header)
        self.status_label.grid(row = 1, column = 0, sticky = "NSE", padx = self.p, pady = self.p)

    def set_status(self, status, comment = "", curr_color = "#FFFFFF"):
        text_to_show = "[..%s..] %s" % (status, comment)
        self.status_label.configure(text = text_to_show, foreground = curr_color)
        self.update()

    def read_profile_descr_file(self, filename):
        strings = pyCOGNAT_basic.read_plain_list(filename)
        destination = dict()
        for string in strings:
            fields = string.split("\t")
            destination[fields[0]] = fields[1]
        return destination

    def export_tax_colors(self):
        svg_filename = tkinter.filedialog.asksaveasfilename(filetypes = (("SVG files", "*.svg"), ("All", "*.*")), title = "Please enter a filename for the SVG scheme", initialfile = "taxonomy_colors")
        if svg_filename == "": # Cancel
            return
        if re.search("\.svg$", svg_filename) == None:
            svg_filename += ".svg"
        try:
            pyCOGNAT_basic.print_simple_legend(self.gbff_colors, self.gbff_order, svg_filename)
        except ImportError:
            print ("Failed to export taxonomy colors, <drawsvg> module cannot be loaded")

    def check_settings(self): # Write
        return True

    def create_new_project(self):
        curr_project_name = self.project_data_tab.get_project_name()
        if curr_project_name == "":
            self.set_status("ERROR", "Enter new project name before creating one", "red")
            return
        if pyCOGNAT_basic.check_filename(curr_project_name) == False:
            self.set_status("ERROR", "Project name contain non-alphanumeric characters, project cannot be created", "red")
            return

        max_gene_number = self.project_data_tab.get_max_gene_number()
        try:
            max_gene_number = int(max_gene_number)
        except ValueError:
            self.set_status("ERROR", "Max. gene number should be integer value", "red")
            return

        new_project_path = os.path.join(self.settings.work_dir, curr_project_name)
        if os.path.isdir(new_project_path):
            self.set_status("ERROR", "Project with this name already exists", "red")
            return
        os.mkdir(new_project_path)

        ids_filename = tkinter.filedialog.askopenfilename(filetypes = (("Plain text", "*.txt"), ("IDs or loci", "*.ids"), ("All", "*.*")), title = "Please select a file with protein IDs or locuses")
        if ids_filename == "": # Cancel
            return

        self.project_data_tab.protein_ids.read_from_file(ids_filename) # Reading file to a <protein_ids> widget
        ids_list = self.project_data_tab.protein_ids.get_content_as_list() # Obtaining list of protein ids
        self.project_data_tab.print_to_console("Obtained %i protein IDs/locuses. Slicing of the COGNAT database started" % len(ids_list))
        self.set_status("WORKING", "Slicing of the COGNAT database started...", "red")
        database_slice_path = os.path.join(new_project_path, "COGNAT_database")
        os.mkdir(database_slice_path)
        gbff_num = pyCOGNAT_basic.obtain_database_slice(self.settings.cognat_database, database_slice_path, self.settings.script_dir, ids_list, max_gene_number, self.project_data_tab.initital_progress, self.project_data_tab.reading_progress, self.project_data_tab.loading_progress, self.parent)
        if gbff_num == None:
            self.project_data_tab.print_to_console("Slicing of the COGNAT database failed")
            self.set_status("ERROR", "Slicing of the COGNAT database failed", "red")
        else:
            self.project_data_tab.fix_status(ids_filename, self.settings.cognat_database)
            self.project_data_tab.enable_analysis()
            self.project_data_tab.print_to_console("Slicing of the COGNAT database finished, obtained %i assemblies" % gbff_num)
            self.save_project()
            self.set_status("DONE", "Slicing of the COGNAT database finished!", "green")

    def load_project(self, project_dir = None):
        if project_dir == None:
            project_dir = tkinter.filedialog.askdirectory(initialdir = self.settings.work_dir, title = "Please select the directory to load a project")
            if project_dir == "": # Cancel
                return
        if not os.path.isdir(os.path.join(project_dir, "COGNAT_database")): # no COGNAT database is found in this directory
            print ("Given directory '%s' is not a pyCOGNAT project directory" % project_dir)
            self.set_status("ERROR", "Wrong directory was provided, see the console for more details", "red")
            return
        try:
            self.project_data_tab.clear()
            self.project_data_tab.load(project_dir)

            for tab in self.clade_tabs:
                self.tabs.forget(tab)
            self.clade_tabs = list()
            tabs_to_load = list()
            project_files = os.listdir(project_dir)
            for f in project_files:
                match = re.match("([\w ]+)\.([\w ]+)\.ids$", f)
                if match != None:
                    clade_name = match.group(2)
                    color_filename = "%s.%s.colors" % (match.group(1), match.group(2))
                    if color_filename in project_files:
                        color_filename = os.path.join(project_dir, color_filename)
                    else:
                        color_filename = None
                    tabs_to_load.append([clade_name, os.path.join(project_dir, f), color_filename])
            i = 0
            for (tab_name, tab_ids_filename, tab_color_filename) in tabs_to_load:
                new_clade_frame = CladeFrame.CladeFrame(self.tabs, self)
                new_clade_frame.load(tab_name, tab_ids_filename, tab_color_filename)
                self.tabs.insert(i, new_clade_frame, text = tab_name)
                self.clade_tabs.insert(i, new_clade_frame)
                i += 1
            self.set_status("DONE", "Project successfully loaded!", "green")
        except Exception as e:
            self.set_status("ERROR", "Failed to load the project, see the console for more details", "red")
            raise e

    def save_project(self):
        curr_project_name = self.project_data_tab.get_project_name()
        if curr_project_name == "":
            self.set_status("ERROR", "Enter project name before saving one", "red")
            return
        try:
            self.project_data_tab.save()
            project_dir = os.path.join(self.settings.work_dir, curr_project_name)
            for filename in os.listdir(project_dir):
                if re.match("%s\.\w+\.ids", filename): # old clade file
                    file_path = os.path.join(project_dir, filename)
                    os.remove(file_path)
            for clade_tab in self.clade_tabs:
                clade_tab.save()
            self.set_status("DONE", "Project successfully saved!", "green")
        except Exception as e:
            self.set_status("ERROR", "Failed to save the project, see the console for more details", "red")
            print ("Double-check that project name and/or clade names contain only alpha-numeric characters (A-Za-z0-9) and '_'")
            #print ("#, %, &, {, }, \\, <, >, *, ?, /, $, !, ', \", :, @, +, `, |, =")
            raise e

    def add_clade_tab(self):
        clade_frame = CladeFrame.CladeFrame(self.tabs, self)
        self.tabs.insert(0, clade_frame, text = "New clade")
        self.clade_tabs.insert(0, clade_frame)
        self.tabs.select(0)

    def close_app(self):
        tkinter.Tk.destroy(self.parent)

class MainMenubar(tkinter.Menu):
    def __init__(self, parent, host):
        tkinter.Menu.__init__(self, parent, tearoff = 0)
        self.host = host
        self.create_UI()

    def create_UI(self):
        menu_file = tkinter.Menu(self, tearoff = 0)
        menu_file.add_command(label = "New project", command = self.host.create_new_project)
        menu_file.add_command(label = "Load project", command = self.host.load_project)
        menu_file.add_command(label = "Save project", command = self.host.save_project)
        menu_file.add_separator()
        menu_file.add_command(label = "Exit", command = self.host.close_app)
        self.add_cascade(menu = menu_file, label = "File")

        menu_database = tkinter.Menu(self, tearoff = 0)
        menu_database.add_command(label = "Select database", command = self.select_database)
        menu_database.add_command(label = "Select taxonomy color file", command = self.select_tax_color)
        self.add_cascade(menu = menu_database, label = "Database")

        menu_help = tkinter.Menu(self, tearoff = 0)
        menu_help.add_command(label = "Help", command = self.show_help)
        menu_help.add_command(label = "About", command = self.show_about)
        self.add_cascade(menu = menu_help, label = "Help")

    def create_help_window(self, title, h, w, filename):
        help_win = tkinter.Toplevel()
        help_win.title(title)
        help_win.iconbitmap(ICON_FILENAME)
        help_win.wm_attributes("-topmost", 1)

        help_win.grid_rowconfigure(0, weight = 1)
        help_win.grid_columnconfigure(0, weight = 1)
        helptext = tkinter.Text(help_win, height = h, width = w, font = ("Courier New", 10))
        helptext.grid(row = 0, column = 0, sticky = "NSEW")
        text_scr_y = tkinter.Scrollbar(help_win, command = helptext.yview)
        helptext.configure(yscrollcommand = text_scr_y.set)
        text_scr_y.grid(row = 0, column = 1, sticky="NS")
        text_of_help = "* This section is under construction *\nSoon a lot of helpfull information will be added. But not yet."
        helptext.insert(tkinter.END, text_of_help)
        for header in re.findall("\*[\w ]+\*", text_of_help):
            hit_position = helptext.search(header, "1.0", stopindex = tkinter.END)
            end_position = "%s+%dc" % (hit_position, len(header))
            helptext.tag_add("header", hit_position, end_position)
        helptext.tag_config("header", foreground = self.host.header)
        helptext.configure(state = tkinter.DISABLED)

    def show_help(self):
        self.create_help_window("pyCOGNAT main help topics", 35, 150, "dasha.txt")

    def show_about(self):
        logo_win = tkinter.Toplevel()
        logo_win.title("About pyCOGNAT (version %s)" % curr_version)
        logo_win.iconbitmap(ICON_FILENAME)
        logo_win.wm_attributes("-topmost", 1)
        logo_win.grid_rowconfigure(0, weight = 1)
        logo_win.grid_columnconfigure(0, weight = 1)

        logo_image = tkinter.PhotoImage(file = LOGO_FILENAME)
        logo_image = logo_image.subsample(10, 10)
        logo = tkinter.Label(logo_win, image = logo_image)
        logo.image = logo_image
        logo.grid(row = 0, column = 0, rowspan = 3, sticky = "NSEW", padx = self.host.p, pady = self.host.p)

        label1 = tkinter.Label(logo_win, text = "pyCOGNAT is a software for gene neighborhood analysis.", foreground = "#000000", font = ("Arial", 12), anchor = "w", justify = tkinter.LEFT)
        label1.grid(row = 0, column = 1, sticky = "NSW", padx = self.host.p, pady = self.host.p)
        label2 = tkinter.Label(logo_win, text = "Please visit http://boabio.belozersky.msu.ru/ for more information.", foreground = "#000000", font = ("Arial", 12), anchor = "w", justify = tkinter.LEFT)
        label2.grid(row = 1, column = 1, sticky = "NSW", padx = self.host.p, pady = self.host.p)
        label3 = tkinter.Label(logo_win, text = "(c) Dr. Daria V. Dibrova, Lomonosov MSU", foreground = "#000000", font = ("Arial", 12), anchor = "w", justify = tkinter.LEFT)
        label3.grid(row = 2, column = 1, sticky = "NSW", padx = self.host.p, pady = self.host.p)

        logo_win.resizable(False, False)

    def select_tax_color(self):
        filename = tkinter.filedialog.askopenfilename(filetypes = (("Plain text", "*.txt"), ("All", "*.*")), title = "Please select a tab-separated, two-column or four-column taxonomy color file")
        if filename == "": # Cancel
            return
        self.host.settings.tax_colors = filename
        (self.host.gbff_colors, self.host.gbff_order) = pyCOGNAT_basic.read_taxonomy_colors(filename)

    def select_database(self):
        database = tkinter.filedialog.askdirectory(initialdir = self.host.settings.cognat_database, title = "Please select COGNAT database folder")
        if database == "": # Cancel
            return
        self.host.settings.cognat_database = database

root = tkinter.Tk()
root.title("pyCOGNAT (version %s)" % curr_version)
root.iconbitmap(ICON_FILENAME)
root.grid_rowconfigure(0, weight = 1)
root.grid_columnconfigure(0, weight = 1)

def set_proper_size(main_window):
    main_window.update_idletasks()
    w = main_window.winfo_width()
    h = main_window.winfo_height()
    screen_w = main_window.parent.winfo_screenwidth()
    screen_h = main_window.parent.winfo_screenheight()
    main_window.parent.geometry("%dx%d+%d+%d" % (w, screen_h - 100, (screen_w - w)/2, 0))

def rus_copy(event):
    #print ("Crtl + C (russian) pressed!")
    if type(event.widget) == tkinter.Text:
        if event.widget.tag_ranges("sel"):
            content = event.widget.get(tkinter.SEL_FIRST, tkinter.SEL_LAST)
            root.clipboard_clear()
            root.clipboard_append(content)

def rus_paste(event):
    #print ("Crtl + V (russian) pressed!")
    if type(event.widget) == tkinter.Text:
        if event.widget.tag_ranges("sel"):
            event.widget.delete(tkinter.SEL_FIRST, tkinter.SEL_LAST)
        event.widget.insert("insert", event.widget.selection_get(selection='CLIPBOARD'))

def rus_cut(event):
    #print ("Crtl + X (russian) pressed!")
    if type(event.widget) == tkinter.Text:
        if event.widget.tag_ranges("sel"):
            rus_copy(event)
            event.widget.delete(tkinter.SEL_FIRST, tkinter.SEL_LAST)

def select_all(event):
    if type(event.widget) == tkinter.Text:
        event.widget.tag_add("sel", 1.0, "%s-%dc" % (tkinter.END, 1))

def keypress(e):
    if e.keycode == 86 and e.keysym != 'v':
        rus_paste(e)
    elif e.keycode == 67 and e.keysym != 'c':
        rus_copy(e)
    elif e.keycode == 88 and e.keysym != 'x':
        rus_cut(e)
    elif e.keycode == 65 and e.keysym != 'a':
        select_all(e)
    #print ("Pressed: '%s', %s" % (e.keycode, e.keysym))

root.bind("<Control-KeyPress>", keypress)
main_window = pyCOGNAT(root, INI_FILENAME)
main_window.grid(row = 0, column = 0, sticky = "NSEW")
menubar = MainMenubar(root, main_window)
main_window.menubar = menubar
root["menu"] = menubar

set_proper_size(main_window)

if platform.system() == "Windows":
    root.wm_state("zoomed")
root.mainloop()
#!/usr/bin/env python

"""
    pygenefinder GUI.
    Created Nov 2019
    Copyright (C) Damien Farrell

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 3
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""

from __future__ import absolute_import, print_function
import sys,os,subprocess,glob,platform
import pickle
import threading,time
try:
    from tkinter import *
    from tkinter.ttk import *
except:
    from Tkinter import *
    from ttk import *
if (sys.version_info > (3, 0)):
    from tkinter import filedialog, messagebox, simpledialog
    from tkinter.scrolledtext import ScrolledText
else:
    import tkFileDialog as filedialog
    import tkSimpleDialog as simpledialog
    import tkMessageBox as messagebox

import matplotlib
matplotlib.use('TkAgg', warn=False)
import pandas as pd
import numpy as np
from pandastable import Table
from . import tools, app, images, dialogs, tables

home = os.path.expanduser("~")
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module

class Progress(Frame):
    def __init__(self, parent):
        Frame.__init__(self, height=30)
        self.main = self.master
        self.progressbar = Progressbar(parent, orient=HORIZONTAL,
                                       mode='indeterminate', length=200)
        self.progressbar.pack(side=LEFT)

    def start(self):
        self.t = threading.Thread()
        self.t.__init__(target = self.progressbar.start, args = ())
        self.t.start()
        return

    def end(self):
        if self.t.isAlive() == False:
            self.progressbar.stop()
            self.t.join()

class AMRFinderApp(Frame):
    """Application using tkinter widgets.
        Args:
            parent: parent tkinter Frame, default None
    """

    def __init__(self, parent=None, filenames=[], project=None):
        """Initialize the application. """

        self.parent=parent
        if not self.parent:
            Frame.__init__(self)
            self.main=self.master
        else:
            self.main=Toplevel()
            self.master=self.main

        icon = os.path.join(module_path,'logo.gif')
        img = PhotoImage(file=icon)
        self.main.tk.call('wm', 'iconphoto', self.main._w, img)

        self.main.title('pygenefinder')
        self.create_menu_bar()
        self.filenames = filenames
        self.inputs = []
        self.annotations = {}
        self.outputdir = os.path.join(home,'amr_results')
        self.sheets = {}
        self.setup_gui()
        self.main.protocol('WM_DELETE_WINDOW',self.quit)
        self.main.lift()
        if project != None:
            self.load_project(project)
        self.load_fasta_table()
        return

    def setup_gui(self):
        """Add all GUI elements"""

        #progress and info pane
        self.set_geometry()
        bottom = Frame(self.main)
        bottom.pack(side=BOTTOM,fill=BOTH)
        self.progress = Progress(bottom)
        self.progress.pack(side=RIGHT,pady=10,padx=10)
        self.outdirvar = StringVar()
        self.outdirvar.set(self.outputdir)
        outlabel = Label(bottom, textvariable=self.outdirvar)
        outlabel.pack(side=LEFT,pady=10,padx=10)

        self.m = PanedWindow(self.main, orient=HORIZONTAL)
        self.m.pack(fill=BOTH,expand=1)
        fr = self.options_frame()
        self.m.add(fr)

        left = PanedWindow(self.main, orient=VERTICAL)
        #top table
        f1 = Frame(left)
        left.add(f1)
        t = self.fasta_table = tables.FilesTable(f1, app=self, dafaframe=None,
                                                    showtoolbar=0, showstatusbar=1)
        t.model.df = pd.DataFrame()
        t.show()

        #bottom table
        self.nb = Notebook(self.main)
        left.add(self.nb)
        self.m.add(left)

        #right panedwindow
        right = PanedWindow(self.main, orient=VERTICAL)
        self.st = self.text = ScrolledText(right, bg='white', fg='black')
        right.add(self.st)
        self.plotfr = Frame(right)
        self.fig, self.canvas = dialogs.addFigure(self.plotfr)
        self.ax = self.fig.add_subplot(111)
        right.add(self.plotfr)
        self.m.add(right)
        return

    def options_frame(self):

        self.opts = AppOptions(parent=self)
        w = self.opts.showDialog(self.main, layout='vertical')
        dialogs.addButton(w, 'Run', self.run)
        return w

    def create_menu_bar(self):
        """Create the menu bar for the application. """

        self.menu = Menu(self.main)
        file_menu = Menu(self.menu,tearoff=0)
        filemenuitems = {'01Load Fasta Files':{'cmd': lambda: self.load_fasta_files()},
                         '02Load Test Files':{'cmd': lambda: self.load_test()},
                         '03Set Output Folder':{'cmd': lambda: self.set_output_folder()},
                         '04sep':'',
                         '05Load Project':{'cmd': lambda: self.load_project_dialog()},
                         '06Save Project':{'cmd': lambda: self.save_project()},
                         '07Quit':{'cmd':self.quit}}
        self.file_menu = self.create_pulldown(self.menu, filemenuitems, var=file_menu)
        self.menu.add_cascade(label='File',menu=self.file_menu['var'])

        analysis_menu = Menu(self.menu,tearoff=0)
        analysismenuitems = {'01Run':{'cmd': self.run},
                             '02Get sequences for gene':{'cmd': self.show_fasta_sequences},
                             '03Show alignment for gene':{'cmd': self.show_gene_alignment},
                             '04Annotate contigs':{'cmd': self.annotate_all},
                            }
        self.analysis_menu = self.create_pulldown(self.menu, analysismenuitems, var=analysis_menu)
        self.menu.add_cascade(label='Analysis',menu=self.analysis_menu['var'])

        tables_menu = Menu(self.menu,tearoff=0)
        tablesmenuitems = {'01Delete current table':{'cmd': self.delete_table},
                            }
        self.tables_menu = self.create_pulldown(self.menu, tablesmenuitems, var=tables_menu)
        self.menu.add_cascade(label='Tables',menu=self.tables_menu['var'])

        self.help_menu={'01Online Help':{'cmd':self.online_documentation},
                        '02About':{'cmd':self.about}}
        self.help_menu=self.create_pulldown(self.menu,self.help_menu)
        self.menu.add_cascade(label='Help',menu=self.help_menu['var'])

        self.main.config(menu=self.menu)
        return

    def save_project(self, filename='test.pygf'):
        """Save project"""

        data={}
        data['inputs'] = self.inputs
        data['sheets'] = self.sheets
        #data['meta'] = self.saveMeta(table)
        pickle.dump(data, open(filename,'wb'))
        return

    def load_project(self, filename):
        """Load project"""

        #filename='test.pygfr'
        self.clear_project()
        data = pickle.load(open(filename,'rb'))
        self.inputs = data['inputs']
        if 'sheets' in data:
            self.sheets = data['sheets']
        for s in self.sheets:
            self.add_table(s, self.sheets[s])
        ft=self.fasta_table
        ft.model.df = self.inputs
        #ft.expandColumns(150)
        self.fasta_table.redraw()
        self.filename = filename
        self.main.title('pygenefinder: %s' %filename)
        return

    def load_project_dialog(self):

        filename = filedialog.askopenfilename(defaultextension='.dexpl"',
                                                initialdir=os.getcwd(),
                                                filetypes=[("project","*.pygf"),
                                                           ("All files","*.*")],
                                                parent=self.main)
        if not filename:
            return
        if not os.path.exists(filename):
            print ('no such file')
        self.load_project(filename)
        return

    def set_output_folder(self):

        dir = filedialog.askdirectory(initialdir=home,
                                      parent=self.main)
        if not dir:
            return
        self.outputdir = dir
        self.outdirvar.set(dir)
        return

    def load_test(self):
        """Load test_files"""

        files = glob.glob(os.path.join(app.datadir, '*.fa'))
        self.filenames = files
        self.clear_project()
        self.load_fasta_table()
        return

    def load_fasta_files(self):
        """Load fasta files"""

        filenames = filedialog.askopenfilenames(defaultextension='.dexpl"',
                                                initialdir=os.getcwd(),
                                                filetypes=[("fasta","*.fa"),("fasta","*.fasta"),
                                                           ("fasta","*.fna"),
                                                           ("All files","*.*")],
                                                parent=self.main)
        if not filenames:
            return
        self.filenames = filenames
        self.load_fasta_table()
        self.clear_project()
        return

    def clear_project(self):
        """Clear all loaded inputs and results"""

        self.inputs=[]
        self.fasta_table.model.df = pd.DataFrame()
        self.fasta_table.redraw()
        self.sheets={}
        for n in self.nb.tabs():
            self.nb.forget(n)
        self.fig.clear()
        self.canvas.draw()
        return

    def add_table(self, name, df, kind='results'):
        """Add a table"""

        self.sheets[name] = df
        self.show_table(name, df, kind)
        return

    def show_table(self, name, df, kind='results'):
        """Add a table to the results notebook"""

        f = Frame(height=500)
        if kind == 'results':
            t = tables.GenesTable(f, dataframe=df, app=self, showtoolbar=0, showstatusbar=1,
                                    editable=False, height=500)
        elif kind == 'features':
            t = tables.FeaturesTable(f, dataframe=df, app=self, showtoolbar=0, showstatusbar=1,
                                    editable=False, height=500)
        t.show()
        self.nb.add(f, text=name)
        return

    def delete_table(self):
        """Delete current table"""

        s = self.nb.index(self.nb.select())
        name = self.nb.tab(s, 'text')
        self.nb.forget(s)
        del self.sheets[name]
        return

    def get_selected_table(self):
        """current table"""

        s = self.nb.index(self.nb.select())
        name = self.nb.tab(s, 'text')
        t = self.sheets[name]
        return t

    def load_fasta_table(self):
        """Load fasta inputs into table"""

        if self.filenames is None or len(self.filenames) == 0:
            return
        names = [os.path.splitext(os.path.basename(i))[0] for i in self.filenames]
        self.inputs = pd.DataFrame({'label':names,'filename':self.filenames})
        pt = self.fasta_table
        pt.model.df = self.inputs
        #pt.expandColumns(50)
        pt.redraw()
        return

    def write(self, string):
        """for stdout redirect to scrolledtext"""

        self.st.insert(END, string)
        self.st.see(END)

    def flush(self, string=None):
        return

    def get_selected_file(self):
        """Get selected file from files table"""

        df = self.fasta_table.model.df
        row = self.fasta_table.getSelectedRow()
        data = df.iloc[row]
        return data

    def show_fasta(self):
        """SHow selected input fasta file"""

        df = self.fasta_table.model.df
        row = self.fasta_table.getSelectedRow()
        data = df.iloc[row]
        from Bio import SeqIO
        seqs = SeqIO.parse(data.filename, 'fasta')
        w = Toplevel(self)
        w.grab_set()
        w.transient(self)
        w.title(data.filename)
        ed = dialogs.SimpleEditor(w, height=25)
        ed.pack(in_=w, fill=BOTH, expand=Y)

        for s in seqs:
            ed.text.insert(END, s.format("fasta"))
        return

    def show_annotation(self):
        """Show annotation results for an input file in table"""

        row = self.get_selected_file()
        name = row.label

        if name in self.sheets:
            return
        genbank = row.genbank
        print(row.genbank)
        featsdf = tools.genbank_to_dataframe(row.genbank)

        self.add_table(name, featsdf, kind='features')
        return

    def annotate_file(self):

        self.opts.applyOptions()
        kwds = self.opts.kwds
        df = self.fasta_table.model.df
        row = self.fasta_table.getSelectedRow()
        data = df.iloc[row]
        name = data.label
        if not os.path.exists(self.outputdir):
            os.makedirs(self.outputdir)
        outfile = os.path.join(self.outputdir, data.label + '.gbk')
        self.progress.start()
        def fun1(a):
            for i in range(5):
                print (i)
                time.sleep(1)
            return
        thread = threading.Thread(target = fun1, kwargs={'a':2})
        #thread.start()
        print (outfile)
        #feats = app.annotate_contigs(data.filename, outfile, threads=kwds['threads'])

        #self.progressbar.stop()
        df.loc[row,'genbank'] = outfile
        self.fasta_table.redraw()
        return

    def annotate_all(self):
        """Run gene annotation for input files"""

        self.opts.applyOptions()
        kwds = self.opts.kwds
        self.inputs = self.fasta_table.model.df
        files = self.inputs.filename
        self.annotations = {}
        for file in files:
            #app.annotate_contigs(file, threads=kwds['threads'])
            self.annotate_file(file, **self.opts.kwds)
            sys.stdout.flush()

        return

    def run(self):
        """Run pipeline"""

        #self.progress.start()
        sys.stdout = self
        if len(self.inputs) == 0:
            print('you need to load fasta files')
            return
        self.opts.applyOptions()
        app.check_databases()
        kwds = self.opts.kwds
        db = kwds['db']
        app.fetch_sequence_from_url(db)
        #update inputs from table
        self.inputs = self.fasta_table.model.df
        files = self.inputs.filename
        print ('running %s files' %len(files))
        app.make_target_database(files)
        targfile = os.path.join(app.tempdir, 'targets.fasta')
        bl = app.find_genes(targfile, db, **kwds)
        #bl.to_csv('%s_results.csv' %db)
        print ('found %s genes' %len(bl.gene.unique()))
        m = app.pivot_blast_results(bl)
        self.bl = bl
        #print (m)
        self.fig.clear()
        app.plot_heatmap(m.T, fig=self.fig, title='results matrix')
        self.canvas.draw()
        #m.to_csv('%s_matrix.csv' %db)
        print ('done')

        self.show_table(db, bl)
        self.show_table('matrix',m.T.reset_index())
        return

    def get_selected_gene(self):
        """get selected gene"""

        t = self.get_selected_table()
        row = t.getSelectedRow()
        data = self.bl.iloc[row]
        gene = data.gene
        return gene

    def show_fasta_sequences(self):

        files = self.inputs.filename
        gene = self.get_selected_gene()
        seqs = app.get_gene_hits(self.bl, gene, files, db='card')
        print()
        for s in seqs:
            print(s.format("fasta"))
        return

    def show_gene_alignment(self):

        files = self.inputs.filename
        gene = self.get_selected_gene()
        seqs = app.get_gene_hits(self.bl, gene, files, db='card')
        #maaft_alignment(seqfile)
        print ('alignments for gene: %s' %gene)
        app.get_alignment(seqs)
        return

### utility methods

    def get_best_geometry(self):
        """Calculate optimal geometry from screen size"""

        ws = self.main.winfo_screenwidth()
        hs = self.main.winfo_screenheight()
        self.w = w = ws/1.4; h = hs*0.7
        x = (ws/2)-(w/2); y = (hs/2)-(h/2)
        g = '%dx%d+%d+%d' % (w,h,x,y)
        return g

    def set_geometry(self):
        self.winsize = self.get_best_geometry()
        self.main.geometry(self.winsize)
        return

    def create_pulldown(self, menu, dict, var=None):
        """Create pulldown menu, returns a dict.
        Args:
            menu: parent menu bar
            dict: dictionary of the form -
            {'01item name':{'cmd':function name, 'sc': shortcut key}}
            var: an already created menu
        """

        if var is None:
            var = Menu(menu,tearoff=0)
        #dialogs.applyStyle(var)
        items = list(dict.keys())
        items.sort()
        for item in items:
            if item[-3:] == 'sep':
                var.add_separator()
            else:
                command = dict[item]['cmd']
                label = '%-25s' %(item[2:])
                if 'img' in dict[item]:
                    img = dict[item]['img']
                else:
                    img = None
                if 'sc' in dict[item]:
                    sc = dict[item]['sc']
                    #bind command
                    #self.main.bind(sc, command)
                else:
                    sc = None
                var.add('command', label=label, command=command, image=img,
                        compound="left")#, accelerator=sc)
        dict['var'] = var
        return dict

    def set_styles(self):
        """Set theme and widget styles"""

        style = self.style = Style(self)
        available_themes = self.style.theme_names()
        plf = util.checkOS()
        if plf == 'linux':
            style.theme_use('default')
        elif plf == 'darwin':
            style.theme_use('clam')

        self.bg = bg = self.style.lookup('TLabel.label', 'background')
        style.configure('Horizontal.TScale', background=bg)
        #set common background style for all widgets because of color issues
        #if plf in ['linux','darwin']:
        #    self.option_add("*background", bg)
        dialogs.applyStyle(self.menu)
        return

    def _check_snap(self):
        if os.environ.has_key('SNAP_USER_COMMON'):
            print ('running inside snap')
            return True
        return False

    def about(self):
        """About dialog"""

        abwin = Toplevel()
        x,y,w,h = dialogs.get_parent_geometry(self.main)
        abwin.geometry('+%d+%d' %(x+w/2-200,y+h/2-200))
        abwin.title('About')
        abwin.transient(self)
        abwin.grab_set()
        abwin.resizable(width=False, height=False)
        #abwin.configure(background=self.bg)
        logo = images.logo_image()
        label = Label(abwin,image=logo,anchor=CENTER)
        label.image = logo
        label.grid(row=0,column=0,sticky='ew',padx=4,pady=4)
        style = Style()
        style.configure("BW.TLabel", font='arial 11')
        from . import __version__
        pandasver = pd.__version__
        pythonver = platform.python_version()
        mplver = matplotlib.__version__
        if self._check_snap == True:
            snap='(snap)'
        else:
            snap=''

        text='pygenefinder GUI\n'\
                +'version '+__version__+snap+'\n'\
                +'Copyright (C) Damien Farrell 2019-\n'\
                +'This program is free software; you can redistribute it and/or\n'\
                +'modify it under the terms of the GNU General Public License\n'\
                +'as published by the Free Software Foundation; either version 3\n'\
                +'of the License, or (at your option) any later version.\n'\
                +'Using Python v%s\n' %pythonver\
                +'pandas v%s, matplotlib v%s' %(pandasver,mplver)

        row=1
        #for line in text:
        tmp = Label(abwin, text=text, style="BW.TLabel")
        tmp.grid(row=row,column=0,sticky='news',pady=2,padx=4)
        return

    def online_documentation(self,event=None):
        """Open the online documentation"""
        import webbrowser
        link='https://github.com/dmnfarrell/pygenefinder'
        webbrowser.open(link,autoraise=1)
        return

class AppOptions(dialogs.TkOptions):
    """Class for provinding options"""
    def __init__(self, parent=None):
        """Setup variables"""

        self.parent = parent
        import multiprocessing
        cpus = multiprocessing.cpu_count()
        threadlist = list(range(1,cpus+1))
        dbs = app.db_names
        self.groups = {'options':['db','identity','coverage','threads']}
        self.opts = {'db':{'type':'combobox','default':'card',
                    'items':dbs,'label':'database'},
                    'identity':{'type':'entry','default':90},
                    'coverage':{'type':'entry','default':50},
                    'threads':{'type':'entry','default':4},
                    #'best hit':{'type':'checkbutton','default':True,
                    #'label':'keep best hit only'}
                    }
        return

def main():
    "Run the application"

    import sys, os
    from argparse import ArgumentParser
    parser = ArgumentParser(description='AMRfinder tool')
    parser.add_argument("-f", "--fasta", dest="filenames",default=[],
                        help="input fasta file", metavar="FILE")
    parser.add_argument("-p", "--proj", dest="project",default=None,
                        help="load .pygf project file", metavar="FILE")
    args = vars(parser.parse_args())
    #if args['filenames'] != None:
    #    app = AMRFinderApp(filenames=args['filenames'])
    #elif args['project'] != None:
    app = AMRFinderApp(**args)

    app.mainloop()

if __name__ == '__main__':
    main()

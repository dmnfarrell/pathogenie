#!/usr/bin/env python

"""
    pyamrfinder GUI.
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
from . import tools, app, images, dialogs

home = os.path.expanduser("~")
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
datadir = os.path.join(module_path, 'data')


class AMRFinderApp(Frame):
    """Application using tkinter widgets.
        Args:
            parent: parent tkinter Frame, default None
    """

    def __init__(self, parent=None, filenames=[]):
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

        self.main.title('pyAMRfinder')
        self.create_menu_bar()
        self.setup_gui()
        self.filenames = filenames
        self.inputs = []
        self.main.protocol('WM_DELETE_WINDOW',self.quit)
        self.main.lift()
        self.load_fasta_table()
        return

    def setup_gui(self):
        """Add all GUI elements"""

        self.m = PanedWindow(self.main, orient=HORIZONTAL)
        self.m.pack(fill=BOTH,expand=1)
        fr = self.options_frame()
        self.m.add(fr)

        left = PanedWindow(self.main, orient=VERTICAL)
        #top table
        f1 = Frame(left)
        left.add(f1)
        t = self.fasta_table = Table(f1, dafaframe=None, showtoolbar=0, showstatusbar=1, editable=False)
        t.model.df = pd.DataFrame()
        t.show()

        #bottom table
        self.nb = Notebook(self.main)
        f2 = Frame()
        t = self.results_table = Table(f2, showtoolbar=0, showstatusbar=1, height=600)
        t.model.df = pd.DataFrame()
        t.show()
        self.nb.add(f2, text='results')
        f3 = Frame()
        t = self.matrix_table = Table(f3, showtoolbar=0, showstatusbar=1, height=600)
        t.model.df = pd.DataFrame()
        t.show()
        self.nb.add(f3, text='matrix')

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
        self.set_geometry()
        return

    def options_frame(self):

        f = Frame(self.main)
        self.opts = AppOptions(parent=self)
        w = self.opts.showDialog(self.main, layout='vertical')
        return w

    def create_menu_bar(self):
        """Create the menu bar for the application. """

        self.menu = Menu(self.main)
        file_menu = Menu(self.menu,tearoff=0)
        filemenuitems = {'01Load Fasta':{'cmd': lambda: self.load_fasta_files()},
                         '02Load Test Files':{'cmd': lambda: self.load_test()},
                    '04sep':'',
                    '06Quit':{'cmd':self.quit}}
        self.file_menu = self.create_pulldown(self.menu, filemenuitems, var=file_menu)
        self.menu.add_cascade(label='File',menu=self.file_menu['var'])

        analysis_menu = Menu(self.menu,tearoff=0)
        analysismenuitems = {'01Run':{'cmd': self.run}}
        self.analysis_menu = self.create_pulldown(self.menu, analysismenuitems, var=analysis_menu)
        self.menu.add_cascade(label='Analysis',menu=self.analysis_menu['var'])

        self.help_menu={'01Online Help':{'cmd':self.online_documentation},
                        '02About':{'cmd':self.about}}
        self.help_menu=self.create_pulldown(self.menu,self.help_menu)
        self.menu.add_cascade(label='Help',menu=self.help_menu['var'])

        self.main.config(menu=self.menu)
        return

    def read_file_info(self):

        return

    def load_test(self):
        """Load test_files"""

        files = glob.glob(os.path.join(datadir, '*.fa'))
        self.filenames = files
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
        return

    def load_fasta_table(self):
        """Load fasta inputs into table"""

        if self.filenames is None or len(self.filenames) == 0:
            return
        names = [os.path.basename(i) for i in self.filenames]
        self.inputs = pd.DataFrame({'label':names,'filename':self.filenames})
        pt = self.fasta_table
        pt.model.df = self.inputs
        pt.adjustColumnWidths()
        pt.redraw()
        return

    def load_results(self):
        #print (self.bl)
        self.results_table.model.df = self.bl
        self.results_table.redraw()
        return

    def write(self, string):
        """for stdout redirect to scrolledtext"""

        self.st.insert(END, string)

    def flush(self, string):
        return

    def run(self):
        """Run pipeline"""

        sys.stdout = self
        #app.run(self.filenames, 'card')
        if len(self.inputs) == 0:
            print('you need to load fasta files')
            return
        self.opts.applyOptions()
        print (self.opts.kwds)
        db = self.opts.kwds['db']
        app.fetch_sequence_db(db)
        #update inputs from table
        self.inputs = self.fasta_table.model.df
        files = self.inputs.filename
        print ('running %s files' %len(files))
        app.make_blast_database(files)
        bl = app.run_blast('out.fasta', db)
        bl.to_csv('%s_results.csv' %db)
        m = app.pivot_blast_results(bl)
        self.bl = bl
        print (m)
        app.plot_heatmap(m, ax=self.ax)
        self.canvas.draw()
        t=self.matrix_table
        t.model.df = m.T.reset_index()
        t.redraw()
        t.setWrap()
        m.to_csv('%s_matrix.csv' %db)
        print ('done')
        self.load_results()
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

        text='pyamrfinder GUI\n'\
                +'version '+__version__+snap+'\n'\
                +'Copyright (C) Damien Farrell 2014-\n'\
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
        link='https://github.com/dmnfarrell/pyamrfinder'
        webbrowser.open(link,autoraise=1)
        return

class AppOptions(dialogs.TkOptions):
    """Class for provinding options"""
    def __init__(self, parent=None):
        """Setup variables"""

        self.parent = parent
        dbs = ['card','resfinder','arg-annot','vfdb']
        self.groups = {'options':['db']}
        self.opts = {'db':{'type':'combobox','default':'card',
                    'items':dbs,
                    'label':'database'},}
        return

def main():
    "Run the application"

    import sys, os
    from argparse import ArgumentParser
    parser = ArgumentParser(description='AMRfinder tool')
    parser.add_argument("-f", "--fasta", dest="filenames",default=[],
                        help="input fasta file", metavar="FILE")
    args = vars(parser.parse_args())

    app = AMRFinderApp(filenames=args['filenames'])
    app.mainloop()

if __name__ == '__main__':
    main()

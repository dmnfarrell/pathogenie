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
from . import tools, app, images

home = os.path.expanduser("~")
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
datadir = os.path.join(module_path, 'data')

def get_parent_geometry(parent):
    x = parent.winfo_rootx()
    y = parent.winfo_rooty()
    w = parent.winfo_width()
    h = parent.winfo_height()
    return x,y,w,h

class AMRFinderApp(Frame):
    """Application using tkinter widgets.
        Args:
            parent: parent tkinter Frame, default None
    """

    def __init__(self, parent=None):
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
        self.filenames = []
        self.inputs = []
        self.main.protocol('WM_DELETE_WINDOW',self.quit)
        self.main.lift()
        return

    def setup_gui(self):
        """Add all GUI elements"""

        self.m = PanedWindow(self.main, orient=HORIZONTAL)
        self.m.pack(fill=BOTH,expand=1)
        tables = PanedWindow(self.main, orient=VERTICAL)
        f1 = Frame(tables)
        tables.add(f1)
        t = self.fasta_table = Table(f1, dafaframe=None, showtoolbar=0, showstatusbar=1, editable=False)
        t.model.df = pd.DataFrame()
        t.show()
        f2 = Frame(tables)
        t = self.results_table = Table(f2, showtoolbar=0, showstatusbar=1, height=600)
        t.model.df = pd.DataFrame()
        tables.add(f2)
        t.show()
        self.m.add(tables)
        self.st = self.text = ScrolledText(self.main, bg='white', fg='black')
        self.m.add(self.st)
        self.set_geometry()
        return

    def create_menu_bar(self):
        """Create the menu bar for the application. """

        self.menu = Menu(self.main)
        file_menu = Menu(self.menu,tearoff=0)
        filemenuitems = {'01New Project':{'cmd': self.new_project},
                    '02Load Fasta':{'cmd': lambda: self.load_fasta_files()},
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

    def new_project(self):
        return

    def close_project(self):
        return

    def load_fasta_files(self):
        """Load fasta files"""
        filenames = filedialog.askopenfilenames(defaultextension='.dexpl"',
                                                initialdir=os.getcwd(),
                                                filetypes=[("fasta","*.fa"),("fasta","*.fasta"),
                                                           ("All files","*.*")],
                                                parent=self.main)
        if not filenames:
            return
        self.filenames = filenames
        self.inputs = pd.DataFrame({'name':self.filenames})
        self.load_fasta_table()
        return

    def write(self, string):
        """"""
        self.st.insert(END, string)

    def run(self):

        sys.stdout = self
        #app.run(self.filenames, 'card')
        if len(self.inputs) == 0:
            print('you need to load fasta files')
            return
        db='card'
        app.fetch_sequence_db(db)
        app.make_blast_database(self.filenames)
        bl = app.run_blast('out.fasta', db)
        bl.to_csv('%s_results.csv' %db)
        m = app.pivot_blast_results(bl)
        self.bl = bl
        print (m)
        app.plot_heatmap(m)
        m.to_csv('%s_matrix.csv' %db)
        print ('done')
        self.load_table()
        return

    def load_fasta_table(self):

        self.fasta_table.model.df = self.inputs
        self.fasta_table.redraw()
        return

    def load_table(self):
        print (self.bl)
        self.results_table.model.df = self.bl
        self.results_table.redraw()
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
        x,y,w,h = get_parent_geometry(self.main)
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

def main():
    "Run the application"

    import sys, os
    from argparse import ArgumentParser
    parser = ArgumentParser(description='AMRfinder tool')
    parser.add_argument("-f", "--fasta", dest="filename",
                        help="input fasta file", metavar="FILE")
    args = vars(parser.parse_args())

    app = guiAMRFinderApp()
    app.mainloop()

if __name__ == '__main__':
    main()

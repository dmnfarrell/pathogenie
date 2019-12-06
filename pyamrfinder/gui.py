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
import sys,os,subprocess,glob
try:
    from tkinter import *
    from tkinter.ttk import *
except:
    from Tkinter import *
    from ttk import *
if (sys.version_info > (3, 0)):
    from tkinter import filedialog, messagebox, simpledialog
else:
    import tkFileDialog as filedialog
    import tkSimpleDialog as simpledialog
    import tkMessageBox as messagebox

import pandas as pd
import numpy as np
from . import tools

home = os.path.expanduser("~")
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
datadir = os.path.join(module_path, 'data')

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

        #icon = os.path.join(module_path,'icon.gif')
        #img = PhotoImage(file=icon)
        #self.main.tk.call('wm', 'iconphoto', self.main._w, img)

        self.main.title('pyAMRfinder')
        self.create_menu_bar()
        self.setup_gui()

        self.main.protocol('WM_DELETE_WINDOW',self.quit)
        self.main.lift()
        return

    def setup_gui(self):
        """Add all GUI elements"""

        self.m = PanedWindow(self.main, orient=HORIZONTAL)
        self.m.pack(fill=BOTH,expand=1)
        self.nb = Notebook(self.main)
        self.m.add(self.nb)
        self.set_geometry()
        return

    def create_menu_bar(self):
        """Create the menu bar for the application. """

        self.menu = Menu(self.main)
        file_menu = Menu(self.menu,tearoff=0)
        filemenuitems = {'01New Project':{'cmd': self.new_project},
                    '02Load Fasta':{'cmd': lambda: self.load_fasta()},
                    '04sep':'',
                    '06Quit':{'cmd':self.quit}}
        self.file_menu = self.create_pulldown(self.menu, filemenuitems, var=file_menu)
        self.menu.add_cascade(label='File',menu=self.file_menu['var'])
        self.main.config(menu=self.menu)

        analysis_menu = Menu(self.menu,tearoff=0)
        analysismenuitems = {'01Run':{'cmd': self.run}}
        self.analysis_menu = self.create_pulldown(self.menu, analysismenuitems, var=analysis_menu)
        self.menu.add_cascade(label='Analysis',menu=self.analysis_menu['var'])
        return

    def new_project(self):
        return

    def close_project(self):
        return

    def load_fasta(self):
        """Load fasta files"""

        return

    def run(self):

        return

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

def main():
    "Run the application"

    import sys, os
    from argparse import ArgumentParser
    parser = ArgumentParser(description='AMRfinder tool')
    parser.add_argument("-f", "--fasta", dest="filename",
                        help="input fasta file", metavar="FILE")
    args = vars(parser.parse_args())

    app = AMRFinderApp()
    app.mainloop()

if __name__ == '__main__':
    main()

#!/usr/bin/env python

"""
    pyamrfinder pandastable sub-classes.
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
import sys,os,platform
try:
    from tkinter import *
    from tkinter.ttk import *
except:
    from Tkinter import *
    from ttk import *

import pandas as pd
import numpy as np
from pandastable import Table

class FilesTable(Table):
    """
      Custom table class inherits from Table.
      You can then override required methods
     """
    def __init__(self, parent=None, app=None, **kwargs):
        Table.__init__(self, parent, **kwargs)
        self.app = app
        return

    def popupMenu(self, event, rows=None, cols=None, outside=None):
        """Custom right click menu"""

        defaultactions = {
            "Preferences" : self.showPreferences,
            }
        popupmenu = Menu(self, tearoff = 0)
        def popupFocusOut(event):
            popupmenu.unpost()
        popupmenu.add_command(label='Show Sequences', command=self.app.show_fasta)
        popupmenu.add_command(label='Annotate File', command=self.app.annotate_file)
        #popupmenu.add_command(label='Annotate all', command=self.app.annotate)
        popupmenu.add_command(label='Show Annotation', command=self.app.show_annotation)
        popupmenu.bind("<FocusOut>", popupFocusOut)
        popupmenu.focus_set()
        popupmenu.post(event.x_root, event.y_root)
        return popupmenu

class GenesTable(Table):
    """
      Custom table class inherits from Table.
      You can then override required methods
     """
    def __init__(self, parent=None, app=None, **kwargs):
        Table.__init__(self, parent, **kwargs)
        self.app = app
        return

    def popupMenu(self, event, rows=None, cols=None, outside=None):
        """Custom right click menu"""

        popupmenu = Menu(self, tearoff = 0)
        def popupFocusOut(event):
            popupmenu.unpost()
        popupmenu.add_command(label='Show gene sequences', command=self.app.show_fasta_sequences)
        popupmenu.add_command(label='Show gene alignment', command=self.app.show_gene_alignment)
        popupmenu.bind("<FocusOut>", popupFocusOut)
        popupmenu.focus_set()
        popupmenu.post(event.x_root, event.y_root)
        return popupmenu

class FeaturesTable(Table):
    """
      Custom table class inherits from Table.
      You can then override required methods
     """
    def __init__(self, parent=None, app=None, **kwargs):
        Table.__init__(self, parent, **kwargs)
        self.app = app
        return

    def popupMenu(self, event, rows=None, cols=None, outside=None):
        """Custom right click menu"""

        popupmenu = Menu(self, tearoff = 0)
        def popupFocusOut(event):
            popupmenu.unpost()
        #popupmenu.add_command(label=' ', command=self.app.func)

        popupmenu.bind("<FocusOut>", popupFocusOut)
        popupmenu.focus_set()
        popupmenu.post(event.x_root, event.y_root)
        return popupmenu

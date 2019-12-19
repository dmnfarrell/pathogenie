#!/usr/bin/env python

"""
    pygenefinder dialog widgets.
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

def get_parent_geometry(parent):
    x = parent.winfo_rootx()
    y = parent.winfo_rooty()
    w = parent.winfo_width()
    h = parent.winfo_height()
    return x,y,w,h

def getDictfromTkVars(opts, tkvars, widgets):
    kwds = {}
    for i in opts:
        if not i in tkvars:
            continue
        if opts[i]['type'] == 'listbox':
            items = widgets[i].curselection()
            kwds[i] = [widgets[i].get(j) for j in items]
            #print (items, kwds[i])
        else:
            try:
                kwds[i] = int(tkvars[i].get())
            except:
                kwds[i] = tkvars[i].get()
    return kwds

def addButton(frame, name, callback, img=None,
              side=TOP, compound=None, width=None, padding=2):
    """Add a button with image, toolip to a tkinter frame"""

    if img==None:
        b = Button(frame, text=name, command=callback)
    else:
        b = Button(frame, text=name, command=callback, width=width,
                         image=img, compound=compound, padding=padding)
    b.image = img
    b.pack(side=side,fill=X,pady=padding)
    return

def dialogFromOptions(parent, opts, groups=None, callback=None,
                        sticky='news',  layout='horizontal'):
    """Auto create tk vars and widgets for corresponding options and
       and return the enclosing frame"""

    tkvars = {}
    widgets = {}
    dialog = Frame(parent)
    if groups == None:
        groups = {'options': opts.keys()}
    c=0
    row=0
    for g in groups:
        if g == 'hidden':
            continue
        if layout=='horizontal':
            row=0; c+=1
            side=LEFT
            fill=Y
        else:
            c=0; row+=1
            side=TOP
            fill=X
        frame = LabelFrame(dialog, text=g)
        #frame.grid(row=row,column=c,sticky=sticky)
        frame.pack(side=side,fill=fill,expand=False)

        for i in groups[g]:
            w=None
            opt = opts[i]
            if opt['type'] == 'entry':
                if 'label' in opt:
                    label=opt['label']
                else:
                    label=i
                if 'width' in opt:
                    w=opt['width']
                else:
                    w=6
                Label(frame,text=label).pack()
                if type(opts[i]['default']) is int:
                    tkvars[i] = v = IntVar()
                else:
                    tkvars[i] = v = StringVar()
                v.set(opts[i]['default'])
                w = Entry(frame,textvariable=v, width=w, command=callback)
            elif opt['type'] == 'scrolledtext':
                w = ScrolledText(frame, width=20, wrap=WORD)
                tkvars[i] = None
            elif opt['type'] == 'checkbutton':
                tkvars[i] = v = IntVar()
                v.set(opts[i]['default'])
                w = Checkbutton(frame,text=opt['label'],
                         variable=v)
            elif opt['type'] == 'combobox':
                if 'label' in opt:
                   label=opt['label']
                else:
                    label = i
                if 'width' in opt:
                    w=opt['width']
                else:
                    w=16
                Label(frame,text=label).pack()
                tkvars[i] = v = StringVar()
                v.set(opts[i]['default'])
                w = Combobox(frame, values=opt['items'],
                         textvariable=v,width=w,
                         validatecommand=callback,validate='key')
                w.set(opt['default'])
                if 'tooltip' in opt:
                    ToolTip.createToolTip(w, opt['tooltip'])
            elif opt['type'] == 'spinbox':
                if 'label' in opt:
                   label=opt['label']
                else:
                    label = i
                Label(frame,text=label).pack()
                tkvars[i] = v = StringVar()
                w = Spinbox(frame, values=opt['items'],
                         textvariable=v,width=w,
                         validatecommand=callback,validate='key')
                w.set(opt['default'])
                #if 'tooltip' in opt:
                #    ToolTip.createToolTip(w, opt['tooltip'])
            elif opt['type'] == 'listbox':
                if 'label' in opt:
                   label=opt['label']
                else:
                    label = i
                Label(frame,text=label).pack()
                w,v = addListBox(frame, values=opt['items'],width=12)
                tkvars[i] = v #add widget instead of var
            elif opt['type'] == 'radio':
                Label(frame,text=label).pack()
                if 'label' in opt:
                   label=opt['label']
                else:
                    label = i
                Label(frame,text=label).pack()
                tkvars[i] = v = StringVar()
                for item in opt['items']:
                    w = Radiobutton(frame, text=item, variable=v, value=item).pack()
            elif opt['type'] == 'scale':
                fr,to=opt['range']
                tkvars[i] = v = DoubleVar()
                v.set(opts[i]['default'])
                w = Scale(frame,label=opt['label'],
                         from_=fr,to=to,
                         orient='horizontal',
                         resolution=opt['interval'],
                         variable=v)
            elif opt['type'] == 'colorchooser':
                tkvars[i] = v = StringVar()
                clr = opts[i]['default']
                v.set(clr)
                def func(var):
                    clr=var.get()
                    new=pickColor(parent,clr)
                    if new != None:
                        var.set(new)
                w = Button(frame, text=opt['label'], command=lambda v=v : func(v))

            if w != None:
                w.pack(fill=BOTH,expand=1,pady=1)
                widgets[i] = w
            row+=1

    return dialog, tkvars, widgets

def applyStyle(w):
    """Apply style to individual widget to prevent widget color issues on linux"""

    plf = util.checkOS()
    if plf in ['linux','darwin']:
        bg = Style().lookup('TLabel.label', 'background')
        w.configure(fg='black', bg=bg,
                     activeforeground='white', activebackground='#0174DF')
    return

def setWidgetStyles(widgets):
    """set styles of list of widgets"""

    style = Style()
    bg = style.lookup('TLabel.label', 'background')
    for w in widgets:
        try:
            w.configure(fg='black', bg=bg)
        except:
            pass
    return

def addFigure(parent, figure=None, resize_callback=None):
    """Create a tk figure and canvas in the parent frame"""

    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg#, NavigationToolbar2TkAgg
    from matplotlib.figure import Figure

    if figure == None:
        figure = Figure(figsize=(5,4), dpi=80, facecolor='white')

    canvas = FigureCanvasTkAgg(figure, master=parent, resize_callback=resize_callback)
    canvas.draw()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
    canvas.get_tk_widget().configure(highlightcolor='gray75',
                                   highlightbackground='gray75')
    #toolbar = NavigationToolbar2TkAgg(canvas, parent)
    #toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    return figure, canvas

def addListBox(parent, values=[], width=10, height=6, label=''):
    """Add an EasyListBox"""

    frame=Frame(parent)
    Label(frame, text=label).grid(row=0)
    yScroll = Scrollbar(frame, orient = VERTICAL)
    yScroll.grid(row = 1, column = 1, sticky = N+S)
    listItemSelected = lambda index: index
    lbx = Listbox(frame, width, height, yScroll.set, listItemSelected)
    lbx.grid(row = 1, column = 0, sticky = N+S+E+W)
    frame.columnconfigure(0, weight = 1)
    frame.rowconfigure(0, weight = 1)
    yScroll["command"] = lbx.yview
    for i in values:
        lbx.insert(END, i)
    return frame, lbx

def getListBoxSelection(w):
    items = w.curselection()
    return [w.get(j) for j in items]

class TkOptions(object):
    """Class to generate tkinter widget dialog for dict of options"""
    def __init__(self, parent=None):
        """Setup variables"""

        self.parent = parent
        return

    def applyOptions(self):
        """Set the plot kwd arguments from the tk variables"""

        kwds = {}
        for i in self.opts:
            if not i in self.tkvars:
                continue
            #print (i, self.opts[i]['type'], self.widgets[i])
            if self.opts[i]['type'] == 'listbox':
                items = self.widgets[i].curselection()
                kwds[i] = [self.widgets[i].get(j) for j in items]
            elif self.opts[i]['type'] == 'scrolledtext':
                kwds[i] = self.widgets[i].get('1.0',END)
            else:
                kwds[i] = self.tkvars[i].get()
        self.kwds = kwds
        return

    def setWidgetStyles(self):

        style = Style()
        bg = style.lookup('TLabel.label', 'background')
        for i in self.widgets:
            try:
                self.widgets[i].configure(fg='black', bg=bg)
            except:
                pass
        return

    def apply(self):
        self.applyOptions()
        if self.callback != None:
            self.callback()
        return

    def showDialog(self, parent, layout='horizontal'):
        """Auto create tk vars, widgets for corresponding options and
           and return the frame"""

        dialog, self.tkvars, self.widgets = dialogFromOptions(parent,
                                                              self.opts, self.groups,
                                                              layout=layout)
        self.setWidgetStyles()
        return dialog

    def updateFromOptions(self, kwds=None):
        """Update all widget tk vars using plot kwds dict"""

        if kwds != None:
            self.kwds = kwds
        elif hasattr(self, 'kwds'):
            kwds = self.kwds
        else:
            return
        if self.tkvars == None:
            return
        for i in kwds:
            if i in self.tkvars and self.tkvars[i]:
                self.tkvars[i].set(kwds[i])
        return

class SimpleEditor(Frame):
    """Simple text editor"""

    def __init__(self, parent=None, width=100, height=40, font=None):

        Frame.__init__(self, parent)
        st = self.text = ScrolledText(self, width=width, height=height, bg='white', fg='black')
        st.pack(in_=self, fill=BOTH, expand=1)
        if font == None:
            if 'Windows' in platform.system():
                font = ('Courier New',10)
            else:
                font = 'monospace 10'
        st.config(font=font)
        btnform = Frame(self)
        btnform.pack(fill=BOTH)
        Button(btnform, text='Save',  command=self.onSave).pack(side=LEFT)
        #Button(frm, text='Cut',   command=self.onCut).pack(side=LEFT)
        #Button(frm, text='Paste', command=self.onPaste).pack(side=LEFT)
        Button(btnform, text='Find',  command=self.onFind).pack(side=LEFT)
        Button(btnform, text='Clear',  command=self.onClear).pack(side=LEFT)
        self.target=''
        return

    def onSave(self):
        """Save text"""

        filename = filedialog.asksaveasfilename(defaultextension='.txt',
                                    initialdir=os.path.expanduser('~'),
                                     filetypes=(('Text files', '*.txt'),
                                                ('All files', '*.*')))
        if filename:
            with open(filename, 'w') as stream:
                stream.write(self.text.get('1.0',END))
        return

    def onClear(self):
        """Clear text"""
        self.text.delete('1.0',END)
        return

    def onFind(self):
        self.target = simpledialog.askstring('SimpleEditor', 'Search String?',
                                initialvalue=self.target)
        if self.target:
            where = self.text.search(self.target, INSERT, END, nocase=True)
            if where:
                pastit = '{}+{}c'.format(where, len(self.target))
                self.text.tag_add(SEL, where, pastit)
                self.text.mark_set(INSERT, pastit)
                self.text.see(INSERT)
                self.text.focus()


class Progress():
    """ threaded progress bar for tkinter gui """
    def __init__(self, parent, side=LEFT):
        import threading
        self.maximum = 100
        self.interval = 10
        self.progressbar = Progressbar(parent, orient=HORIZONTAL,
                                           mode="indeterminate",
                                           maximum=self.maximum)
        #self.progressbar.grid(row=row, column=column,
        #                      columnspan=columnspan, sticky="we")
        self.progressbar.pack(fill=X,side=side)
        self.thread = threading.Thread()
        self.thread.__init__(target=self.progressbar.start(self.interval),
                             args=())
        self.thread.start()

    def pb_stop(self):
        """ stops the progress bar """
        if not self.thread.isAlive():
            VALUE = self.progressbar["value"]
            self.progressbar.stop()
            self.progressbar["value"] = VALUE

    def pb_start(self):
        """ starts the progress bar """
        if not self.thread.isAlive():
            VALUE = self.progressbar["value"]
            self.progressbar.configure(mode="indeterminate",
                                       maximum=self.maximum,
                                       value=VALUE)
            self.progressbar.start(self.interval)

    def pb_clear(self):
        """ stops the progress bar """
        if not self.thread.isAlive():
            self.progressbar.stop()
            self.progressbar.configure(mode="determinate", value=0)

    def pb_complete(self):
        """ stops the progress bar and fills it """
        if not self.thread.isAlive():
            self.progressbar.stop()
            self.progressbar.configure(mode="determinate",
                                       maximum=self.maximum,
                                       value=self.maximum)

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
import sys,os,traceback,subprocess
import glob,platform,shutil
import pickle
import threading,time
from PySide2 import QtCore
from PySide2.QtWidgets import *
from PySide2.QtGui import *

import pandas as pd
import numpy as np
from Bio import SeqIO
from . import tools, app, widgets, tables

home = os.path.expanduser("~")
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
logoimg = os.path.join(module_path, 'logo.png')

class pygenefinderApp(QMainWindow):
    """GUI Application using PySide2 widgets"""
    def __init__(self, filenames=[], project=None):

        QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("pygenefinder")

        self.setWindowIcon(QIcon(logoimg))
        self.create_menu()
        self.main = QSplitter(self)
        #screen_resolution = QDesktopWidget().screenGeometry()
        screen = QGuiApplication.screens()[0]
        screen_resolution  = screen.geometry()
        if screen_resolution.height() > 1080:
            fac=0.8
            x=150; y=150
        else:
            fac=1
            x=0; y=0
        width, height = screen_resolution.width()*fac, screen_resolution.height()*fac
        self.setGeometry(QtCore.QRect(x, y, width, height))
        self.main.setFocus()
        self.setCentralWidget(self.main)
        self.setup_gui()
        self.clear_project()

        if project != None:
            self.load_project(project)
        self.threadpool = QtCore.QThreadPool()

        return

    def setup_gui(self):
        """Add all GUI elements"""

        self.m = QSplitter(self.main)
        #mainlayout = QHBoxLayout(self.m)
        left = QWidget(self.m)
        #mainlayout.addWidget(left)
        self.opts = AppOptions(parent=self.m)
        dialog = self.opts.showDialog(left, wrap=1, section_wrap=1)
        left.setFixedWidth(250)
        center = QWidget(self.m)#, orientation=QtCore.Qt.Vertical)
        #mainlayout.addWidget(center)
        l = QVBoxLayout(center)
        self.fasta_table = tables.FilesTable(center, app=self, dataframe=pd.DataFrame())
        l.addWidget(self.fasta_table)
        self.fasta_table.setColumnWidth(0,200)
        self.fasta_table.setColumnWidth(1,400)

        self.tabs = QTabWidget(center)
        self.tabs.setTabsClosable(True)
        self.tabs.tabCloseRequested.connect(self.close_tab)
        l.addWidget(self.tabs)

        self.right = right = QWidget(self.m)
        l2 = QVBoxLayout(right)
        #mainlayout.addWidget(right)
        self.right_tabs = QTabWidget(right)
        self.right_tabs.setTabsClosable(True)
        self.right_tabs.tabCloseRequested.connect(self.close_right_tab)
        l2.addWidget(self.right_tabs)
        self.info = QTextEdit(right, readOnly=True)
        #self.info.setStyleSheet("font-family: monospace; font-size: 12px;")
        font = QFont("Monospace")
        font.setPointSize(9)
        font.setStyleHint(QFont.TypeWriter)
        self.info.setFont(font)
        #l2.addWidget(self.info)
        self.right_tabs.addTab(self.info, 'log')
        self.info.setText("Welcome to pygenefinder")
        self.m.setSizes([50,200,150])
        self.m.setStretchFactor(1,0)

        self.statusBar = QStatusBar()
        from . import __version__
        self.projectlabel = QLabel('')
        self.statusBar.addWidget(self.projectlabel, 1)
        self.outdirLabel = QLabel("")
        self.statusBar.addWidget(self.outdirLabel, 1)
        self.progressbar = QProgressBar()
        self.progressbar.setRange(0,1)
        self.statusBar.addWidget(self.progressbar, 2)
        self.setStatusBar(self.statusBar)
        return

    @QtCore.Slot(int)
    def close_tab(self, index):
        """Close current tab"""

        #index = self.tabs.currentIndex()
        name = self.tabs.tabText(index)
        self.tabs.removeTab(index)
        del self.sheets[name]
        return

    @QtCore.Slot(int)
    def close_right_tab(self, index):
        """Close right tab"""

        name = self.right_tabs.tabText(index)
        if name == 'log':
            return
        self.right_tabs.removeTab(index)
        return

    def options_frame(self):
        """Load options dialog"""

        self.opts = AppOptions(parent=self)
        w = self.opts.showDialog(self.main, layout='vertical')
        dialogs.addButton(w, 'Run', self.run)
        return w

    def create_menu(self):
        """Create the menu bar for the application. """

        self.file_menu = QMenu('&File', self)
        #self.file_menu.addAction('&New', self.newProject,
        #        QtCore.Qt.CTRL + QtCore.Qt.Key_N)
        self.file_menu.addAction('&Add Fasta Files', self.load_fasta_files_dialog,
                QtCore.Qt.CTRL + QtCore.Qt.Key_F)
        self.file_menu.addAction('&Load Test Files', self.load_test,
                QtCore.Qt.CTRL + QtCore.Qt.Key_T)
        self.menuBar().addSeparator()
        self.file_menu.addAction('&New Project', self.new_project,
                QtCore.Qt.CTRL + QtCore.Qt.Key_N)
        self.file_menu.addAction('&Open Project', self.load_project_dialog,
                QtCore.Qt.CTRL + QtCore.Qt.Key_O)
        self.file_menu.addAction('&Save Project', self.save_project,
                QtCore.Qt.CTRL + QtCore.Qt.Key_S)
        self.file_menu.addAction('&Save Project As', self.save_project_dialog)
        self.file_menu.addAction('&Quit', self.quit,
                QtCore.Qt.CTRL + QtCore.Qt.Key_Q)
        self.menuBar().addMenu(self.file_menu)

        self.analysis_menu = QMenu('&Analysis', self)
        self.menuBar().addSeparator()
        self.menuBar().addMenu(self.analysis_menu)
        self.analysis_menu.addAction('&Run Blast',
            lambda: self.run_threaded_process(self.run_gene_finder, self.find_genes_completed))
        self.analysis_menu.addAction('&Annotate Selected',
            lambda: self.run_threaded_process(self.annotate_files, self.annotation_completed))
        self.analysis_menu.addAction('&Gene Presence/Absence Matrix', self.show_gene_matrix)

        self.settings_menu = QMenu('&Settings', self)
        self.menuBar().addMenu(self.settings_menu)
        self.settings_menu.addAction('&Set Output Folder', self.set_output_folder)
        self.settings_menu.addAction('&Add Blast Sequences', self.add_sequences_db)
        self.settings_menu.addAction('&Add Trusted Proteins', self.add_trusted_proteins)

        self.help_menu = QMenu('&Help', self)
        self.menuBar().addMenu(self.help_menu)
        self.help_menu.addAction('&Help', self.online_documentation)
        self.help_menu.addAction('&About', self.about)
        return

    def show_info(self, msg):

        self.info.append(msg)
        self.info.verticalScrollBar().setValue(
            self.info.verticalScrollBar().maximum())
        return

    def save_project(self):
        """Save project"""

        if self.proj_file == None:
            self.save_project_dialog()

        filename = self.proj_file
        data={}
        data['inputs'] = self.fasta_table.getDataFrame()
        data['sheets'] = self.sheets
        data['annotations'] = self.annotations
        data['outputdir'] = self.outputdir
        #data['meta'] = self.saveMeta(table)
        self.projectlabel.setText(filename)
        pickle.dump(data, open(filename,'wb'))
        return

    def save_project_dialog(self):
        """Save as project"""

        options = QFileDialog.Options()
        filename, _ = QFileDialog.getSaveFileName(self,"Save Project",
                                                  "","Project files (*.pygf);;All files (*.*)",
                                                  options=options)
        if filename:
            if not os.path.splitext(filename)[1] == '.pygf':
                filename += '.pygf'
            self.proj_file = filename
            self.save_project()
        return

    def new_project(self):
        """New project"""

        self.clear_project(ask=True)
        #self.set_output_folder()
        return

    def clear_project(self, ask=False):
        """Clear all loaded inputs and results"""

        reply=None
        if ask == True:
            reply = QMessageBox.question(self, 'Confirm', "This will clear the current project.\nAre you sure?",
                                        QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.No:
            return False

        self.outputdir = None
        self.sheets = {}
        self.proj_file = None
        self.fasta_table.setDataFrame(pd.DataFrame({'name':[]}))
        self.tabs.clear()
        self.annotations = {}
        self.projectlabel.setText('')
        self.outdirLabel.setText(self.outputdir)
        return

    def load_project(self, filename=None):
        """Load project"""

        self.clear_project()
        data = pickle.load(open(filename,'rb'))
        inputs = data['inputs']
        keys = ['sheets','annotations','outputdir']
        for k in keys:
            if k in data:
                self.__dict__[k] = data[k]

        for s in self.sheets:
            df = self.sheets[s]['data']
            kind = self.sheets[s]['kind']
            self.add_table(s, df, kind)
        ft = self.fasta_table
        ft.setDataFrame(inputs)
        ft.resizeColumns()
        self.proj_file = filename
        self.projectlabel.setText(self.proj_file)
        self.outdirLabel.setText(self.outputdir)
        return

    def load_project_dialog(self):
        """Load project"""

        filename, _ = QFileDialog.getOpenFileName(self, 'Open Project', './',
                                        filter="Project Files(*.pygf);;All Files(*.*)")
        if not filename:
            return
        if not os.path.exists(filename):
            print ('no such file')
        self.load_project(filename)
        return

    def load_test(self):
        """Load test_files"""

        reply = self.clear_project(ask=True)
        if reply == False:
            return
        filenames = glob.glob(os.path.join(app.datadir, '*.fa'))
        self.load_fasta_table(filenames)
        return

    def load_fasta_table(self, filenames):
        """Append/Load fasta inputs into table"""

        if filenames is None or len(filenames) == 0:
            return
        #names = [os.path.splitext(os.path.basename(i))[0] for i in filenames]
        info = [tools.get_fasta_info(f) for f in filenames]
        new = pd.DataFrame(info)
        df = self.fasta_table.model.df
        if len(df)>0:
            new = pd.concat([df,new],sort=False).reset_index(drop=True)
        self.fasta_table.setDataFrame(new)
        self.fasta_table.resizeColumns()
        return

    def load_fasta_files_dialog(self):
        """Load fasta files"""

        options = QFileDialog.Options()
        filenames, _ = QFileDialog.getOpenFileNames(self, 'Open File', './',
                                                    filter="Fasta Files(*.fa *.fna *.fasta);;All Files(*.*)")
        if not filenames:
            return
        self.load_fasta_table(filenames)
        return

    def add_table(self, name, df, kind='default'):
        """Add a table to the results notebook"""

        if kind == 'results':
            t = tables.ResultsTable(self.tabs,  app=self, dataframe=df)
        elif kind == 'features':
            t = tables.FeaturesTable(self.tabs, app=self, dataframe=df)
        else:
            t = tables.DefaultTable(self.tabs, app=self, dataframe=df)
        i = self.tabs.addTab(t, name)
        self.tabs.setCurrentIndex(i)
        self.sheets[name] = {'data': df, 'kind':kind}
        return

    def delete_table(self):
        """Delete current table"""

        name = self.nb.tab(s, 'text')
        del self.sheets[name]
        return

    def get_selected_table(self):
        """current table"""

        #s = self.nb.index(self.nb.select())
        #name = self.nb.tab(s, 'text')
        #t = self.sheets[name]
        return t

    def show_file_info(self, row=None):
        """Show fasta file info"""

        df = self.fasta_table.model.df
        data = df.iloc[row]
        info  = tools.get_fasta_info(data.filename)
        self.info.append(str(info))
        return

    def show_fasta(self):
        """Show selected input fasta file"""

        df = self.fasta_table.model.df
        row = self.fasta_table.getSelectedRow()
        data = df.iloc[row]
        from Bio import SeqIO
        seqs = SeqIO.parse(data.filename, 'fasta')

        for s in seqs:
            ed.text.insert(END, s.format("fasta"))
        return

    def show_gene_matrix(self):
        """Presence/absence matrix of gene features across samples"""

        recs = self.annotations
        x = []
        for name in self.annotations:
            recs = self.annotations[name]
            featsdf = tools.records_to_dataframe(recs)
            featsdf['sample'] = name
            x.append(featsdf)
        x = pd.concat(x)
        x.pivot_table
        m = pd.pivot_table(x, index='sample', columns=['gene','product'], values='start')
        m[m.notnull()] = 1
        m = m.fillna(0)
        m = m.T.reset_index()
        self.add_table('feature matrix', m, kind='default')
        self.feature_matrix = m
        return

    def plot_gene_matrix(self):

        m = self.feature_matrix
        return

    def plot_feature_summary(self):
        """Summary of features"""

        df = self.fasta_table.model.df
        row = self.fasta_table.getSelectedRows()[0]
        data = df.iloc[row]
        name = data.label
        w = widgets.PlotViewer(self)
        import pylab as plt
        fig,ax = plt.subplots(2,1, figsize=(7,5), dpi=65, facecolor=(1,1,1), edgecolor=(0,0,0))
        axs=ax.flat
        w.show_figure(fig)
        recs = self.annotations[name]
        featsdf = tools.records_to_dataframe(recs)
        #featsdf.length.hist(ax=ax)
        featsdf['category'] = featsdf['product'].apply(tools.apply_cat)
        featsdf = featsdf[featsdf.category!='other']
        cats = featsdf.category.value_counts()
        cats.plot.bar(ax=axs[0])
        i = self.right_tabs.addTab(w, name)
        self.right_tabs.setCurrentIndex(i)
        return w

    def show_annotation_report(self):

        return

    def show_feature_table(self):
        """Show features for input file in table"""

        df = self.fasta_table.model.df
        row = self.fasta_table.getSelectedRows()[0]
        data = df.iloc[row]
        name = data.label
        if name not in self.annotations:
            self.info.append('no annotation for this file')
            return
        recs = self.annotations[name]
        if name in self.sheets:
            return
        featsdf = tools.records_to_dataframe(recs)
        featsdf['category'] = featsdf['product'].apply(tools.apply_cat)
        self.add_table(name, featsdf, kind='features')
        return

    def show_genbank_file(self):
        """Show genbank file contents"""

        df = self.fasta_table.model.df
        row = self.fasta_table.getSelectedRows()[0]
        data = df.iloc[row]
        name = data.label
        if name not in self.annotations:
            self.info.append('no annotation for this file')
            return
        w = widgets.FileViewer(self)
        recs = self.annotations[name]
        w.show_records(recs)
        return

    def show_gff_file(self):
        """Show gff file contents"""

        df = self.fasta_table.model.df
        row = self.fasta_table.getSelectedRows()[0]
        data = df.iloc[row]
        name = data.label
        recs = self.annotations[name]
        w = widgets.FileViewer(self)
        w.show_records(recs, format='gff')
        return

    def plot_feature(self, row):
        """Show feature details in viewer"""

        index = self.tabs.currentIndex()
        name = self.tabs.tabText(index)
        table = self.tabs.widget(index)
        df = table.model.df
        data = df.iloc[row]
        recs = self.annotations[name]

        s = widgets.SeqFeaturesViewer(self)
        s.load_records(recs)
        s.set_record(recname=data['id'])
        s.redraw(start = data.start-2000, end=data.end+2000)
        s.show()
        return

    def add_annotatation(self):
        """Add an annotation from an external file"""

        filename, _ = QFileDialog.getOpenFileName(self, 'Open File', './',
                                        filter="Genbank Files(*.gb *.gbk);;All Files(*.*)")
        if not filename:
            return

        df = self.fasta_table.model.df
        row = self.fasta_table.getSelectedRows()[0]
        data = df.iloc[row]
        name = data.label
        index = df.index[row]
        #print (df)
        df.at[index,'genbank'] = filename
        self.fasta_table.refresh()
        recs = list(SeqIO.parse(filename,'gb'))
        self.annotations[name] = recs
        return

    def annotate_files(self, progress_callback):
        """Run gene annotation for input files.
        progress_callback: signal for indicating progress in gui
        """

        retval = self.check_output_folder()
        if retval == 0:
            return
        self.running = True
        self.opts.applyOptions()
        kwds = self.opts.kwds
        overwrite = kwds['overwrite']
        kingdom = kwds['kingdom']
        inputs = self.fasta_table.model.df
        rows = self.fasta_table.getSelectedRows()
        df = inputs.loc[rows]
        #files = inputs.filename
        msg = 'Running genome annotation..\nThis may take some time.'
        progress_callback.emit(msg)
        for i,row in df.iterrows():
            if not os.path.exists(self.outputdir):
                os.makedirs(self.outputdir)
            outfile = os.path.join(self.outputdir, row.label + '.gbk')
            if overwrite == False and row.label in self.annotations:
                msg = 'Annotation exists, skipping'
                progress_callback.emit(msg)
                continue
            progress_callback.emit(row.filename)
            fdf,recs = app.run_annotation(row.filename, threads=int(kwds['threads']), kingdom=kingdom)
            tools.recs_to_genbank(recs, outfile)
            inputs.loc[i,'genbank'] = outfile
            self.fasta_table.refresh()
            self.annotations[row.label] = recs
            featcount = int(sum([len(rec.features) for rec in recs]))
            progress_callback.emit('Found %s genes' %featcount)
            progress_callback.emit('Wrote output to %s' %outfile)
            inputs.loc[i,'genes'] = featcount
        return

    def annotation_completed(self):
        """Gene blasting/finding completed"""

        self.info.append("finished")
        self.progressbar.setRange(0,1)
        df = self.fasta_table.getDataFrame()
        self.fasta_table.refresh()
        self.running = False
        return

    def run_gene_finder(self, progress_callback):
        """Run gene blasting function"""

        self.blast_results = None
        retval = self.check_output_folder()
        if retval == 0:
            return
        inputs = self.fasta_table.getDataFrame()
        if len(inputs) == 0:
            msg = 'you need to load fasta files'
            progress_callback.emit(msg)
            return
        self.opts.applyOptions()
        app.check_databases()
        kwds = self.opts.kwds
        #kwds={'db':'card'}
        db = kwds['db']
        app.fetch_sequence_from_url(db)
        #update inputs from table

        files = inputs.filename

        msg = 'blasting %s files against %s sequences' %(len(files), db)
        progress_callback.emit(msg)
        app.make_target_database(files)
        targfile = os.path.join(app.tempdir, 'targets.fasta')
        bl = app.find_genes(targfile, db, **kwds)
        #bl.to_csv('%s_results.csv' %db)
        msg = 'found %s genes' %len(bl.gene.unique())
        progress_callback.emit(msg)
        m = app.pivot_blast_results(bl)
        self.blast_results = {db: bl}
        return

    def get_tab_names(self):
        return {self.tabs.tabText(index):index for index in range(self.tabs.count())}

    def find_genes_completed(self):
        """Gene blasting/finding completed"""

        self.progressbar.setRange(0,1)
        if self.blast_results == None:
            return
        print("done")
        curr = self.get_tab_names()
        for db in self.blast_results:
            if db in curr:
                self.tabs.removeTab(curr[db])
            self.add_table(db, self.blast_results[db], kind='results')
        self.info.append('done')

        #self.fig.clear()
        #app.plot_heatmap(m.T, fig=self.fig, title='results matrix')
        #self.canvas.draw()
        #m.to_csv('%s_matrix.csv' %db)
        return

    def run_threaded_process(self, process, on_complete):
        """Execute a function in the background with a worker"""

        worker = Worker(fn=process)
        self.threadpool.start(worker)
        worker.signals.finished.connect(on_complete)
        worker.signals.progress.connect(self.progress_fn)
        self.progressbar.setRange(0,0)
        return

    def progress_fn(self, msg):

        print (msg)
        self.info.append(msg)
        self.info.verticalScrollBar().setValue(1)
        return

    def get_selected_gene(self, row):
        """get selected gene"""

        return data

    def show_fasta_sequences(self, row):
        """Show fasta sequences across all samples in blast result"""

        inputs = self.fasta_table.getDataFrame()
        files = inputs.filename
        index = self.tabs.currentIndex()
        name = self.tabs.tabText(index)
        df = self.sheets[name]['data']
        data = df.iloc[row]
        gene = data.gene
        seqs = app.get_gene_hits(df, gene, files, db=name)
        for s in seqs:
            self.show_info(s.format("fasta"))
            #self.info.append(s.format("fasta"))
        return

    def show_protein_sequences(self, row):
        inputs = self.fasta_table.getDataFrame()
        files = inputs.filename
        index = self.tabs.currentIndex()
        name = self.tabs.tabText(index)
        df = self.sheets[name]['data']
        data = df.iloc[row]
        gene = data.gene
        seqs = app.get_gene_hits(df, gene, files, db=name)
        prots = [s.translate() for s in seqs]
        for s in prots:
            self.info.append(s.format("fasta"))
        return

    def show_gene_alignment(self, row):

        inputs = self.fasta_table.getDataFrame()
        files = inputs.filename
        index = self.tabs.currentIndex()
        name = self.tabs.tabText(index)
        df = self.sheets[name]['data']
        data = df.iloc[row]
        gene = data.gene
        seqs = app.get_gene_hits(df, gene, files, db=name)
        #maaft_alignment(seqfile)
        aln = app.get_alignment(seqs)
        #self.info.append(aln.format('clustal'))
        self.show_info(aln.format('clustal'))
        #from Bio import Phylo
        #tree = Phylo.read('temp.dnd', "newick")
        #Phylo.draw_ascii(tree)
        return aln, gene

    def show_phylo_tree(self, row):

        inputs = self.fasta_table.getDataFrame()
        files = inputs.filename
        index = self.tabs.currentIndex()
        name = self.tabs.tabText(index)
        df = self.sheets[name]['data']
        data = df.iloc[row]

        return

    def save_gene_alignment(self, row):

        aln, gene = self.show_gene_alignment(row)
        from Bio import AlignIO
        fname = os.path.join(self.outputdir,'%s.aln' %gene)
        AlignIO.write(aln,fname,'clustal')
        self.info.append('alignment saved to %s' %fname)
        return

    def show_alignment_viewer(self, row):

        inputs = self.fasta_table.getDataFrame()
        files = inputs.filename
        index = self.tabs.currentIndex()
        name = self.tabs.tabText(index)
        df = self.sheets[name]['data']
        data = df.iloc[row]
        gene = data.gene
        seqs = app.get_gene_hits(df, gene, files, db=name)
        aln = app.get_alignment(seqs)
        if not hasattr(self,'alignviewer'):
            self.alignviewer = widgets.AlignmentViewer(self)
        else:
            self.alignviewer.show()
        self.alignviewer.show_alignment(aln, gene)
        return

    def set_output_folder(self):
        """Set the output folder"""

        selected_directory = QFileDialog.getExistingDirectory()
        if selected_directory:
            self.outputdir = selected_directory
        #check it's empty?
        self.outdirLabel.setText(self.outputdir)
        return

    def check_output_folder(self):
        """check if we have an output dir"""

        if self.outputdir == None:
            #QMessageBox.warning(self, 'No output folder set',
            #    'You should set an output folder from the Settings menu')
            self.show_info('You should set an output folder from the Settings menu')
            return 0
        return 1

    def add_sequences_db(self):
        """Add custom search sequences"""

        path = app.customdbdir
        if not os.path.exists(path):
            os.makedirs(path)
        filename, _ = QFileDialog.getOpenFileName(self, 'Open Fasta', './',
                                        filter="Fasta Files(*.fa *.fna *.fasta);;All Files(*.*)")
        if not filename:
            return
        shutil.copy(filename, os.path.join(path, os.path.basename(filename)))
        app.get_files_in_path(app.customdbdir)
        return

    def add_trusted_proteins(self):
        """"""

        path = app.customdbdir
        if not os.path.exists(path):
            os.makedirs(path)
        filename, _ = QFileDialog.getOpenFileName(self, 'Open Protein Fasta', './',
                                    filter="Fasta Files(*.fa *.faa *.fasta);;All Files(*.*)")
        if not filename:
            return
        #app.add_protein_db(filename)
        #options = {'name':{'type':'entry','default':90}}
        #dlg = widgets.DynamicDialog(self, options=options, title='Add Protein File')
        #dlg.show()
        name = os.path.basename(filename)
        #self.trusted = pd.DataFrame({'name':name,'path':filename})
        #df = self.trusted
        #self.opts.widgets['trusted'].addItems([name])
        return

    def _check_snap(self):
        if os.environ.has_key('SNAP_USER_COMMON'):
            print ('running inside snap')
            return True
        return False

    def quit(self):
        self.close()
        return

    def closeEvent(self, event=None):

        if self.proj_file != None and event != None:
            reply = QMessageBox.question(self, 'Confirm', "Save the current project?",
                                            QMessageBox.Cancel | QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if reply == QMessageBox.Cancel:
                event.ignore()
                return
            elif reply == QMessageBox.Yes:
                self.save_project()
        #self.close()
        event.accept()

    def _check_snap(self):
        if os.environ.has_key('SNAP_USER_COMMON'):
            print ('running inside snap')
            return True
        return False

    def online_documentation(self,event=None):
        """Open the online documentation"""

        import webbrowser
        link='https://github.com/dmnfarrell/pygenefinder'
        webbrowser.open(link,autoraise=1)

        #web = QWebEngineView(QUrl('https://github.com/dmnfarrell/pygenefinder'))
        #web.load(QUrl())
        #web.show()
        return

    def about(self):

        from . import __version__
        import matplotlib
        import PySide2
        pandasver = pd.__version__
        pythonver = platform.python_version()
        mplver = matplotlib.__version__
        qtver = PySide2.QtCore.__version__
        if self._check_snap == True:
            snap='(snap)'
        else:
            snap=''

        text='pygenefinder GUI\n'\
                +'version '+__version__+snap+'\n'\
                +'Copyright (C) Damien Farrell 2019-\n'\
                +'This program is free software; you can redistribute it and/or '\
                +'modify it under the terms of the GNU General Public License '\
                +'as published by the Free Software Foundation; either version 3 '\
                +'of the License, or (at your option) any later version.\n'\
                +'Using Python v%s, PySide2 v%s\n' %(pythonver, qtver)\
                +'pandas v%s, matplotlib v%s' %(pandasver,mplver)

        msg = QMessageBox.about(self, "About", text)
        return

#https://www.learnpyqt.com/courses/concurrent-execution/multithreading-pyqt-applications-qthreadpool/
class Worker(QtCore.QRunnable):
    """Worker thread for running background tasks."""

    def __init__(self, fn, *args, **kwargs):
        super(Worker, self).__init__()
        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()
        self.kwargs['progress_callback'] = self.signals.progress

    @QtCore.Slot()
    def run(self):
        try:
            result = self.fn(
                *self.args, **self.kwargs,
            )
        except:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        else:
            self.signals.result.emit(result)
        finally:
            self.signals.finished.emit()

class WorkerSignals(QtCore.QObject):
    """
    Defines the signals available from a running worker thread.
    Supported signals are:
    finished
        No data
    error
        `tuple` (exctype, value, traceback.format_exc() )
    result
        `object` data returned from processing, anything
    """
    finished = QtCore.Signal()
    error = QtCore.Signal(tuple)
    result = QtCore.Signal(object)
    progress = QtCore.Signal(str)

class AppOptions(widgets.BaseOptions):
    """Class to provide a dialog for global plot options"""

    def __init__(self, parent=None):
        """Setup variables"""
        self.parent = parent
        self.kwds = {}
        dbs = app.db_names
        kingdom = ['bacteria','viruses','archaea']
        trusted = app.get_files_in_path(app.trustedproteindir)
        cpus = [str(i) for i in range(1,os.cpu_count()+1)]
        self.groups = {'general':['threads','overwrite'],
                       'blast':['db','identity','coverage','multiple hits'],
                       'annotation':['kingdom','trusted','hmmer']}
        self.opts = {'threads':{'type':'combobox','default':4,'items':cpus},
                    'overwrite':{'type':'checkbox','default':True},
                    'db':{'type':'combobox','default':'card',
                    'items':dbs,'label':'database'},
                    'identity':{'type':'entry','default':90},
                    'coverage':{'type':'entry','default':50},
                    'multiple hits':{'type':'checkbox','default':False},
                    'kingdom':{'type':'combobox','default':'bacteria',
                    'items':kingdom,'label':'kingdom'},
                    'trusted':{'type':'combobox','default':'',
                    'items':trusted,'label':'use trusted'},
                    'hmmer':{'type':'checkbox','default':True},
                    }
        return

def main():
    "Run the application"

    import sys, os
    from argparse import ArgumentParser
    parser = ArgumentParser(description='pygenefinder gui tool')
    parser.add_argument("-f", "--fasta", dest="filenames",default=[],
                        help="input fasta file", metavar="FILE")
    parser.add_argument("-p", "--proj", dest="project",default=None,
                        help="load .pygf project file", metavar="FILE")
    args = vars(parser.parse_args())

    app = QApplication(sys.argv)
    aw = pygenefinderApp(**args)
    aw.show()
    app.exec_()

if __name__ == '__main__':
    main()

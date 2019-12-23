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
import sys,os,traceback,subprocess,glob,platform
import pickle
import threading,time
from PySide2 import QtCore
from PySide2.QtWidgets import *
from PySide2.QtGui import *

import matplotlib
import pandas as pd
import numpy as np
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
        screen_resolution = QDesktopWidget().screenGeometry()
        width, height = screen_resolution.width()*0.7, screen_resolution.height()*.7
        self.setGeometry(QtCore.QRect(200, 200, width, height))
        center = QDesktopWidget().availableGeometry().center()

        self.filenames = filenames
        #self.inputs = []
        self.annotations = {}
        self.outputdir = os.path.join(home,'amr_results')
        self.sheets = {}
        self.annotations = {}
        self.setup_gui()
        #self.new_project()

        self.main.setFocus()
        self.setCentralWidget(self.main)
        self.statusBar = QStatusBar()
        from . import __version__
        self.statusLabel = QLabel("pygenefinder %s" %__version__)
        self.statusBar.addWidget(self.statusLabel, 1)
        self.progressbar = QProgressBar()
        self.progressbar.setRange(0,1)
        self.statusBar.addWidget(self.progressbar, 2)
        #self.progressbar.setValue(10)
        self.setStatusBar(self.statusBar)
        if project != None:
            self.load_project(project)
        self.threadpool = QtCore.QThreadPool()
        return

    def setup_gui(self):
        """Add all GUI elements"""

        self.m = QSplitter(self.main)
        mainlayout = QHBoxLayout(self.m)
        left = QWidget(self.m)
        mainlayout.addWidget(left)
        self.opts = AppOptions(parent=self.m)
        dialog = self.opts.showDialog(left, wrap=2)

        center = QSplitter(self.m, orientation=QtCore.Qt.Vertical)
        mainlayout.addWidget(center)
        l = QVBoxLayout(center)
        self.fasta_table = tables.FilesTable(center, app=self, dataframe=pd.DataFrame())
        self.load_test()
        l.addWidget(self.fasta_table)
        self.fasta_table.setColumnWidth(0,200)
        self.fasta_table.setColumnWidth(1,400)
        self.tabs = QTabWidget(center)
        self.tabs.setTabsClosable(True)
        self.tabs.tabCloseRequested.connect(self.close_tab)
        l.addWidget(self.tabs)
        center.setSizes([50,100])

        right = QWidget(self.m)
        mainlayout.addWidget(right)
        self.info = QTextEdit(right, readOnly=True)
        l = QVBoxLayout(right)
        l.addWidget(self.info)
        self.info.setText("Welcome to pygenefinder")
        self.m.setSizes([45,200,100])
        return

    @QtCore.Slot(int)
    def close_tab(self, index):
        """Close current tab"""

        #index = self.tabs.currentIndex()
        name = self.tabs.tabText(index)
        self.tabs.removeTab(index)
        del self.sheets[name]
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
        self.file_menu.addAction('&Load Fasta Files', self.load_fasta_files,
                QtCore.Qt.CTRL + QtCore.Qt.Key_F)
        self.file_menu.addAction('&Load Test Files', self.load_test,
                QtCore.Qt.CTRL + QtCore.Qt.Key_T)
        self.file_menu.addAction('&Open Project', self.load_project,
                QtCore.Qt.CTRL + QtCore.Qt.Key_O)
        self.file_menu.addAction('&Save Project', self.save_project,
                QtCore.Qt.CTRL + QtCore.Qt.Key_S)
        self.file_menu.addAction('&Quit', self.quit,
                QtCore.Qt.CTRL + QtCore.Qt.Key_Q)
        self.menuBar().addMenu(self.file_menu)

        self.analysis_menu = QMenu('&Analysis', self)
        self.menuBar().addSeparator()
        self.menuBar().addMenu(self.analysis_menu)
        self.analysis_menu.addAction('&Run',
            lambda: self.run_threaded_process(self.run_gene_finder, self.find_genes_completed))
        self.analysis_menu.addAction('&Annotate',
            lambda: self.run_threaded_process(self.annotate_files, self.annotation_completed))

        self.help_menu = QMenu('&Help', self)
        self.menuBar().addSeparator()
        self.menuBar().addMenu(self.help_menu)
        self.help_menu.addAction('&Help', self.online_documentation)
        self.help_menu.addAction('&About', self.about)
        return

    def save_project(self, filename='test.pygf'):
        """Save project"""

        data={}
        data['inputs'] = self.fasta_table.getDataFrame()
        data['sheets'] = self.sheets
        #data['meta'] = self.saveMeta(table)
        pickle.dump(data, open(filename,'wb'))
        return

    def save_as_project(self):
        """Save as project"""

        options = QFileDialog.Options()
        filename, _ = QFileDialog.getSaveFileName(self,"Save Project",
                                                  "","Dxpl Files (*.dexpl);;All files (*.*)",
                                                  options=options)

    def load_project(self, filename=None):
        """Load project"""

        filename='test.pygf'
        self.clear_project()
        data = pickle.load(open(filename,'rb'))
        inputs = data['inputs']
        if 'sheets' in data:
            self.sheets = data['sheets']
        for s in self.sheets:
            self.add_table(s, self.sheets[s])
        ft = self.fasta_table
        ft.setDataFrame(inputs)

        #self.fasta_table.redraw()
        self.filename = filename
        #self.main.title('pygenefinder: %s' %filename)
        return

    def load_project_dialog(self):
        """Load project"""

        filename, _ = QFileDialog.getOpenFileName(self, 'Open File', './',
                                        filter="All Files(*.*);;Fasta Files(*.fa)")
        if not filename:
            return
        if not os.path.exists(filename):
            print ('no such file')
        self.load_project(filename)
        return

    def clear_project(self):
        """Clear all loaded inputs and results"""

        #self.inputs=[]
        self.fasta_table.df = pd.DataFrame()
        self.fasta_table.refresh()
        self.sheets={}
        for i in range(self.tabs.count()):
            self.tabs.removeTab(i)
        #self.fig.clear()
        #self.canvas.draw()
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
        #self.clear_project()
        self.load_fasta_table()
        return

    def load_fasta_table(self):
        """Load fasta inputs into table"""

        if self.filenames is None or len(self.filenames) == 0:
            return
        names = [os.path.splitext(os.path.basename(i))[0] for i in self.filenames]
        df = pd.DataFrame({'label':names,'filename':self.filenames})
        self.fasta_table.setDataFrame(df)
        return

    def load_fasta_files(self, filenames=None):
        """Load fasta files"""

        #self.filenames = filenames
        self.filenames = ['test_files/RF15A.fa','test_files/RF15B.fa']
        self.load_fasta_table()
        self.clear_project()
        return

    def load_fasta_files_dialog(self):
        """Load fasta files"""

        options = QFileDialog.Options()
        #from PySide2.QtWidgets import QFileDialog
        filenames, _ = QFileDialog.getOpenFileNames(self, 'Open File', './',
                                                    filter="All Files(*.*);;Fasta Files(*.fa)")
        if not filenames:
            return
        self.load_fasta_files(filenames)
        return

    def add_table(self, name, df, kind='results'):
        """Add a table to the results notebook"""

        if kind == 'results':
            t = tables.ResultsTable(self.tabs,  app=self, dataframe=df)
        elif kind == 'features':
            t = tables.FeaturesTable(self.tabs, app=self, dataframe=df)
        i = self.tabs.addTab(t, name)
        self.tabs.setCurrentIndex(i)
        self.sheets[name] = t.model.df
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
        """SHow selected input fasta file"""

        df = self.fasta_table.model.df
        row = self.fasta_table.getSelectedRow()
        data = df.iloc[row]
        from Bio import SeqIO
        seqs = SeqIO.parse(data.filename, 'fasta')

        for s in seqs:
            ed.text.insert(END, s.format("fasta"))
        return

    def show_feature_table(self):
        """Show features for input file in table"""

        df = self.fasta_table.model.df
        row = self.fasta_table.getSelectedRows()[0]
        data = df.iloc[row]
        name = data.label
        if name in self.sheets:
            return
        #genbank = data.genbank
        featsdf = tools.genbank_to_dataframe(data.genbank)

        self.add_table(name, featsdf, kind='features')
        return

    def show_genbank_file(self):
        """Show genbank file contents"""

        df = self.fasta_table.model.df
        row = self.fasta_table.getSelectedRows()[0]
        data = df.iloc[row]
        name = data.label
        w = widgets.SeqFeaturesViewer(self)
        w.show_records(data.genbank)
        return

    def annotate_files(self, progress_callback):
        """Run gene annotation for input files.
        progress_callback: signal"""

        self.opts.applyOptions()
        kwds = self.opts.kwds
        #self.inputs = self.fasta_table.getDataFrame()
        inputs = self.fasta_table.model.df
        rows = self.fasta_table.getSelectedRows()
        df = inputs.loc[rows]
        files = inputs.filename
        msg = 'Running genome annotation..\nThis may take some time.'
        progress_callback.emit(msg)
        for i,row in df.iterrows():
            if not os.path.exists(self.outputdir):
                os.makedirs(self.outputdir)
            outfile = os.path.join(self.outputdir, row.label + '.gbk')
            inputs.loc[i,'genbank'] = outfile
            if os.path.exists(outfile):
                msg = 'Genbank file exists, skipping'
                progress_callback.emit(msg)
                continue
            progress_callback.emit(row.filename)
            feats = app.annotate_contigs(row.filename, outfile, threads=kwds['threads'])
            self.annotations[row.label] = feats
            progress_callback.emit('Found %s genes' %len(feats))
            progress_callback.emit('Wrote output to %s' %outfile)
        return

    def annotation_completed(self):
        """Gene blasting/finding completed"""

        self.info.append("finished")
        self.progressbar.setRange(0,1)
        df = self.fasta_table.getDataFrame()
        self.fasta_table.refresh()
        return

    def run_gene_finder(self, progress_callback):
        """Run gene blasting function"""

        #print (self.inputs)
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
        self.results = {db: bl}

        return

    def get_tab_names(self):
        return {self.tabs.tabText(index):index for index in range(self.tabs.count())}

    def find_genes_completed(self):
        """Gene blasting/finding completed"""

        print("done")
        self.progressbar.setRange(0,1)
        curr = self.get_tab_names()
        for db in self.results:
            if db in curr:
                self.tabs.removeTab(curr[db])
            self.add_table(db, self.results[db])
        self.info.append('done')
        #self.show_table('matrix',m.T.reset_index())

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

    def get_selected_gene(self):
        """get selected gene"""

        t = self.get_selected_table()
        row = t.getSelectedRow()
        df = t.model.df
        data = df.iloc[row]
        gene = data.gene
        return

    def show_fasta_sequences(self):

        inputs = self.fasta_table.getDataFrame()
        files = inputs.filename
        t = self.get_selected_table()
        row = t.getSelectedRow()
        df = t.model.df
        data = df.iloc[row]
        gene = data.gene

        seqs = app.get_gene_hits(df, gene, files, db='card')
        print()
        for s in seqs:
            print(s.format("fasta"))
        return

    def show_gene_alignment(self):

        inputs = self.fasta_table.getDataFrame()
        files = inputs.filename
        gene = self.get_selected_gene()
        seqs = app.get_gene_hits(self.bl, gene, files, db='card')
        #maaft_alignment(seqfile)
        print ('alignments for gene: %s' %gene)
        app.get_alignment(seqs)
        return

    def _check_snap(self):
        if os.environ.has_key('SNAP_USER_COMMON'):
            print ('running inside snap')
            return True
        return False

    def quit(self):
        self.close()

    def closeEvent(self, ce):
        self.quit()

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
                +'This program is free software; you can redistribute it and/or\n'\
                +'modify it under the terms of the GNU General Public License\n'\
                +'as published by the Free Software Foundation; either version 3\n'\
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

#!/usr/bin/env python

"""
    dataframe table widget and sub-classes.
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
import pandas as pd
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from PySide2 import QtCore, QtGui
from PySide2.QtCore import QObject
from PySide2.QtWidgets import *
from PySide2.QtGui import *
from . import tools, widgets

class ColumnHeader(QHeaderView):
    def __init__(self):
        super(QHeaderView, self).__init__()
        return

class DataFrameWidget(QWidget):
    """Widget containing a tableview and toolbars"""
    def __init__(self, parent=None, dataframe=None, toolbar=False, *args):

        super(DataFrameWidget, self).__init__()
        l = self.layout = QGridLayout()
        l.setSpacing(2)
        self.setLayout(self.layout)
        self.table = DataFrameTable(self, dataframe)
        l.addWidget(self.table, 1, 1)
        if toolbar==True:
            self.createToolbar()
        self.pf = None
        return

class DataFrameTable(QTableView):
    """
    QTableView with pandas DataFrame as model.
    """
    def __init__(self, parent=None, dataframe=None, fontsize=12, *args):

        QTableView.__init__(self)
        self.clicked.connect(self.showSelection)
        #self.doubleClicked.connect(self.handleDoubleClick)
        self.setSelectionBehavior(QTableView.SelectRows)
        #self.setSelectionBehavior(QTableView.SelectColumns)
        #self.horizontalHeader = ColumnHeader()
        header = self.horizontalHeader()
        #header.setResizeMode(QHeaderView.ResizeToContents)
        vh = self.verticalHeader()
        vh.setVisible(True)
        vh.setDefaultSectionSize(28)
        hh = self.horizontalHeader()
        hh.setVisible(True)
        #hh.setStretchLastSection(True)
        #hh.setSectionResizeMode(QHeaderView.Interactive)
        hh.setSectionsMovable(True)
        hh.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        hh.customContextMenuRequested.connect(self.columnHeaderMenu)
        hh.sectionClicked.connect(self.columnClicked)
        self.setDragEnabled(True)
        #self.setSortingEnabled(True)
        self.viewport().setAcceptDrops(True)
        self.setDropIndicatorShown(True)
        self.resizeColumnsToContents()
        self.setCornerButtonEnabled(True)
        self.setSortingEnabled(True)

        self.font = QFont("Arial", fontsize)
        self.setFont(self.font)
        tm = DataFrameModel(dataframe)
        self.setModel(tm)
        self.model = tm
        #self.resizeColumnsToContents()
        self.setWordWrap(False)
        return

    def createToolbar(self):

        self.toolbar = ToolBar(self)
        self.layout.addWidget(self.toolbar, 1, 2)

    def setDataFrame(self, df):

        tm = DataFrameModel(df)
        self.setModel(tm)
        self.model = tm
        return

    def getDataFrame(self):
        return self.model.df

    def load(self, filename=None):
        return

    def save(self):
        return

    def copy(self):
        return

    def paste(self):
        return

    def zoomIn(self):

        s = self.font.pointSize()
        self.font.setPointSize(s+1)
        self.setFont(self.font)
        vh = self.verticalHeader()
        h = vh.defaultSectionSize()
        vh.setDefaultSectionSize(h+2)
        return

    def zoomOut(self):

        s = self.font.pointSize()
        self.font.setPointSize(s-1)
        self.setFont(self.font)
        vh = self.verticalHeader()
        h = vh.defaultSectionSize()
        vh.setDefaultSectionSize(h-2)
        return

    def importFile(self, filename=None, dialog=False, **kwargs):

        if dialog is True:
            options = QFileDialog.Options()
            filename, _ = QFileDialog.getOpenFileName(self,"Import File",
                                                      "","All Files (*);;Text Files (*.txt);;CSV files (*.csv)",
                                                      options=options)
            df = pd.read_csv(filename, **kwargs)
            self.table.model.df = df
        return

    def info(self):

        buf = io.StringIO()
        self.table.model.df.info(verbose=True,buf=buf,memory_usage=True)
        td = dialogs.TextDialog(self, buf.getvalue(), 'Info')
        return

    def showSelection(self, item):

        cellContent = item.data()
        print(cellContent)  # test
        row = item.row()
        model = item.model()
        columnsTotal= model.columnCount(None)
        return

    def getSelectedRows(self):
        rows=[]
        for idx in self.selectionModel().selectedRows():
            rows.append(idx.row())
        return rows

    def getSelectedColumns(self):
        #indexes = self.selectionModel().selectedIndexes()
        return self.selectionModel().selectedColumns()

    def getSelectedDataFrame(self):

        df = self.model.df
        #print (self.selectionModel().selectedRows())
        rows=[]
        for idx in self.selectionModel().selectedRows():
            rows.append(idx.row())
        cols=[]
        #print (self.selectionModel().selectedColumns())
        return df.iloc[rows]

    def handleDoubleClick(self, item):

        cellContent = item.data()
        if item.column() != 0:
            return
        return

    def columnClicked(self, col):

        hheader = self.horizontalHeader()
        df = self.model.df
        self.model.df = df.sort_values(df.columns[col])
        return

    def storeCurrent(self):
        """Store current version of the table before a major change is made"""

        self.prevdf = self.model.df.copy()
        return

    def deleteCells(self, rows, cols, answer=None):
        """Clear the cell contents"""

        if answer == None:
            answer = QMessageBox.question(self, 'Delete Cells?',
                             'Are you sure?', QMessageBox.Yes, QMessageBox.No)
        if not answer:
            return
        self.storeCurrent()
        print (rows, cols)
        self.model.df.iloc[rows,cols] = np.nan
        return

    def editCell(self, item):
        return

    def setRowColor(self, rowIndex, color):
        for j in range(self.columnCount()):
            self.item(rowIndex, j).setBackground(color)

    def columnHeaderMenu(self, pos):

        hheader = self.horizontalHeader()
        column = hheader.logicalIndexAt(hheader.mapFromGlobal(pos))
        #print (column)
        model = self.model
        menu = QMenu(self)
        setIndexAction = menu.addAction("Set as Index")
        #deleteColumnAction = menu.addAction("Delete Column")
        #sortAction = menu.addAction("Sort By")
        action = menu.exec_(self.mapToGlobal(pos))
        if action == setIndexAction:
            self.setIndex()
        #elif action == sortAction:
        #    self.sort(column)
        return

    def keyPressEvent(self, event):

        rows = self.getSelectedRows()
        cols = self.getSelectedColumns()
        if event.key() == QtCore.Qt.Key_Delete:
            self.deleteCells(rows, cols)

    def contextMenuEvent(self, event):
        """Reimplemented to create context menus for cells and empty space."""

        # Determine the logical indices of the cell where click occured
        hheader, vheader = self.horizontalHeader(), self.verticalHeader()
        position = event.globalPos()
        row = vheader.logicalIndexAt(vheader.mapFromGlobal(position))
        column = hheader.logicalIndexAt(hheader.mapFromGlobal(position))
        if row == -1:
            return
        # Show a context menu for empty space at bottom of table...
        self.menu = QMenu(self)
        self.addActions(event, row)
        return

    def addActions(self, event, row):

        menu = self.menu
        copyAction = menu.addAction("Copy")
        action = menu.exec_(self.mapToGlobal(event.pos()))
        if action == copyAction:
            self.copy()

    def setIndex(self):
        return

    def copy(self):

        self.model.df
        return

    def refresh(self):

        self.model.beginResetModel()
        self.model.dataChanged.emit(0,0)
        self.model.endResetModel()

    def importFile(self):
        dialogs.ImportDialog(self)
        return

    def exportTable(self):
        filename, _ = QFileDialog.getSaveFileName(self, "Export",
                                                  "","csv files (*.csv);;All files (*.*)")
        if filename:
            self.model.df.to_csv(filename)
        return

    def addColumn(self):
        """Add a  column"""

        df = self.model.df
        name, ok = QInputDialog().getText(self, "Enter Column Name",
                                             "Name:", QLineEdit.Normal)
        if ok and name:
            if name in df.columns:
                return
            df[name] = pd.Series()
            self.refresh()
        return

    def deleteRows(self):

        rows = self.getSelectedRows()
        answer = QMessageBox.question(self, 'Delete Rows?',
                             'Are you sure?', QMessageBox.Yes, QMessageBox.No)
        if not answer:
            return
        idx = self.model.df.index[rows]
        self.model.df = self.model.df.drop(idx)
        self.refresh()
        return

class DataFrameModel(QtCore.QAbstractTableModel):
    def __init__(self, dataframe=None, *args):
        super(DataFrameModel, self).__init__()
        if dataframe is None:
            self.df = pd.DataFrame()
        else:
            self.df = dataframe
        self.bg = '#F4F4F3'
        return

    def update(self, df):
        print('Updating Model')
        self.df = df

    def rowCount(self, parent=QtCore.QModelIndex()):
        return len(self.df.index)

    def columnCount(self, parent=QtCore.QModelIndex()):
        return len(self.df.columns.values)

    def data(self, index, role=QtCore.Qt.DisplayRole):

        i = index.row()
        j = index.column()
        if role == QtCore.Qt.DisplayRole:
            value = self.df.iloc[i, j]
            if type(value) != str and np.isnan(value):
                return ''
            else:
                return '{0}'.format(value)
        elif (role == QtCore.Qt.EditRole):
            value = self.df.iloc[i, j]
            if type(value) != str and np.isnan(value):
                return ''
            else:
                return value
        elif role == QtCore.Qt.BackgroundRole:
            return QColor(self.bg)

    def headerData(self, col, orientation, role):

        if orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            return self.df.columns[col]
        if orientation == QtCore.Qt.Vertical and role == QtCore.Qt.DisplayRole:
            return self.df.index[col]
        return None

    def sort(self, Ncol, order):
        """Sort table by given column number """

        self.layoutAboutToBeChanged.emit()
        col = self.df.columns[Ncol]
        self.df = self.df.sort_values(col)
        self.layoutChanged.emit()
        return

    def setData(self, index, value, role=QtCore.Qt.EditRole):
        """Set data upon edits"""

        #print (index)
        i = index.row()
        j = index.column()
        curr = self.df.iloc[i,j]
        print (curr, value)
        self.df.iloc[i,j] = value
        #self.dataChanged.emit()
        return True

    def onDataChanged(self):
        print (self.df)

    def flags(self, index):
            return Qt.ItemIsSelectable|Qt.ItemIsEnabled|Qt.ItemIsEditable

class DefaultTable(DataFrameTable):
    """
    QTableView with pandas DataFrame as model.
    """
    def __init__(self, parent=None, app=None, dataframe=None, *args):
        DataFrameTable.__init__(self, parent, dataframe)
        self.app = app
        self.setWordWrap(False)

    def addActions(self, event, row):

        menu = self.menu
        #showSequencesAction = menu.addAction("Show sequences")
        #action = menu.exec_(self.mapToGlobal(event.pos()))
        #if action == showSequencesAction:
        #    self.app.show_fasta_sequences(row)

        return

class FilesTable(DataFrameTable):
    """
    QTableView for files view.
    """
    def __init__(self, parent=None, app=None, dataframe=None, fontsize=12, *args):
        #super(DataFrameTable, self).__init__()
        DataFrameTable.__init__(self, parent, dataframe, fontsize)
        self.app = app
        self.setWordWrap(False)
        header = self.horizontalHeader()

    def addActions(self, event, row):

        menu = self.menu
        fileinfoAction = menu.addAction("File Summary")
        showFeaturesAction = menu.addAction("Show Feature Table")
        showGenbankAction = menu.addAction("Show Genbank File")
        showGFFAction = menu.addAction("Show GFF File")
        #plotSummaryAction = menu.addAction("Plot Summary")
        addAnnotationAction = menu.addAction("Add Annotation From File")
        removeAction = menu.addAction("Remove Selected")
        addColumnAction = menu.addAction("Add Column")
        exportAction = menu.addAction("Export Table")

        action = menu.exec_(self.mapToGlobal(event.pos()))
        # Map the logical row index to a real index for the source model
        #model = self.model

        if action == fileinfoAction:
            #print (row)
            self.app.show_file_info(row)
        #elif action == annotateAction:
            #self.app.annotate_files()
            #lambda: self.run_threaded_process(self.annotate_files, self.annotation_completed))
        elif action == addAnnotationAction:
            self.app.add_annotatation()
        elif action == showFeaturesAction:
            self.app.show_feature_table()
        elif action == showGenbankAction:
            self.app.show_genbank_file()
        elif action == showGFFAction:
            self.app.show_gff_file()
        elif action == removeAction:
            self.deleteRows()
        elif action == exportAction:
            self.exportTable()
        elif action == addColumnAction:
            self.addColumn()
        return

    def edit(self, index, trigger, event):
        """Override edit to disable editing of first two columns"""

        if index.column() < 5:
            return False
        else:
            QTableView.edit(self, index, trigger, event)
        return True

    def refresh(self):
        DataFrameTable.refresh(self)

    def resizeColumns(self):
        self.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)


class ResultsTable(DataFrameTable):
    """
    QTableView with pandas DataFrame as model.
    """
    def __init__(self, parent=None, app=None, dataframe=None, *args):
        DataFrameTable.__init__(self, parent, dataframe)
        self.app = app
        self.setWordWrap(False)

    def addActions(self, event, row):

        menu = self.menu
        showSequencesAction = menu.addAction("Show sequences")
        #showProteinSeqAction = menu.addAction("Protein sequences")
        showAlignmentAction = menu.addAction("Alignment for all samples")
        #showAlignmentAction = menu.addAction("Alignment for all hits")
        #showAlignmentViewerAction = menu.addAction("Show alignment in viewer")
        #showTreeAction = menu.addAction("Phylo tree for gene")
        action = menu.exec_(self.mapToGlobal(event.pos()))
        if action == showSequencesAction:
            self.app.show_fasta_sequences(row)
        elif action == showAlignmentAction:
            self.app.show_gene_alignment(row)
        #elif action == showTreeAction:
        #    self.app.show_phylo_tree(row)
        return

class FeaturesTable(DataFrameTable):
    """
    QTableView with pandas DataFrame as model.
    """
    def __init__(self, parent=None, app=None, dataframe=None,
                 fontsize=12, name=None, *args):

        DataFrameTable.__init__(self, parent, dataframe, fontsize)
        self.app = app
        self.name = name
        self.setWordWrap(False)

    def addActions(self, event, row):

        menu = self.menu
        showfeatureAction = menu.addAction("Draw Feature(s)")
        copysequenceAction = menu.addAction("Copy Sequence")
        showproteinsequencesAction = menu.addAction("Show Protein Sequences (fasta)")
        shownucleotidesequenceAction = menu.addAction("Show Nucleotide Sequence (fasta)")
        findorthologsAction = menu.addAction("Find Orthologs")
        exportAction = menu.addAction("Export Table")
        action = menu.exec_(self.mapToGlobal(event.pos()))
        if action == showfeatureAction:
            self.plot_features()
        elif action == copysequenceAction:
            self.copy_sequence(row)
        elif action == shownucleotidesequenceAction:
            self.show_nucleotide_sequence()
        elif action == showproteinsequencesAction:
            self.show_protein_sequences()
        elif action == findorthologsAction:
            self.find_orthologs(row)
        elif action == exportAction:
            self.exportTable()

    def plot_features(self):

        name = self.name
        rows = self.getSelectedRows()
        data = self.model.df.iloc[rows]
        start = data.start.min()
        end = data.end.max()
        recs = self.app.annotations[name]

        s = widgets.SeqFeaturesViewer(self, title=name)
        s.load_records(recs)
        s.set_record(recname=data.iloc[0]['id'])
        s.slider.setValue(start)
        s.redraw(start, end=end)
        s.show()
        return

    def copy_sequence(self, row):

        clip = QApplication.clipboard()
        data = self.model.df.iloc[row]
        clip.setText(data.translation)
        return

    def show_nucleotide_sequence(self):
        """Show sequence for selected region"""

        name = self.name
        rows = self.getSelectedRows()
        data = self.model.df.iloc[rows]
        start = data.start.min()
        end = data.end.max()
        chrom = data.id.iloc[0]
        print (start, end)
        #reference to annotations
        recs = SeqIO.to_dict(self.app.annotations[name])
        rec = recs[chrom]
        new = rec[start:end]
        new.id = name
        new.description = chrom+':'+str(start)+'-'+str(end)
        self.app.show_info('------------------------------------------------')
        self.app.show_info(new.format('fasta'))
        return

    def show_protein_sequences(self, format='fasta'):
        """Show selected protein sequences as fasta"""

        name = self.name
        rows = self.getSelectedRows()
        data = self.model.df.iloc[rows]
        recs = tools.dataframe_to_seqrecords(data, idkey='locus_tag',
                        seqkey='translation', desckey='product', alphabet='protein')
        self.app.show_info('------------------------------------------------')
        self.app.show_info(name)
        s=''
        for rec in recs:
            s += rec.format(format)
        self.app.show_info(s)
        return

    def find_orthologs(self, row):
        """Find orthologs for selected protein"""

        data = self.model.df.iloc[row]
        name = data.locus_tag+'_'+data.id
        seq = SeqRecord(Seq(data.translation),id=name)
        self.app.find_orthologs(seq, label=name)
        return

    def edit(self, index, trigger, event):
        """Override edit to disable editing"""

        return

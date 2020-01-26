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
import pandas as pd
import numpy as np
from PySide2 import QtCore, QtGui
from PySide2.QtCore import QObject
from PySide2.QtWidgets import *
from PySide2.QtGui import *

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
    def __init__(self, parent=None, dataframe=None, *args):
        #super(DataFrameTable, self).__init__()
        QTableView.__init__(self)
        self.clicked.connect(self.showSelection)
        self.doubleClicked.connect(self.handleDoubleClick)
        self.setSelectionBehavior(QTableView.SelectRows)
        #self.setSelectionBehavior(QTableView.SelectColumns)
        #self.horizontalHeader = ColumnHeader()
        header = self.horizontalHeader()
        #header.setResizeMode(QHeaderView.ResizeToContents)
        vh = self.verticalHeader()
        vh.setVisible(True)
        hh = self.horizontalHeader()
        hh.setVisible(True)
        #hh.setStretchLastSection(True)
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

        font = QFont("Arial", 12)
        self.setFont(font)
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

    def importFile(self, filename=None, dialog=False, **kwargs):

        if dialog is True:
            options = QFileDialog.Options()
            filename, _ = QFileDialog.getOpenFileName(self,"Import File",
                                                      "","All Files (*);;Text Files (*.txt);;CSV files (*.csv)",
                                                      options=options)


            df = pd.read_csv(filename, **kwargs)
            self.table.model.df = df
        return

    def plot(self):
        self.table.pf.replot()
        return

    def createPlotViewer(self, parent):
        """Create a plot widget attached to this table"""

        if self.pf == None:
            self.pf = plotting.PlotViewer(table=self.table, parent=parent)
        return self.pf

    def info(self):

        buf = io.StringIO()
        self.table.model.df.info(verbose=True,buf=buf,memory_usage=True)
        td = dialogs.TextDialog(self, buf.getvalue(), 'Info')
        return

    def showasText(self):
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
        print (item)
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
        print (column)
        model = self.model
        menu = QMenu(self)
        setIndexAction = menu.addAction("Set as Index")
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
        #types = [int, np.int64, float]
        if role == QtCore.Qt.DisplayRole:
            i = index.row()
            j = index.column()
            value = self.df.iloc[i, j]
            if type(value) != str and np.isnan(value):
                return ''
            else:
                return '{0}'.format(value)

        elif role == QtCore.Qt.BackgroundRole:
            return QColor(self.bg)

    def headerData(self, col, orientation, role):
        if orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            return self.df.columns[col]
        if orientation == QtCore.Qt.Vertical and role == QtCore.Qt.DisplayRole:
            return self.df.index[col]
        return None

    def flags(self, index):
        return QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemIsEditable

    def sort(self, Ncol, order):
        """Sort table by given column number """

        self.layoutAboutToBeChanged.emit()
        col = self.df.columns[Ncol]
        self.df = self.df.sort_values(col)
        self.layoutChanged.emit()
        return

    def setData(self):
        self.dataChanged.emit()

class DefaultTable(DataFrameTable):
    """Basic QTableView """
    def __init__(self, parent=None, app=None, dataframe=None, *args):

        DataFrameTable.__init__(self)
        self.app = app
        self.setWordWrap(False)
        header = self.horizontalHeader()

    def refresh(self):
        DataFrameTable.refresh(self)

    def addActions(self, event, row):
        menu = self.menu

class FilesTable(DataFrameTable):
    """
    QTableView for files view.
    """
    def __init__(self, parent=None, app=None, dataframe=None, *args):
        #super(DataFrameTable, self).__init__()
        DataFrameTable.__init__(self)
        self.app = app
        self.setWordWrap(False)
        header = self.horizontalHeader()
        #header.setResizeMode(QHeaderView.ResizeToContents)

    def addActions(self, event, row):

        menu = self.menu
        fileinfoAction = menu.addAction("File Summary")
        addAnnotationAction = menu.addAction("Add Annotation From File")
        showFeaturesAction = menu.addAction("Show Feature Table")
        showGenbankAction = menu.addAction("Show Genbank File")
        showGFFAction = menu.addAction("Show GFF File")
        plotSummaryAction = menu.addAction("Plot Feature Summary")
        removeAction = menu.addAction("Remove Selected")
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
        elif action == plotSummaryAction:
            self.app.plot_feature_summary()
        elif action == removeAction:
            self.deleteRows(row)
        return

    def refresh(self):
        DataFrameTable.refresh(self)

    def resizeColumns(self):
        self.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)

    def deleteRows(self, rows):

        idx = self.model.df.index[rows]
        self.model.df = self.model.df.drop(idx)
        self.refresh()
        return

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
        elif action == showTreeAction:
            self.app.show_phylo_tree(row)
        #elif action == showProteinSeqAction:
        #    self.app.show_protein_sequences(row)
        return

class FeaturesTable(DataFrameTable):
    """
    QTableView with pandas DataFrame as model.
    """
    def __init__(self, parent=None, app=None, dataframe=None, *args):
        #super(DataFrameTable, self).__init__()
        DataFrameTable.__init__(self, parent, dataframe)
        self.app = app
        self.setWordWrap(False)

    def addActions(self, event, row):

        menu = self.menu
        showfeatureAction = menu.addAction("Draw Feature(s)")
        copysequenceAction = menu.addAction("Copy Sequence")
        action = menu.exec_(self.mapToGlobal(event.pos()))
        if action == showfeatureAction:
            self.app.plot_feature(row)
        elif action == copysequenceAction:
            self.app.copy_sequence(row)

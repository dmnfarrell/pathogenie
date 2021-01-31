import sys, os

from PySide2.QtWebEngineWidgets import *
from Bio import AlignIO

home = os.path.expanduser("~")

class App(QMainWindow):
    """GUI Application using PySide2 widgets"""
    def __init__(self, filenames=[], project=None):

        QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("app")

        self.create_menu()
        self.main = QWidget(self)
        self.setGeometry(QtCore.QRect(200, 200, 800, 800))

        mainlayout = QHBoxLayout(self.main)
        self.main.setFocus()

        self.statusBar = QStatusBar()
        center = QSplitter(self.main)
        self.setCentralWidget(self.main)
        mainlayout.addWidget(center)
        #l = QVBoxLayout(center)

        self.tabs = QTabWidget(center)
        self.tabs.setTabsClosable(True)

        webview = QWebEngineView()
        #self.locationEdit = QLineEdit( 'file://land_use.html' )
        #mainlayout.addWidget(self.locationEdit)

        html = '''<html>
            <head>
            <title>A Sample Page</title>
            </head>
            <body>
            <h1>Hello, World!</h1>
            <hr />
            I have nothing to say.
            </body>
            </html>'''

        #webview.setHtml(html)
        #webview.load( QtCore.QUrl.fromLocalFile(path) )
        self.show_tests()
        path = os.path.abspath('test.html')
        webview.load( QtCore.QUrl.fromLocalFile(path) )
        self.tabs.addTab(webview,'feats')

        webview = QWebEngineView()
        path = os.path.abspath('seqs.html')
        webview.load( QtCore.QUrl.fromLocalFile(path) )
        self.tabs.addTab(webview,'seqs')
        self.info = QPlainTextEdit(center)

        return

    def show_tests(self):

        from pybioviz import dashboards, plotters
        from bokeh.plotting import figure, output_file, save

        #p = plotters.test1()
        #gff_file = 'Mbovis_AF212297.gff'
        gff_file = 'temp.gff'
        features = gff_to_features(gff_file)
        p = plotters.plot_features(features)
        output_file("test.html")
        save(p)
        aln = AlignIO.read('temp.aln','clustal')
        p = plotters.plot_sequence_alignment(aln)
        output_file("seqs.html")
        save(p)
        return

    def create_menu(self):
        """Create the menu bar for the application. """

        self.file_menu = QMenu('&File', self)
        #self.file_menu.addAction('&New', self.newProject,
        #        QtCore.Qt.CTRL + QtCore.Qt.Key_N)
        self.file_menu.addAction('&Load Files', self.get_file,
                QtCore.Qt.CTRL + QtCore.Qt.Key_F)
        self.menuBar().addMenu(self.file_menu)

    def get_file(self):
        options = QFileDialog.Options()
        filename, _ = QFileDialog.getSaveFileName(self,"Save Project",
                                                  "","Dxpl Files (*.dexpl);;All files (*.*)",
                                                  options=options)
        return


def gff_to_features(gff_file):
    """Get features from gff file"""

    if gff_file is None or not os.path.exists(gff_file):
        return
    from BCBio import GFF
    in_handle = open(gff_file,'r')
    rec = list(GFF.parse(in_handle))[0]
    in_handle.close()
    return rec.features

app = QApplication(sys.argv)
aw = App()
aw.show()
app.exec_()

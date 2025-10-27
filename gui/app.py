# General imports
import sys
from PyQt5.QtWidgets import QMainWindow, QApplication, QPushButton, QWidget, QAction, QTabWidget, QVBoxLayout
from PyQt5.QtGui import QIcon, QPalette, QColor
from PyQt5.QtCore import pyqtSlot, Qt

# Define the application style sheet
STYLE_SHEET = """
QMainWindow, QWidget {
    background-color: #0a1929;  /* Darker navy blue background */
    color: #ffffff;  /* White text color */
}

QTabWidget::pane {
    border: 1px solid #b71c1c;  /* Deep red border */
    background-color: #0a1929;  /* Match main background */
}

/* Style for all tabs, including nested ones */
QTabBar::tab {
    background-color: #8b0000;  /* Darker red for unselected tabs */
    color: #ffffff;  /* White text color */
    padding: 8px 20px;
    border: 1px solid #b71c1c;  /* Deep red border */
    border-bottom: none;
    border-top-left-radius: 4px;
    border-top-right-radius: 4px;
}

QTabBar::tab:selected {
    background-color: #b71c1c;  /* Brighter red for selected tab */
    border-bottom: none;
}

QTabBar::tab:hover {
    background-color: #a01818;  /* Medium red for hover */
}

/* Ensure nested tabs (tabs within widgets) follow the same style */
QTabWidget QTabWidget::pane {
    border: 1px solid #b71c1c;  /* Deep red border */
    background-color: #0a1929;  /* Match main background */
}

QTabWidget QTabBar::tab {
    background-color: #8b0000;  /* Darker red for unselected tabs */
    border: 1px solid #b71c1c;  /* Deep red border */
}

QPushButton {
    background-color: #1565c0;  /* Lighter blue for buttons */
    color: white;
    border: none;
    padding: 5px 15px;
    border-radius: 3px;
}

QPushButton:hover {
    background-color: #1976d2;  /* Even lighter on hover */
}

QPushButton:pressed {
    background-color: #0d47a1;  /* Darker when pressed */
}

QTableWidget {
    background-color: #0d2137;  /* Slightly lighter than main background */
    color: white;
    gridline-color: #b71c1c;  /* Red grid lines */
}

QTableWidget::item:selected {
    background-color: #1565c0;  /* Lighter blue for selection */
}

QHeaderView::section {
    background-color: #0d2137;  /* Match table background */
    color: white;
    padding: 4px;
    border: 1px solid #b71c1c;  /* Red borders */
}

QSpinBox, QLineEdit {
    background-color: #0d2137;  /* Match table background */
    color: white;
    border: 1px solid #b71c1c;  /* Red borders */
    padding: 2px;
}

/* Style for group boxes (including those in MatrixApp) */
QGroupBox {
    border: 1px solid #b71c1c;  /* Red border to match theme */
    background-color: #0a1929;  /* Dark blue background */
    margin-top: 1.5ex;  /* Space for the title */
    color: white;
}

QGroupBox::title {
    color: white;
    subcontrol-origin: margin;
    subcontrol-position: top center;
    padding: 0 3px;
}

/* Additional container styling */
QScrollArea, QFrame {
    border: 1px solid #b71c1c;  /* Red border */
    background-color: #0a1929;  /* Dark blue background */
}

QLabel {
    color: white;
    border: none;
}
"""

# Import tabs
from gtguiadditions.GT_Calc_Window import GT_Calc_Window
from gtguiadditions.GT_Learning_Window import GT_Learning_Window
from linear_algebra import LinearAlgebraTab
from Glossary.Glossary import *
from Glossary.CollapsibleBox import CollapsibleBox
from LinearAlgebra.matrix_app import MatrixApp
from SetTheory.SetTheory import SetTheoryTabfrom 
from LinearAlgebra.matrixAppLearning import TeachingMatrixApp
class App(QMainWindow):

    def __init__(self):
        super().__init__()
        self.title = 'SageMath GUI'
        self.left = 0
        self.top = 0
        self.width = 1100
        self.height = 900
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        
        self.table_widget = MyTableWidget(self)
        self.setCentralWidget(self.table_widget)
        
        self.show()
    
class MyTableWidget(QWidget):
    
    def __init__(self, parent):
        super(QWidget, self).__init__(parent)
        self.layout = QVBoxLayout(self)
        
        # Initialize tab screen
        self.tabs = QTabWidget()
        self.gtTab = QWidget()
        self.gtLearnTab = QWidget()
        self.laTab = QWidget()
        self.laLearnTab = QWidget()
        
        # Add tabs
        self.tabs.addTab(GT_Calc_Window(self), "Graph Theory")
        self.tabs.addTab(GT_Learning_Window(self), "Graph Theory Learning")
        self.tabs.addTab(MatrixApp(self), "Linear Algebra")
        self.tabs.addTab(MatrixApp(self), "Linear Algebra Learning")
        self.tabs.addTab(SetTheoryTab(self),"Set Theory") # ^^^
        self.tabs.addTab(TeachingMatrixApp(self), "Linear Algebra Learning")

        # Add tabs to widget
        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    
    # Apply the stylesheet to the entire application
    app.setStyleSheet(STYLE_SHEET)
    
    # Set up a dark palette as fallback (for widgets that don't fully respect stylesheets)
    dark_palette = QPalette()
    dark_palette.setColor(QPalette.Window, QColor(43, 43, 43))
    dark_palette.setColor(QPalette.WindowText, Qt.white)
    dark_palette.setColor(QPalette.Base, QColor(50, 50, 50))
    dark_palette.setColor(QPalette.AlternateBase, QColor(53, 53, 53))
    dark_palette.setColor(QPalette.ToolTipBase, Qt.white)
    dark_palette.setColor(QPalette.ToolTipText, Qt.white)
    dark_palette.setColor(QPalette.Text, Qt.white)
    dark_palette.setColor(QPalette.Button, QColor(53, 53, 53))
    dark_palette.setColor(QPalette.ButtonText, Qt.white)
    dark_palette.setColor(QPalette.BrightText, Qt.red)
    dark_palette.setColor(QPalette.Link, QColor(42, 130, 218))
    dark_palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
    dark_palette.setColor(QPalette.HighlightedText, Qt.black)
    app.setPalette(dark_palette)
    
    ex = App()
    sys.exit(app.exec_())
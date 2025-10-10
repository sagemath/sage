# General imports
import sys
from PyQt5.QtWidgets import QMainWindow, QApplication, QPushButton, QWidget, QAction, QTabWidget,QVBoxLayout
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot

# Import custom tabs
from gtguiadditions.GT_Calc_Window import GT_Calc_Window
from linear_algebra import LinearAlgebraTab

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
        self.tabs.addTab(GT_Calc_Window(self), "Graph Theory Learning") # We would use the learner one
        self.tabs.addTab(LinearAlgebraTab(self), "Linear Algebra")
        self.tabs.addTab(LinearAlgebraTab(self), "Linear Algebra Learning") # ^^^
        
        # Add tabs to widget
        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())
import sys
import tempfile
from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QMessageBox, QLabel, QWidget, QVBoxLayout, QHBoxLayout, QLineEdit
from PyQt5.QtGui import QPixmap, QIcon
from PyQt5.QtCore import QSize
from Glossary.Glossary import GlossaryWidget
#from sage.all import set, Set


class SetTheoryTab(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.glossary_window = None  # Store reference
        self.vert_textbox_layout = QHBoxLayout()
        self.vert_textbox = QLineEdit(self)
        self.vert_textbox.setPlaceholderText("ex: {1, 2, 3, 4}")
        self.vert_textbox.setGeometry(240,30,250,40)
        self.vert_textbox_layout.addWidget(self.vert_textbox)

        # New stuff for glossary
        self.glossary_button = QPushButton("Glossary", self)
        self.glossary_button.setGeometry(400, 550, 80, 40)
        self.glossary_button.clicked.connect(self.show_glossary)

        
    def show_glossary(self):
        if self.glossary_window is None:
            self.glossary_window = QMainWindow()
            self.glossary_window.setWindowTitle("Set Theory Glossary")
            self.glossary_window.setCentralWidget(GlossaryWidget("ST"))
        self.glossary_window.show()
        self.glossary_window.raise_()
        self.glossary_window.activateWindow()
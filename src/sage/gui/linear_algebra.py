from PyQt5.QtWidgets import QWidget, QVBoxLayout, QPushButton, QLabel, QMainWindow
from .Glossary.Glossary import GlossaryWidget

class LinearAlgebraTab(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Linear Algebra content"))
        self.show_button = QPushButton("Show linear algebra?")
        layout.addWidget(self.show_button)
        self.show_button.clicked.connect(self.on_show)


                # New stuff for glossary
        self.glossary_button = QPushButton("Glossary")
        layout.addWidget(self.glossary_button)
        self.glossary_button.clicked.connect(self.show_glossary)

        self.glossary_window = None  # Store reference

        self.show_button.clicked.connect(self.on_show)

        # connect a local slot or expose signals for the parent to connect
        self.show_button.clicked.connect(self.on_show)


    def on_show(self):
        print("LinearAlgebraTab: show clicked")

    def show_glossary(self):
        if self.glossary_window is None:
            self.glossary_window = QMainWindow()
            self.glossary_window.setWindowTitle("Linear Algebra Glossary")
            self.glossary_window.setCentralWidget(GlossaryWidget("LA"))
        self.glossary_window.show()
        self.glossary_window.raise_()
        self.glossary_window.activateWindow()

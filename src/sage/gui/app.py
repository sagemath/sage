import sys
from PyQt5.QtWidgets import QMainWindow, QApplication, QWidget, QTabWidget, QVBoxLayout
from PyQt5.QtGui import QPalette, QColor
from PyQt5.QtCore import Qt

# Local package imports (relative)
from .Glossary.CollapsibleBox import CollapsibleBox
from .LinearAlgebra.matrix_app import MatrixApp
from .graph_theory import GraphTheoryTab

STYLE_SHEET = """
QMainWindow, QWidget {
    background-color: #0a1929;
    color: #ffffff;
}
QTabWidget::pane {
    border: 1px solid #b71c1c;
    background-color: #0a1929;
}
QTabBar::tab {
    background-color: #8b0000;
    color: #ffffff;
    padding: 8px 20px;
    border: 1px solid #b71c1c;
    border-bottom: none;
    border-top-left-radius: 4px;
    border-top-right-radius: 4px;
}
QTabBar::tab:selected {
    background-color: #b71c1c;
    border-bottom: none;
}
QPushButton {
    background-color: #1565c0;
    color: white;
    border: none;
    padding: 5px 15px;
    border-radius: 3px;
}
"""


class App(QMainWindow):

    def __init__(self):
        super().__init__()
        self.title = 'SageMath GUI'
        self.left = 0
        self.top = 0
        self._width = 1100
        self._height = 900
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self._width, self._height)
        
        self.table_widget = MyTableWidget(self)
        self.setCentralWidget(self.table_widget)
        
        self.show()
    

class MyTableWidget(QWidget):
    
    def __init__(self, parent):
        super().__init__(parent)
        self._layout = QVBoxLayout(self)

        # Initialize tab screen
        self.tabs = QTabWidget()

        # Graph Theory (packaged tab)
        self.tabs.addTab(GraphTheoryTab(self), "Graph Theory")
        self.tabs.addTab(GraphTheoryTab(self), "Graph Theory Learning")

        # Linear Algebra tabs using the packaged MatrixApp
        self.tabs.addTab(MatrixApp(self), "Linear Algebra")
        self.tabs.addTab(MatrixApp(self), "Linear Algebra Learning")

        # Add tabs to widget
        self._layout.addWidget(self.tabs)
        self.setLayout(self._layout)


def main():
    """Start the GUI application.

    This is the central entrypoint used by the package entry point.
    """
    app = QApplication(sys.argv)
    app.setStyleSheet(STYLE_SHEET)

    # Dark palette fallback
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

    window = App()
    return app.exec_()

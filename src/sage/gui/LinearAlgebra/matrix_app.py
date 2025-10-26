import sys
import numpy as np

# Use package-relative import so Python can find Matrix_gui when
# the package is imported as sage.gui
from .Matrix_gui import Ui_MatrixGui
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QMainWindow
from ..Glossary.Glossary import GlossaryWidget


class MatrixApp(QtWidgets.QWidget, Ui_MatrixGui):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setupUi(self)
        # Connect buttons to their respective functions
        self.rowSpinBox.valueChanged.connect(self.update_tabe_size)
        self.colSpinBox.valueChanged.connect(self.update_tabe_size)
        self.computeButton.clicked.connect(self.compute_diagonalization)
        self.orthogonalityButton.clicked.connect(self.check_orthogonality)
        self.clearButton.clicked.connect(self.clear_all)
        self.glossary_window = None
        self.glossaryButton.clicked.connect(self.show_glossary)
        #  Add SVD button connection
        self.svdButton.clicked.connect(self.compute_svd)
        # Create and add SVD Tab
        self.svdTab = QtWidgets.QWidget()
        self.svdLayout = QtWidgets.QVBoxLayout(self.svdTab)

        self.svdUlabel = QtWidgets.QLabel("Matrix U:")
        self.svdUTable = QtWidgets.QTableWidget()
        self.svdSLabel = QtWidgets.QLabel("Singular Values (S):")
        self.svdSTable = QtWidgets.QTableWidget()
        self.svdVTLabel = QtWidgets.QLabel("Matrix V^T:")
        self.svdVTTable = QtWidgets.QTableWidget()

        self.svdLayout.addWidget(self.svdUlabel)
        self.svdLayout.addWidget(self.svdUTable)
        self.svdLayout.addWidget(self.svdSLabel)
        self.svdLayout.addWidget(self.svdSTable)
        self.svdLayout.addWidget(self.svdVTLabel)
        self.svdLayout.addWidget(self.svdVTTable)

        self.resultsTabs.addTab(self.svdTab, "SVD Result")

    def update_tabe_size(self):
        rows = self.rowSpinBox.value()
        cols = self.colSpinBox.value()
        self.matrixTable.setRowCount(rows)
        self.matrixTable.setColumnCount(cols)

    def read_matrix(self):
        rows = self.matrixTable.rowCount()
        cols = self.matrixTable.columnCount()
        matrix = np.zeros((rows, cols), dtype=float)
        for i in range(rows):
            for j in range(cols):
                item = self.matrixTable.item(i, j)
                if item and item.text():
                    try:
                        matrix[i, j] = float(item.text())
                    except ValueError:
                        matrix[i, j] = 0.0
        return matrix
           
    def display_matrix(self, table_widget, matrix):
        rows, cols = matrix.shape
        table_widget.setRowCount(rows)
        table_widget.setColumnCount(cols)
        for i in range(rows):
            for j in range(cols):
                item = QtWidgets.QTableWidgetItem(str(round(matrix[i, j], 4)))
                table_widget.setItem(i, j, item)

    def compute_diagonalization(self):
        matrix = self.read_matrix()
        try:
            eigvals, eigvecs = np.linalg.eig(matrix)
            D = np.diag(eigvals)
            P = eigvecs
            self.diagStatusLabel.setText("Diagonalization successful.")
            self.display_matrix(self.diagMatrixTable, D)
            self.display_matrix(self.transformMatrixTable, P)
        except np.linalg.LinAlgError:
            self.diagStatusLabel.setText("Diagonalization failed.")

    def check_orthogonality(self):
        matrix = self.read_matrix()
        try:
            product = np.dot(matrix.T, matrix)
            identity = np.eye(matrix.shape[0])
            is_orthogonal = np.allclose(product, identity)
            status = "Matrix is orthogonal." if is_orthogonal else "Matrix is not orthogonal."
            self.orthoStatusLabel.setText(status)
            self.display_matrix(self.orthoProductTable, product)
        except Exception as e:
            self.orthoStatusLabel.setText(f"Error checking orthogonality: {str(e)}")     
    
    def compute_svd(self):
        matrix = self.read_matrix()
        try:
            U, S, VT = np.linalg.svd(matrix)
            self.display_matrix(self.svdUTable, U)
            self.display_matrix(self.svdSTable, np.diag(S))
            self.display_matrix(self.svdVTTable, VT)
        except Exception as e:
            QtWidgets.QMessageBox.warning(self, "SVD Error", str(e))

    def show_glossary(self):
        if self.glossary_window is None:
            self.glossary_window = QMainWindow()
            self.glossary_window.setWindowTitle("Linear Algebra Glossary")
            self.glossary_window.setCentralWidget(GlossaryWidget("LA"))
        self.glossary_window.show()
        self.glossary_window.raise_()
        self.glossary_window.activateWindow()

    def clear_all(self):
        self.matrixTable.clearContents()
        self.diagMatrixTable.clearContents()
        self.transformMatrixTable.clearContents()
        self.orthoProductTable.clearContents()
        self.diagStatusLabel.setText("Diagonalization Status")
        self.orthoStatusLabel.setText("Orthogonality Status")
        self.svdUTable.clearContents()
        self.svdSTable.clearContents()
        self.svdVTTable.clearContents()

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = MatrixApp()
    window.show()
    sys.exit(app.exec())

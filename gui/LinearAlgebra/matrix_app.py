import sys

import numpy as np

# from ui_matrix_gui import Ui_MatrixGui
# Use a package-qualified import so Python can find Matrix_gui when
# `gui/app.py` is executed (working directory is the `gui/` folder).
from LinearAlgebra.Matrix_gui import Ui_MatrixGui
from PyQt5 import QtWidgets


class MatrixApp(QtWidgets.QWidget,Ui_MatrixGui):
    def __init__(self,parent = None):
        super().__init__(parent)
    #   print("MatrixApp class started")   
        self.setupUi(self)
        # Connect buttons to their respective functions
        self.rowSpinBox.valueChanged.connect(self.update_tabe_size)
        self.colSpinBox.valueChanged.connect(self.update_tabe_size)
        self.computeButton.clicked.connect(self.compute_diagonalization)
        self.orthogonalityButton.clicked.connect(self.check_orthogonality)
        self.clearButton.clicked.connect(self.clear_all)
        self.glossaryButton.clicked.connect(self.close)  # Change this to make the glossary button be a glossary button!
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

        self.svdLayout.addWidget(self.svdULabel)
        self.svdLayout.addWidget(self.svdUTable)
        self.svdLayout.addWidget(self.svdSLabel)
        self.svdLayout.addWidget(self.svdSTable)
        self.svdLayout.addWidget(self.svdVTLabel)
        self.svdLayout.addWidget(self.svdVTTable)

        self.resultsTabs.addTab(self.svdTab, "SVD Result")
    def update_tabe_size(self):
        # Get the row and column count for the matrixtable
      # print("update table size")
        rows = self.rowSpinBox.value()
        cols = self.colSpinBox.value()
        self.matrixTable.setRowCount(rows)
        self.matrixTable.setColumnCount(cols)

    def read_matrix(self):
        """Reads matrix input from matrixTable QTableWidget"""
    #   print("Method read_matrix was clicked!")
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
        """Displays a numpy matrix in a QTableWidget"""
    #   print("Method Display_matrix started")
        rows, cols = matrix.shape
        table_widget.setRowCount(rows)
        table_widget.setColumnCount(cols)
        for i in range(rows):
            for j in range(cols):
                item = QtWidgets.QTableWidgetItem(str(round(matrix[i, j], 4)))
                table_widget.setItem(i, j, item)

    def compute_diagonalization(self):
        """Computes diagonalization and displays results"""       
    #   print("Method compute_diagonalization was clicked!")
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
        """Checks if the matrix is orthogonal"""
    #   print("Method check_orthogonality!")
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

    def clear_all(self):
        """Clears all input and output fields"""
     #   print("Clear_all")
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
  # print("if condition-main")
    app = QtWidgets.QApplication(sys.argv)
    window = MatrixApp()
    window.show()
    sys.exit(app.exec())
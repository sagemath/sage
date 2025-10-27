import sys

import numpy as np
from LinearAlgebra.uiTeaching_tab import Ui_MatrixTeachingUI
from PyQt5 import QtCore, QtWidgets, uic
from PyQt5.QtWebEngineWidgets import QWebEngineView


class VideoPopup(QtWidgets.QDialog):
    def __init__(self, video_url, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Instructional Video")
        self.setFixedSize(600, 400)

        layout = QtWidgets.QVBoxLayout(self)
        self.web_view = QWebEngineView()
        self.web_view.setUrl(QtCore.QUrl(video_url))

        close_btn = QtWidgets.QPushButton("Close")
        close_btn.clicked.connect(self.close)

        layout.addWidget(self.web_view)
        layout.addWidget(close_btn)

class TeachingMatrixApp(QtWidgets.QWidget, Ui_MatrixTeachingUI):
    def __init__(self,parent = None):
        super().__init__(parent)
       # uic.loadUi("uiTeaching_tab_only.ui", self)
        self.setupUi(self)
        # Connect spin boxes to update matrix size
        self.rowSpinBox.valueChanged.connect(self.update_matrix_size)
        self.colSpinBox.valueChanged.connect(self.update_matrix_size)
        # Connect buttons
        self.checkAnswerButton.clicked.connect(self.check_answers)
        self.clearButton.clicked.connect(self.clear_all)
        self.exitButton.clicked.connect(self.close)

        # Help buttons
        self.helpDeterminantButton.clicked.connect(lambda: self.show_video_popup("determinant"))
        self.helpEigenvaluesButton.clicked.connect(lambda: self.show_video_popup("eigenvalues"))
        self.helpInverseButton.clicked.connect(lambda: self.show_video_popup("inverse"))

        # Track incorrect attempts
        self.incorrect_attempts = {"determinant": 0, "eigenvalues": 0, "inverse": 0}

        # Video links
        self.video_links = {
            "determinant": "https://www.youtube.com/watch?v=3ROzG6n4yMc",
            "eigenvalues": "https://www.youtube.com/watch?v=COwDHiXAISA",
            "inverse": "https://www.youtube.com/watch?v=kR9rO-6Y2Zk"
        }

    def update_matrix_size(self):
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

    def check_answers(self):
        matrix = self.read_matrix()
      #  print("inside check answer function")
        correct_det = round(np.linalg.det(matrix), 2) if matrix.size > 0 else None
        correct_eigen = [round(val, 2) for val in np.linalg.eigvals(matrix)] if matrix.size > 0 else []
        invertible = "yes" if correct_det and correct_det != 0 else "no"
        print("Correct_det-->",correct_det)
        # Determinant check
        user_det = self.determinantInput.text().strip()
       
        if user_det == str(correct_det):
            QtWidgets.QMessageBox.information(self, "Determinant", "Determinant is Correct!")
        else:
            self.handle_incorrect("determinant")
        print("User_det-->",user_det)
        #sys.exit()
        print("Correct Eigen Value-->",correct_eigen)
        # Eigenvalues check
        user_eigen = self.eigenvaluesInput.text().strip().split(",")
        if sorted([val.strip() for val in user_eigen]) == sorted([str(val) for val in correct_eigen]):
            QtWidgets.QMessageBox.information(self, "Eigenvalues", "Eigenvalues are Correct!")
        else:
            self.handle_incorrect("eigenvalues")
        print ("User_eigen-->", user_eigen)

        # Inverse check
        user_inverse = self.inverseInput.text().strip().lower()
        if user_inverse == invertible:
            QtWidgets.QMessageBox.information(self, "Invertibility", "Inverse is Correct!")
        else:
            self.handle_incorrect("inverse")

    def handle_incorrect(self, topic):
        self.incorrect_attempts[topic] += 1
        if self.incorrect_attempts[topic] == 1:
            QtWidgets.QMessageBox.warning(self, "Hint", f"Review {topic} concept.")
        elif self.incorrect_attempts[topic] == 2:
            QtWidgets.QMessageBox.warning(self, "Hint", f"Use matrix properties for {topic}.")
        else:
            self.show_video_popup(topic)

    def show_video_popup(self, topic):
        video_url = self.video_links.get(topic, "https://www.youtube.com/embed/default_video")
        popup = VideoPopup(video_url, self)
        popup.exec_()

    def clear_all(self):
        self.matrixTable.clearContents()
        self.determinantInput.clear()
        self.eigenvaluesInput.clear()
        self.inverseInput.clear()
        self.incorrect_attempts = {"determinant": 0, "eigenvalues": 0, "inverse": 0}

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = TeachingMatrixApp()
    window.show()
    sys.exit(app.exec_())
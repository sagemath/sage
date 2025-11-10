from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QWidget, QLabel, QLineEdit, QPushButton, QVBoxLayout, QHBoxLayout, QMessageBox
import sys

class SlopeCalculator(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Slope Calculator")
        self.setGeometry(100, 100, 400, 250)
        self.initUI()

    def initUI(self):
        # --- Intro label ---
        self.label_intro = QLabel("Enter two points (x₁, y₁) and (x₂, y₂):")

        # First point inputs
        self.x1_label = QLabel("x₁:")
        self.x1_input = QLineEdit()
        self.y1_label = QLabel("y₁:")
        self.y1_input = QLineEdit()

        # Second point inputs
        self.x2_label = QLabel("x₂:")
        self.x2_input = QLineEdit()
        self.y2_label = QLabel("y₂:")
        self.y2_input = QLineEdit()

        # Compute button
        self.compute_button = QPushButton("Compute Slope")
        self.compute_button.clicked.connect(self.compute_slope)

        # Result label
        self.result_label = QLabel("Slope (m): ---")

        # --- Layout setup ---
        main_layout = QVBoxLayout()
        main_layout.addWidget(self.label_intro)

        # Row 1: x1, y1
        row1_layout = QHBoxLayout()
        row1_layout.addWidget(self.x1_label)
        row1_layout.addWidget(self.x1_input)
        row1_layout.addWidget(self.y1_label)
        row1_layout.addWidget(self.y1_input)
        main_layout.addLayout(row1_layout)

        # Row 2: x2, y2
        row2_layout = QHBoxLayout()
        row2_layout.addWidget(self.x2_label)
        row2_layout.addWidget(self.x2_input)
        row2_layout.addWidget(self.y2_label)
        row2_layout.addWidget(self.y2_input)
        main_layout.addLayout(row2_layout)

        # Add compute button and result
        main_layout.addWidget(self.compute_button)
        main_layout.addWidget(self.result_label)

        self.setLayout(main_layout)

    def compute_slope(self):
        """Compute slope from two points"""
        try:
            x1 = float(self.x1_input.text())
            y1 = float(self.y1_input.text())
            x2 = float(self.x2_input.text())
            y2 = float(self.y2_input.text())

            if x2 == x1:
                QMessageBox.warning(self, "Error", "x₂ and x₁ cannot be the same — slope undefined.")
                return

            slope = (y2 - y1) / (x2 - x1)
            self.result_label.setText(f"Slope (m): {slope:.4f}")

        except ValueError:
            QMessageBox.warning(self, "Input Error", "Please enter valid numbers for all fields.")

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = SlopeCalculator()
    window.show()
    sys.exit(app.exec_())

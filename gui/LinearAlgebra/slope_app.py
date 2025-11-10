import sys
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QWidget, QLabel, QLineEdit, QPushButton, QVBoxLayout, QHBoxLayout, QMessageBox

class SlopeCalculator(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Slope Calculator")
        self.setGeometry(100, 100, 400, 250)
        self.initUI()

    def initUI(self):
        # --- Input fields for two points ---
        self.label_intro = QLabel("Enter two points (x₁, y₁) and (x₂, y₂):")

        # First point
        self.x1_label = QLabel("x₁:")
        self.x1_input = QLineEdit()
        self.y1_label = QLabel("y₁:")
        self.y1_input = QLineEdit()

        # Second point
        self.x2_label = QLabel("x₂:")
        self.x2_input = QLineEdit()
        self.y2_label = QLabel("y₂:")
        self.y2_input = QLineEdit()

        # Button to compute slope
        self.compute_button = QPushButton("Compute Slope")
        self.compute_button.clicked.connect(self.compute_slope)

        # Label to show result
        self.result_label = QLabel("Slope (m): ---")

        # --- Layout setup ---
        layout = QVBoxLayout()

        layout.addWidget(self.label_intro)

        # Row 1: x1, y1
        row1 = QHBoxLayout()
        row1.addWidget(self.x1_label)
        row1.addWidget(self.x1_input)
        row1.addWidget(self.y1_label)
        row1.addWidget(self.y1_input)
        layout.addLayout(row1)

        # Row 2: x2, y2
        row2 = QHBoxLayout()
        row2.addWidget(self.x2_label)
        row2.addWidget(self.x2_input)
        row2.addWidget(self.y2_label)
        row2.addWidget(self.y2_input)
        layout.addLayout(row2)

        layout.addWidget(self.compute_button)
        layout.addWidget(self.result_label)

        self.setLayout(layout)

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
    sys.exit(app.exec())

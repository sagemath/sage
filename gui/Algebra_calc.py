import sys
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QWidget, QLabel, QLineEdit, QPushButton, QVBoxLayout, QHBoxLayout, QMessageBox
import math

class Algebra_calc(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Algebra Calculator")
        self.setGeometry(100, 100, 400, 250)

        main_layout = QVBoxLayout()
    
        # Intro label
        self.label_intro = QLabel("Enter coefficients of quadratic equation ax² + bx + c = 0:")

        # Inputs for a, b, c
        self.a_label = QLabel("a:")
        self.a_input = QLineEdit()
        self.b_label = QLabel("b:")
        self.b_input = QLineEdit()
        self.c_label = QLabel("c:")
        self.c_input = QLineEdit()

        # Button to compute roots
        self.compute_button = QPushButton("Compute Roots")
        self.compute_button.clicked.connect(self.compute_roots)

        # Result label
        self.result_label = QLabel("Roots: ---")

        # Layout setup
        main_layout = QVBoxLayout()
        main_layout.addWidget(self.label_intro)

        # Row for a, b, c inputs
        row_layout = QHBoxLayout()
        row_layout.addWidget(self.a_label)
        row_layout.addWidget(self.a_input)
        row_layout.addWidget(self.b_label)
        row_layout.addWidget(self.b_input)
        row_layout.addWidget(self.c_label)
        row_layout.addWidget(self.c_input)
        main_layout.addLayout(row_layout)

        main_layout.addWidget(self.compute_button)
        main_layout.addWidget(self.result_label)

        self.setLayout(main_layout)
 # --- Input fields for two points ---
        self.slope_intro = QLabel("Enter two points (x₁, y₁) and (x₂, y₂):")

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
        self.slope_result_label = QLabel("Slope (m): ---")

        # --- Layout setup ---
        layout = QVBoxLayout()

        main_layout.addWidget(self.slope_intro)

        # Row 1: x1, y1
        row1 = QHBoxLayout()
        row1.addWidget(self.x1_label)
        row1.addWidget(self.x1_input)
        row1.addWidget(self.y1_label)
        row1.addWidget(self.y1_input)
        main_layout.addLayout(row1)

        # Row 2: x2, y2
        row2 = QHBoxLayout()
        row2.addWidget(self.x2_label)
        row2.addWidget(self.x2_input)
        row2.addWidget(self.y2_label)
        row2.addWidget(self.y2_input)
        main_layout.addLayout(row2)

        main_layout.addWidget(self.compute_button)
        main_layout.addWidget(self.slope_result_label)
        
        self.setLayout(layout)
    def compute_roots(self):
        try:
            a = float(self.a_input.text())
            b = float(self.b_input.text())
            c = float(self.c_input.text())

            if a == 0:
                QMessageBox.warning(self, "Error", "Coefficient 'a' cannot be zero for a quadratic equation.")
                return

            discriminant = b**2 - 4*a*c

            if discriminant < 0:
                real_part = -b / (2 * a)
                imag_part = math.sqrt(-discriminant) / (2 * a)
                root1 = f"{real_part:.4f} + {imag_part:.4f}i"
                root2 = f"{real_part:.4f} - {imag_part:.4f}i"
                self.result_label.setText(f"Roots are complex:\n{root1}, {root2}")
            else:
                root1 = (-b + math.sqrt(discriminant)) / (2*a)
                root2 = (-b - math.sqrt(discriminant)) / (2*a)
                self.result_label.setText(f"Roots:\n{root1:.4f}, {root2:.4f}")

        except ValueError:
            QMessageBox.warning(self, "Input Error", "Please enter valid numbers for all coefficients.")
   
       

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
            self.slope_result_label.setText(f"Slope (m): {slope:.4f}")

        except ValueError:
            QMessageBox.warning(self, "Input Error", "Please enter valid numbers for all fields.")
if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = Algebra_calc()
    window.show()
    sys.exit(app.exec_())

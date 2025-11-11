from PyQt5.QtWidgets import (
    QWidget, QMainWindow, QPushButton, QMessageBox, QLabel,
    QVBoxLayout, QHBoxLayout, QLineEdit, QSizePolicy, QComboBox,
    QGroupBox, QTextEdit
)
from Glossary.Glossary import GlossaryWidget
#from sage.all import Set

class SetTheoryLearningTab(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        # Main layout
        main_layout = QVBoxLayout(self)

        # Top box: display area for a string (read-only, can show multi-line)
        self.display_box = QGroupBox("Display")
        display_layout = QVBoxLayout(self.display_box)
        self.display_view = QTextEdit()
        self.display_view.setReadOnly(True)
        self.display_view.setPlaceholderText("Displayed text will appear here...")
        self.display_view.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        display_layout.addWidget(self.display_view)
        self.display_box.setLayout(display_layout)

        # Bottom box: input baseline with text input and a submit button
        self.input_box = QGroupBox("Input")
        input_layout = QHBoxLayout(self.input_box)
        self.input_field = QLineEdit()
        self.input_field.setPlaceholderText("Type here and press Enter or Submit")
        self.input_field.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.submit_btn = QPushButton("Submit")
        self.submit_btn.setDefault(True)
        input_layout.addWidget(self.input_field)
        input_layout.addWidget(self.submit_btn)
        self.input_box.setLayout(input_layout)

        # Assemble
        main_layout.addWidget(self.display_box, stretch=3)
        main_layout.addWidget(self.input_box, stretch=0)
        self.setLayout(main_layout)

        # Connections
        self.submit_btn.clicked.connect(self._on_submit)
        self.input_field.returnPressed.connect(self._on_submit)

        # Initial example text
        self.set_display_text("Welcome to Set Theory Learning. Enter a string below to display it here.")

    def _on_submit(self):
        text = self.input_field.text().strip()
        if text:
            # append new line to display (keeps history); change to replace if desired
            current = self.display_view.toPlainText()
            if current:
                new_text = current + "\n" + text
            else:
                new_text = text
            self.display_view.setPlainText(new_text)
            self.input_field.clear()
            # Optionally keep focus for quick repeated input
            self.input_field.setFocus()

    def set_display_text(self, text: str):
        """Replace the display with the provided text."""
        self.display_view.setPlainText(text)
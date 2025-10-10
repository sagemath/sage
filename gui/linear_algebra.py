from PyQt5.QtWidgets import QWidget, QVBoxLayout, QPushButton, QLabel

class LinearAlgebraTab(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Linear Algebra content"))
        self.show_button = QPushButton("Show linear algebra?")
        layout.addWidget(self.show_button)
        self.show_button.clicked.connect(self.on_show)

    def on_show(self):
        print("LinearAlgebraTab: show clicked")
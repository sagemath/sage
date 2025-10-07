from PyQt5.QtWidgets import QWidget, QVBoxLayout, QPushButton, QLabel

class GraphTheoryTab(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Graph Theory content"))
        self.show_button = QPushButton("Show graph theory?")
        layout.addWidget(self.show_button)

        # connect a local slot or expose signals for the parent to connect
        self.show_button.clicked.connect(self.on_show)

    def on_show(self):
        print("GraphTheoryTab: show clicked")
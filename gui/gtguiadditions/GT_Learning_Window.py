import sys
import tempfile
from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QMessageBox, QLabel, QWidget, QVBoxLayout, QHBoxLayout, QLineEdit, QComboBox
from PyQt5.QtGui import QPixmap, QIcon
from PyQt5.QtCore import QSize
from Glossary.Glossary import GlossaryWidget
from sage.all import Graph

class GTImageWindow(QWidget):
    def __init__(self, image_path):
        super().__init__()
        self.setWindowTitle("Learning Graph")
        self.resize(400, 400)

        layout = QVBoxLayout()
        label = QLabel()
        pixmap = QPixmap(image_path)
        label.setPixmap(pixmap)
        label.setScaledContents(True)
        layout.addWidget(label)
        self.setLayout(layout)

class GT_Learning_Window(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Learning Graph")
        self.setGeometry(100, 100, 900, 700)

        main_layout = QVBoxLayout()
        main_layout.setSpacing(15)

        # vertex input box
        vert_layout = QHBoxLayout()
        vert_label_box = QVBoxLayout()
        vert_label_box.addWidget(QLabel("Input Vertice Names:"))
        vert_label_box.addWidget(QLabel("(separated by commas)"))
        vert_layout.addLayout(vert_label_box)

        self.vert_textbox = QLineEdit()
        self.vert_textbox.setPlaceholderText("ex: 1, 2, 3, 4")
        vert_layout.addWidget(self.vert_textbox)
        main_layout.addLayout(vert_layout)

        # edge input box
        edge_layout = QHBoxLayout()
        edge_label_box = QVBoxLayout()
        edge_label_box.addWidget(QLabel("Input Edges:"))
        edge_sub = QLabel("(ordered pairs separated by commas)")
        edge_sub.setWordWrap(True)
        edge_label_box.addWidget(edge_sub)
        edge_layout.addLayout(edge_label_box)

        self.edge_textbox = QLineEdit()
        self.edge_textbox.setPlaceholderText("ex: (1,2), (2,3), (1,4)")
        edge_layout.addWidget(self.edge_textbox)
        main_layout.addLayout(edge_layout)

        #button to display graph
        self.display_graph_button = QPushButton("Display Graph")
        self.display_graph_button.clicked.connect(self.on_display_button)
        main_layout.addWidget(self.display_graph_button)

        # degree input box
        degree_layout = QHBoxLayout()
        degree_layout.addWidget(QLabel("Degrees:"))
        self.degree_textbox = QLineEdit()
        self.degree_textbox.setPlaceholderText("ex: 1:2, 2:1, etc")
        degree_layout.addWidget(self.degree_textbox)
        main_layout.addLayout(degree_layout)

        # density input box
        density_layout = QHBoxLayout()
        density_label_box = QVBoxLayout()
        density_label_box.addWidget(QLabel("Density (percentage to two decimal points):"))
        density_layout.addLayout(density_label_box)

        self.density_textbox = QLineEdit()
        self.density_textbox.setPlaceholderText("ex: 50.00")
        density_layout.addWidget(self.density_textbox)
        density_layout.addWidget(QLabel("%"))
        main_layout.addLayout(density_layout)

        # plannar yes/no
        planar_layout = QHBoxLayout()
        planar_layout.addWidget(QLabel("Planar:"))
        self.planar_select = QComboBox()
        self.planar_select.addItems(["Yes", "No"])
        planar_layout.addWidget(self.planar_select)
        main_layout.addLayout(planar_layout)

        # eulerian yes/no
        euler_layout = QHBoxLayout()
        euler_layout.addWidget(QLabel("Eulerian:"))
        self.euler_select = QComboBox()
        self.euler_select.addItems(["Yes", "No"])
        euler_layout.addWidget(self.euler_select)
        main_layout.addLayout(euler_layout)

        # hamiltonian yes/no
        hamilton_layout = QHBoxLayout()
        hamilton_layout.addWidget(QLabel("Hamiltonian:"))
        self.hamilton_select = QComboBox()
        self.hamilton_select.addItems(["Yes", "No"])
        hamilton_layout.addWidget(self.hamilton_select)
        main_layout.addLayout(hamilton_layout)

        # main layout
        self.setLayout(main_layout)

    #create graph from vertices and edges to use in other functions
    def get_graph(self):
        vert_text = self.vert_textbox.text()
        vertices = [v.strip() for v in vert_text.split(',') if v.strip()] #takes input and puts it in correct form for graph
        edge_text = self.edge_textbox.text()
        edge_pairs = []
        for e in edge_text.split('),'): #takes input and puts it in correct form for graph
            e = e.replace('(', '').replace(')', '').strip()
            if e:
                parts = e.split(',')
                if len(parts) == 2:
                    edge_pairs.append((parts[0].strip(), parts[1].strip()))
        G = Graph()
        G.add_vertices(vertices)
        G.add_edges(edge_pairs)
        return G, vertices

    # displays graph
    def on_display_button(self):
        G, _ = self.get_graph()
        tmp_file = tempfile.NamedTemporaryFile(suffix=".png", delete=False) 
        G.plot().save(tmp_file.name) #saves png to temp location
        self.image_window = GTImageWindow(tmp_file.name)
        self.image_window.show()
import sys
import tempfile
from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QMessageBox, QLabel, QWidget, QVBoxLayout, QHBoxLayout, QLineEdit
from PyQt5.QtGui import QPixmap, QIcon
from PyQt5.QtCore import QSize
from Glossary.Glossary import GlossaryWidget
from sage.all import Graph

class GTImageWindow(QWidget):
    def __init__(self, image_path):
        super().__init__()
        self.setWindowTitle("Graph")
        self.resize(400, 400)

        layout = QVBoxLayout()
        label = QLabel()
        pixmap = QPixmap(image_path)
        label.setPixmap(pixmap)
        label.setScaledContents(True)
        layout.addWidget(label)
        self.setLayout(layout)


class GT_Calc_Window(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Graph Calculations")
        self.resize(900, 700)

        self.glossary_window = None

        # set main layout
        main_layout = QVBoxLayout()
        main_layout.setSpacing(15)

        # vertice input
        vert_layout = QHBoxLayout()
        vert_labels = QVBoxLayout()
        vert_labels.addWidget(QLabel("Input Vertice Names:"))
        vert_labels.addWidget(QLabel("(separated by commas)"))
        vert_layout.addLayout(vert_labels)

        self.vert_textbox = QLineEdit()
        self.vert_textbox.setPlaceholderText("ex: 1, 2, 3, 4")
        vert_layout.addWidget(self.vert_textbox)
        main_layout.addLayout(vert_layout)

        # edge input
        edge_layout = QHBoxLayout()
        edge_labels = QVBoxLayout()
        edge_labels.addWidget(QLabel("Input Edges:"))
        edge_sub = QLabel("(ordered pairs separated by commas)")
        edge_sub.setWordWrap(True)
        edge_labels.addWidget(edge_sub)
        edge_layout.addLayout(edge_labels)

        self.edge_textbox = QLineEdit()
        self.edge_textbox.setPlaceholderText("ex: (1,2), (2,3), (1,4)")
        edge_layout.addWidget(self.edge_textbox)
        main_layout.addLayout(edge_layout)

        # buttons
        self.display_graph_button = QPushButton("Display Graph")
        self.display_graph_button.clicked.connect(self.on_display_button)
        main_layout.addWidget(self.display_graph_button)

        button_row_1 = QHBoxLayout() #sets up row for buttons

        self.degree_button = QPushButton("Degrees")
        self.degree_button.clicked.connect(self.on_degree_button)
        self.degreelabel = QLabel()
        self.degreelabel.hide()
        button_row_1.addWidget(self.degree_button)
        button_row_1.addWidget(self.degreelabel)

        self.density_button = QPushButton("Density")
        self.density_button.clicked.connect(self.on_density_button)
        self.densitylabel = QLabel()
        self.densitylabel.hide()
        button_row_1.addWidget(self.density_button)
        button_row_1.addWidget(self.densitylabel)

        main_layout.addLayout(button_row_1)

    
        button_row_2 = QHBoxLayout() #new row for buttons

        self.planar_button = QPushButton("Planar")
        self.planar_button.clicked.connect(self.on_planar_button)
        self.planarlabel = QLabel()
        self.planarlabel.hide()
        button_row_2.addWidget(self.planar_button)
        button_row_2.addWidget(self.planarlabel)

        self.eulerian_button = QPushButton("Eulerian")
        self.eulerian_button.clicked.connect(self.on_eulerian_button)
        self.eulerianlabel = QLabel()
        self.eulerianlabel.hide()
        button_row_2.addWidget(self.eulerian_button)
        button_row_2.addWidget(self.eulerianlabel)

        self.hamiltonian_button = QPushButton("Hamiltonian")
        self.hamiltonian_button.clicked.connect(self.on_hamiltonian_button)
        self.hamiltonianlabel = QLabel()
        self.hamiltonianlabel.hide()
        button_row_2.addWidget(self.hamiltonian_button)
        button_row_2.addWidget(self.hamiltonianlabel)

        main_layout.addLayout(button_row_2)

       
        self.glossary_button = QPushButton("Glossary")
        self.glossary_button.clicked.connect(self.show_glossary)
        main_layout.addWidget(self.glossary_button)

    
        self.setLayout(main_layout)

    # creates the graph
    def get_graph(self):
        vert_text = self.vert_textbox.text()
        vertices = [v.strip() for v in vert_text.split(',') if v.strip()]
        edge_text = self.edge_textbox.text()
        edge_pairs = []
        for e in edge_text.split('),'):
            e = e.replace('(', '').replace(')', '').strip()
            if e:
                parts = e.split(',')
                if len(parts) == 2:
                    edge_pairs.append((parts[0].strip(), parts[1].strip()))
        G = Graph()
        G.add_vertices(vertices)
        G.add_edges(edge_pairs)
        return G, vertices

    def on_display_button(self):
        G, _ = self.get_graph()
        tmp_file = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
        G.plot().save(tmp_file.name) #saves png to temp location
        self.image_window = GTImageWindow(tmp_file.name)
        self.image_window.show()

    
    def on_degree_button(self):
        G, vertices = self.get_graph()
        degrees = G.degree()
        text = "   ".join([f"{v}: {d}" for v, d in zip(vertices, degrees)])
        self.degreelabel.setText(text)
        self.degreelabel.setWordWrap(True)
        self.degreelabel.show()

    def on_density_button(self):
        try:
            G, _ = self.get_graph()
            if len(G.vertices()) < 2: #protects against too few vertices
                self.densitylabel.setText("N/A")
            else:
                density = float(G.density()) * 100
                self.densitylabel.setText(f"{density:.2f}%")
            self.densitylabel.show()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error computing density:\n{e}")

    def on_planar_button(self):
        G, _ = self.get_graph()
        is_planar = G.is_planar()
        self.planarlabel.setText("Yes" if is_planar else "No")
        self.planarlabel.setStyleSheet(f"color: {'green' if is_planar else 'red'};")
        self.planarlabel.show()

    def on_eulerian_button(self):
        G, _ = self.get_graph()
        is_eulerian = G.is_eulerian()
        self.eulerianlabel.setText("Yes" if is_eulerian else "No")
        self.eulerianlabel.setStyleSheet(f"color: {'green' if is_eulerian else 'red'};")
        self.eulerianlabel.show()

    def on_hamiltonian_button(self):
        G, _ = self.get_graph()
        is_hamiltonian = G.is_hamiltonian()
        self.hamiltonianlabel.setText("Yes" if is_hamiltonian else "No")
        self.hamiltonianlabel.setStyleSheet(f"color: {'green' if is_hamiltonian else 'red'};")
        self.hamiltonianlabel.show()

    def show_glossary(self):
        if self.glossary_window is None:
            self.glossary_window = QMainWindow()
            self.glossary_window.setWindowTitle("Graph Theory Glossary")
            self.glossary_window.setCentralWidget(GlossaryWidget("GT"))
        self.glossary_window.show()
        self.glossary_window.raise_()
        self.glossary_window.activateWindow()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = GT_Calc_Window()
    window.show()
    sys.exit(app.exec_())

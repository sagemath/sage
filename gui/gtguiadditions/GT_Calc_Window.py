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
        self.setGeometry(400, 200, 200, 200)

        gtimagelayout = QVBoxLayout()
        gtimagelabel = QLabel(self)
        pixmap = QPixmap(image_path)

        gtimagelabel.setPixmap(pixmap)
        gtimagelabel.setScaledContents(True)  # Scales the image to fit the window

        gtimagelayout.addWidget(gtimagelabel)
        self.setLayout(gtimagelayout)
        self.layout = QHBoxLayout()

        

class GT_Calc_Window(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.textlayout = QHBoxLayout()
        self.setWindowTitle("Testing")
        self.setGeometry(100, 100, 900, 700)  # x, y, width, height

        self.vert_label = QLabel("Input Vertice Names:", self) #create label for vertice inputs
        self.vert_label.setGeometry(52,30,125,40)
        self.vert_sublabel = QLabel ("(seperated by commas)", self)
        self.vert_sublabel.setGeometry(45,55,160,20)
        

        self.vert_textbox_layout = QHBoxLayout()
        self.vert_textbox = QLineEdit(self)
        self.vert_textbox.setPlaceholderText("ex: 1, 2, 3, 4")
        self.vert_textbox.setGeometry(240,30,250,40)
        self.vert_textbox_layout.addWidget(self.vert_textbox)
        

        self.edge_label = QLabel("Input Edges:", self) #create label for vertice inputs
        self.edge_label.setGeometry(75,120,125,40)
        self.edge_sublabel = QLabel ("(ordered pairs seperated by commas)", self)
        self.edge_sublabel.setGeometry(35,145,160,40)
        self.edge_sublabel.setWordWrap(True)
        
        

        self.edge_textbox_layout = QHBoxLayout()
        self.edge_textbox = QLineEdit(self)
        self.edge_textbox.setPlaceholderText("ex: (1,2), (2,3), (1,4)")
        self.edge_textbox.setGeometry(240,125,250,40)
        self.edge_textbox_layout.addWidget(self.edge_textbox)

        
        # Create a button show graph
        self.display_graph_button = QPushButton("Display Graph", self)
        self.display_graph_button.setGeometry(575, 70, 150, 50)# x, y, width, height
        self.display_graph_button.clicked.connect(self.on_display_button)
        gtimagelayout = QVBoxLayout()
        gtimagelayout.addWidget(self.display_graph_button)

        # Create a button density
        self.density_button = QPushButton("Density", self)
        self.density_button.setGeometry(500, 200, 150, 50) # x, y, width, height
        self.density_button.clicked.connect(self.on_density_button)
        self.densitylabel = QLabel("",self)
        self.densitylabel.hide()
        self.textlayout.addWidget(self.density_button)
        self.textlayout.addWidget(self.densitylabel)


        # Create a button degree
        self.degree_button = QPushButton("Degrees", self)
        self.degree_button.setGeometry(200, 200, 150, 50) # x, y, width, height
        self.degree_button.clicked.connect(self.on_degree_button)
        self.degreelabel = QLabel("",self)
        self.degreelabel.hide()
        self.textlayout.addWidget(self.degree_button)
        self.textlayout.addWidget(self.degreelabel)

        # Create a button planar
        self.planar_button = QPushButton("Planar", self)
        self.planar_button.setGeometry(52, 400, 150, 50) # x, y, width, height
        self.planar_button.clicked.connect(self.on_planar_button)
        self.planarlabel = QLabel("",self)
        self.planarlabel.hide()
        self.textlayout.addWidget(self.planar_button)
        self.textlayout.addWidget(self.planarlabel)
        
        # Create a button eulerian
        self.eulerian_button = QPushButton("Eulerian", self)
        self.eulerian_button.setGeometry(350, 400, 150, 50) # x, y, width, height
        self.eulerian_button.clicked.connect(self.on_eulerian_button)
        self.eulerianlabel = QLabel("",self)
        self.eulerianlabel.hide()
        self.textlayout.addWidget(self.eulerian_button)
        self.textlayout.addWidget(self.eulerianlabel)
     
        # Create a button hamiltonian 
        self.hamiltonian_button = QPushButton("Hamiltonian", self)
        self.hamiltonian_button.setGeometry(650, 400, 150, 50) # x, y, width, height
        self.hamiltonian_button.clicked.connect(self.on_hamiltonian_button)
        self.hamiltonianlabel = QLabel("",self)
        self.hamiltonianlabel.hide()
        self.textlayout.addWidget(self.hamiltonian_button)
        self.textlayout.addWidget(self.hamiltonianlabel)

        # New stuff for glossary
        self.glossary_button = QPushButton("Glossary")
        self.layout().addWidget(self.glossary_button)  # or add to your specific layout
        self.glossary_button.clicked.connect(self.show_glossary)

    

    def on_display_button(self): #creates pop up window of graph picture
         # Get vertices
        vert_text = self.vert_textbox.text()
        vertices = [v.strip() for v in vert_text.split(',') if v.strip()]

        # Get edges
        edge_text = self.edge_textbox.text()
        edge_pairs = []
        for e in edge_text.split('),'):
            e = e.replace('(', '').replace(')', '').strip()
            if e:
                parts = e.split(',')
                if len(parts) == 2:
                    edge_pairs.append((parts[0].strip(), parts[1].strip()))

        # Build Sage graph
        G = Graph()
        G.add_vertices(vertices)
        G.add_edges(edge_pairs)

        # Save graph to temp PNG
        tmp_file = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
        G.plot().save(tmp_file.name)

        self.image_window = GTImageWindow(tmp_file.name)
        self.image_window.show()

    def on_degree_button(self):
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
    
        degrees = G.degree()
        degree_text = "   ".join([f"{v}: {d}" for v, d in zip(vertices, degrees)])
    
        self.degreelabel.setText(degree_text)
        self.degreelabel.setGeometry(200,275,150,50) 
        self.degreelabel.setWordWrap(True)
        self.degreelabel.show()

   
    def on_planar_button(self):
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
        
        if G.is_planar():
            self.planarlabel.setText("Yes")
            self.planarlabel.setStyleSheet("color: green;")
            self.planarlabel.setGeometry(120, 475, 150, 50)
        else:
            self.planarlabel.setText("No")
            self.planarlabel.setStyleSheet("color: red;")
            self.planarlabel.setGeometry(120, 475, 150, 50)
        self.planarlabel.show()

    def on_eulerian_button(self):
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
        if G.is_eulerian():
            self.eulerianlabel.setText("Yes")
            self.eulerianlabel.setStyleSheet("color: green;")
            self.eulerianlabel.setGeometry(420, 475, 150, 50)
        else:
            self.eulerianlabel.setText("No")
            self.eulerianlabel.setStyleSheet("color: red;")
            self.eulerianlabel.setGeometry(420, 475, 150, 50)
        self.eulerianlabel.show()

    def on_hamiltonian_button(self):
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
        if G.is_hamiltonian():
            self.hamiltonianlabel.setText("Yes")
            self.hamiltonianlabel.setStyleSheet("color: green;")
            self.hamiltonianlabel.setGeometry(720, 475, 150, 50)
        else:
            self.hamiltonianlabel.setText("No")
            self.hamiltonianlabel.setStyleSheet("color: red;")
            self.hamiltonianlabel.setGeometry(720, 475, 150, 50)
        self.hamiltonianlabel.show()
    
    def show_glossary(self):
        if self.glossary_window is None:
            self.glossary_window = QMainWindow()
            self.glossary_window.setWindowTitle("Graph Theory Glossary")
            self.glossary_window.setCentralWidget(GlossaryWidget("GT"))
        self.glossary_window.show()
        self.glossary_window.raise_()
        self.glossary_window.activateWindow()

    def on_density_button(self):
        try:
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

        # Debug prints (visible in console) to help diagnose malformed input
            print("DEBUG: vertices =", vertices)
            print("DEBUG: edge_pairs =", edge_pairs)

            G = Graph()
            G.add_vertices(vertices)
            G.add_edges(edge_pairs)

        # Guard against too few vertices (density undefined / division by zero)
            if len(G.vertices()) < 2:
                self.densitylabel.setText("N/A")
            else:
            # Convert to float to avoid strange Sage types when formatting
                density = float(G.density())
                density_perc = density * 100.0
                self.densitylabel.setText(f"{density_perc:.2f}%")

            self.densitylabel.setGeometry(555, 275, 150, 50)
            self.densitylabel.show()

        except Exception as exc:
        # Show a message and print full traceback to console for debugging
            QMessageBox.critical(self, "Error", f"Error computing density:\n{exc}")
            import traceback
            traceback.print_exc()

     
       


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = GT_Calc_Window()
    window.show()
    sys.exit(app.exec_())

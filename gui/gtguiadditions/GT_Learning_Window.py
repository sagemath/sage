import sys
import tempfile
from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QMessageBox, QLabel, QWidget, QVBoxLayout, QHBoxLayout, QLineEdit, QComboBox, QDialog
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtMultimedia import QMediaPlayer, QMediaContent
from PyQt5.QtMultimediaWidgets import QVideoWidget
from PyQt5.QtGui import QPixmap, QIcon
from PyQt5.QtCore import QSize, QUrl
from Glossary.Glossary import GlossaryWidget
#from sage.all import Graph

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

        self.glossary_window = None

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
        degree_layout.addWidget(QLabel("Degrees (list in order seperated by commas):"))
        self.degree_textbox = QLineEdit()
        self.degree_textbox.setPlaceholderText("ex: 2,4,1,2 etc")
        degree_layout.addWidget(self.degree_textbox)
        main_layout.addLayout(degree_layout)

        self.degree_check_button = QPushButton("Check Degrees")
        self.degree_check_button.clicked.connect(self.on_degree_check)
        degree_layout.addWidget(self.degree_check_button)


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

        self.density_check_button = QPushButton("Check Density")
        self.density_check_button.clicked.connect(self.on_density_check)
        density_layout.addWidget(self.density_check_button)

        # plannar yes/no
        planar_layout = QHBoxLayout()
        planar_layout.addWidget(QLabel("Planar:"))
        self.planar_select = QComboBox()
        self.planar_select.addItems(["Select","Yes", "No"])
        planar_layout.addWidget(self.planar_select)
        main_layout.addLayout(planar_layout)

        self.planar_check_button = QPushButton("Check Planar")
        self.planar_check_button.clicked.connect(self.on_planar_check)
        planar_layout.addWidget(self.planar_check_button)

        # eulerian yes/no
        euler_layout = QHBoxLayout()
        euler_layout.addWidget(QLabel("Eulerian:"))
        self.euler_select = QComboBox()
        self.euler_select.addItems(["Select","Yes", "No"])
        euler_layout.addWidget(self.euler_select)
        main_layout.addLayout(euler_layout)

        self.euler_check_button = QPushButton("Check Eulerian")
        self.euler_check_button.clicked.connect(self.on_euler_check)
        euler_layout.addWidget(self.euler_check_button)

        # hamiltonian yes/no
        hamilton_layout = QHBoxLayout()
        hamilton_layout.addWidget(QLabel("Hamiltonian:"))
        self.hamilton_select = QComboBox()
        self.hamilton_select.addItems(["Select","Yes", "No"])
        hamilton_layout.addWidget(self.hamilton_select)
        main_layout.addLayout(hamilton_layout)

        self.hamilton_check_button = QPushButton("Check Hamiltonian")
        self.hamilton_check_button.clicked.connect(self.on_hamilton_check)
        hamilton_layout.addWidget(self.hamilton_check_button)

        self.glossary_button = QPushButton("Glossary")
        self.glossary_button.clicked.connect(self.show_glossary)
        main_layout.addWidget(self.glossary_button)

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

    def show_youtube_hint(self, youtube_url):
        #shows pop up of youtube video
        dialog = QDialog(self)
        dialog.setWindowTitle("INCORRECT Here's a Video to Help")

        layout = QVBoxLayout(dialog)

        view = QWebEngineView()
        layout.addWidget(view)

        view.load(QUrl(youtube_url))

        dialog.resize(800, 450)
        dialog.exec_()


    def on_degree_check(self):
        if not hasattr(self, "incorrect_attempts"):
            self.incorrect_attempts = 0

        G, _ = self.get_graph()
        actual_degree = G.degree() 
        input_degree_text = self.degree_textbox.text()
        # Messages for first and second incorrect attempts
        messages = [
            "Incorrect. Recall to find the degree of a vertex, count how many edges are incident to it.",
            "Still not quite right. Make sure you aren't double counting any edges."
        ]

        # change input text to integers
        try:
            input_degree = [int(x.strip()) for x in input_degree_text.split(",")]
        except ValueError:
            QMessageBox.warning(self, "Invalid Input", "Please enter degrees as comma-separated integers.")
            return
        # Compare input to actual degrees
        if input_degree == actual_degree:
            QMessageBox.information(self, "Correct", f"Correct! The Graph degrees are {actual_degree}")
            self.incorrect_attempts = 0  # reset on correct
        else:
            if self.incorrect_attempts == 0:
                QMessageBox.information(self, "Incorrect", messages[0])
            elif self.incorrect_attempts == 1:
                QMessageBox.information(self, "Incorrect", messages[1])
            elif self.incorrect_attempts == 2:
            # Third attempt — show the video instead of a message
                self.show_youtube_hint("https://youtu.be/C4s5j2-Hos4?si=rXe2XdMMF09sl0hQ") 
            else:
            # for anymore inocrrect answers
                QMessageBox.information(self, "Incorrect", "Keep trying!")

        self.incorrect_attempts += 1

    
    def on_density_check(self):
        if not hasattr(self, "incorrect_attempts"):
            self.incorrect_attempts = 0

        G, _ = self.get_graph()
        actual_density = float(G.density()) * 100
        input_density = float(self.density_textbox.text())
        #messages for the first and second attempts
        messages = [
            "Incorrect. Recall to find density, you divide the number of edges present by all possible edges (aka edges in the complete graph)",
            "Still not quite right. Make sure you are multiplying by 100 to express density as a percentage."
        ]

        if abs(actual_density - input_density) < 0.001:
            QMessageBox.information(self, "Correct", f"Correct! The Graph density is {actual_density:.2f}%")
            self.incorrect_attempts = 0  # reset on correct
        else:
            if self.incorrect_attempts == 0:
                QMessageBox.information(self, "Incorrect", messages[0])
            elif self.incorrect_attempts == 1:
                QMessageBox.information(self, "Incorrect", messages[1])
            elif self.incorrect_attempts == 2:
            # Third attempt — show the video instead of a message
                self.show_youtube_hint("https://youtu.be/42ZYhknJhwM?si=JaCxvluQrXUKb1yB")
            else:
            # If they keep getting it wrong after the third time
                QMessageBox.information(self, "Incorrect", "Keep trying!")

            self.incorrect_attempts += 1

    def on_planar_check(self):
        G, _ = self.get_graph()
        actual_planar = G.is_planar()
        input_planar = self.planar_select.currentText()
        if (actual_planar and input_planar == "Yes") or (not actual_planar and input_planar == "No"):
            QMessageBox.information(self, "Correct", "Correct! Your planar answer is right.")
        elif input_planar == "Select":
            QMessageBox.warning(self, "Incomplete", "Please select 'Yes' or 'No'.")
        else:
            self.show_youtube_hint("https://youtu.be/LSkB6jR44aE?si=v_6A-kfw9r8PMBpw")

    def on_euler_check(self):
        G, _ = self.get_graph()
        actual_euler = G.is_eulerian()
        input_euler = self.euler_select.currentText()
        if (actual_euler and input_euler == "Yes") or (not actual_euler and input_euler == "No"):
            QMessageBox.information(self, "Correct", "Correct! Your Eulerian answer is right.")
        elif input_euler == "Select":
            QMessageBox.warning(self, "Incomplete", "Please select 'Yes' or 'No'.")
        else:
            self.show_youtube_hint("https://youtu.be/sw79Z34v0dQ?si=m6Nzq3J8BBUb2Byk")

    def on_hamilton_check(self):
        G, _ = self.get_graph()
        actual_hamilton = G.is_hamiltonian()
        input_hamilton = self.hamilton_select.currentText()
        if (actual_hamilton and input_hamilton == "Yes") or (not actual_hamilton and input_hamilton == "No"):
            QMessageBox.information(self, "Correct", "Correct! Your Hamiltonian answer is right.")
        elif input_hamilton == "Select":
            QMessageBox.warning(self, "Incomplete", "Please select 'Yes' or 'No'.")
        else:
            self.show_youtube_hint("https://youtu.be/2UczS2hQLsI?si=WbEqbrLVXjjw4wlG")

    def show_glossary(self):
        if self.glossary_window is None:
            self.glossary_window = QMainWindow()
            self.glossary_window.setWindowTitle("Graph Theory Glossary")
            self.glossary_window.setCentralWidget(GlossaryWidget("GT"))
        self.glossary_window.show()
        self.glossary_window.raise_()
        self.glossary_window.activateWindow()
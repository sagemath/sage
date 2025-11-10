import sys
import tempfile
from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QMessageBox, QLabel, QWidget, QVBoxLayout, QHBoxLayout, QLineEdit, QComboBox, QDialog
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtMultimedia import QMediaPlayer, QMediaContent
from PyQt5.QtMultimediaWidgets import QVideoWidget
from PyQt5.QtGui import QPixmap, QIcon
from PyQt5.QtCore import QSize, QUrl, Qt
from Glossary.Glossary import GlossaryWidget
from sage.all import Graph

class GTImageWindow(QWidget):
    def __init__(self, image_path):
        super().__init__()
        self.setWindowTitle("Learning Graph")
        self.resize(400, 400)

        self.label = QLabel(alignment=Qt.AlignCenter)
        self.layout = QVBoxLayout(self)
        self.layout.addWidget(self.label)

        # Load the pixmap once
        self.original_pixmap = QPixmap(image_path)
        self.update_scaled_pixmap()

    def resizeEvent(self, event):
        # Smoothly rescale the image when window size changes
        self.update_scaled_pixmap()
        super().resizeEvent(event)

    def update_scaled_pixmap(self):
        if not self.original_pixmap.isNull():
        # Find how much space we actually have inside the label
            available_size = self.label.size()
            original_size = self.original_pixmap.size()

        # Only scale down (never up)
            target_width = min(available_size.width(), original_size.width())
            target_height = min(available_size.height(), original_size.height())

            scaled = self.original_pixmap.scaled(
                target_width,
                target_height,
                Qt.KeepAspectRatio,
                Qt.SmoothTransformation
        )
        self.label.setPixmap(scaled)


class GT_Learning_Window(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Learning Graph")
        self.setGeometry(100, 100, 900, 700)

        self.glossary_window = None

        self.main_layout = QVBoxLayout()
        self.main_layout.setSpacing(15)

        # vertex input box
        self.vert_layout = QHBoxLayout()
        self.vert_label_box = QVBoxLayout()
        self.vert_label_box.addWidget(QLabel("Input Vertice Names:"))
        self.vert_label_box.addWidget(QLabel("(separated by commas)"))
        self.vert_layout.addLayout(self.vert_label_box)

        self.vert_textbox = QLineEdit()
        self.vert_textbox.setPlaceholderText("ex: 1, 2, 3, 4")
        self.vert_layout.addWidget(self.vert_textbox)
        self.main_layout.addLayout(self.vert_layout)

        # edge input box
        self.edge_layout = QHBoxLayout()
        self.edge_label_box = QVBoxLayout()
        self.edge_label_box.addWidget(QLabel("Input Edges:"))
        self.edge_sub = QLabel("(ordered pairs separated by commas)")
        self.edge_sub.setWordWrap(True)
        self.edge_label_box.addWidget(self.edge_sub)
        self.edge_layout.addLayout(self.edge_label_box)

        self.edge_textbox = QLineEdit()
        self.edge_textbox.setPlaceholderText("ex: (1,2), (2,3), (1,4)")
        self.edge_layout.addWidget(self.edge_textbox)
        self.main_layout.addLayout(self.edge_layout)

        quiz_layout = QHBoxLayout()
        self.quiz_vert_edge_button = QPushButton ("Quiz Me with Vertices and Edges")
        self.quiz_vert_edge_button.clicked.connect(self.on_quiz_vert_edge)
        quiz_layout.addWidget(self.quiz_vert_edge_button)

        self.quiz_graph_button = QPushButton ("Quiz Me with a Graph")
        self.quiz_graph_button.clicked.connect(self.on_quiz_graph)
        quiz_layout.addWidget(self.quiz_graph_button)
        self.main_layout.addLayout(quiz_layout)


        #button to display graph
        self.display_graph_button = QPushButton("Display Graph")
        self.display_graph_button.clicked.connect(self.on_display_button)
        self.main_layout.addWidget(self.display_graph_button)

        # degree input box
        degree_layout = QHBoxLayout()
        degree_layout.addWidget(QLabel("Degrees (list in order seperated by commas):"))
        self.degree_textbox = QLineEdit()
        self.degree_textbox.setPlaceholderText("ex: 2,4,1,2 etc")
        degree_layout.addWidget(self.degree_textbox)
        self.main_layout.addLayout(degree_layout)

        self.degree_help_button = QPushButton("Degree Example")
        self.degree_help_button.clicked.connect(self.on_degree_help)
        degree_layout.addWidget(self.degree_help_button)

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
        self.main_layout.addLayout(density_layout)

        self.density_help_button = QPushButton("Density Example")
        self.density_help_button.clicked.connect(self.on_density_help)
        density_layout.addWidget(self.density_help_button)

        self.density_check_button = QPushButton("Check Density")
        self.density_check_button.clicked.connect(self.on_density_check)
        density_layout.addWidget(self.density_check_button)

        # plannar yes/no
        planar_layout = QHBoxLayout()
        planar_layout.addWidget(QLabel("Planar:"))
        self.planar_select = QComboBox()
        self.planar_select.addItems(["Select","Yes", "No"])
        planar_layout.addWidget(self.planar_select)
        self.main_layout.addLayout(planar_layout)

        self.planar_help_button = QPushButton("Planar Example")
        self.planar_help_button.clicked.connect(self.on_planar_help)
        planar_layout.addWidget(self.planar_help_button)

        self.planar_check_button = QPushButton("Check Planar")
        self.planar_check_button.clicked.connect(self.on_planar_check)
        planar_layout.addWidget(self.planar_check_button)

        # eulerian yes/no
        euler_layout = QHBoxLayout()
        euler_layout.addWidget(QLabel("Eulerian:"))
        self.euler_select = QComboBox()
        self.euler_select.addItems(["Select","Yes", "No"])
        euler_layout.addWidget(self.euler_select)
        self.main_layout.addLayout(euler_layout)

        self.eulerian_help_button = QPushButton("Eulerian Example")
        self.eulerian_help_button.clicked.connect(self.on_eulerian_help)
        euler_layout.addWidget(self.eulerian_help_button)

        self.euler_check_button = QPushButton("Check Eulerian")
        self.euler_check_button.clicked.connect(self.on_euler_check)
        euler_layout.addWidget(self.euler_check_button)

        # hamiltonian yes/no
        hamilton_layout = QHBoxLayout()
        hamilton_layout.addWidget(QLabel("Hamiltonian:"))
        self.hamilton_select = QComboBox()
        self.hamilton_select.addItems(["Select","Yes", "No"])
        hamilton_layout.addWidget(self.hamilton_select)
        self.main_layout.addLayout(hamilton_layout)

        self.hamilton_help_button = QPushButton("Hamiltonian Example")
        self.hamilton_help_button.clicked.connect(self.on_hamilton_help)
        hamilton_layout.addWidget(self.hamilton_help_button)

        self.hamilton_check_button = QPushButton("Check Hamiltonian")
        self.hamilton_check_button.clicked.connect(self.on_hamilton_check)
        hamilton_layout.addWidget(self.hamilton_check_button)

        self.glossary_button = QPushButton("Glossary")
        self.glossary_button.clicked.connect(self.show_glossary)
        self.main_layout.addWidget(self.glossary_button)

        # main layout
        self.setLayout(self.main_layout)

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
    
    def clear_layout(self, layout):
     while layout.count():
        item = layout.takeAt(0)
        widget = item.widget()
        if widget:
            widget.deleteLater()
        elif item.layout():
            self.clear_layout(item.layout())
            item.layout().deleteLater()

    def on_quiz_vert_edge(self):
        for layout in [self.vert_layout, self.edge_layout]:
            self.clear_layout(layout)
            self.main_layout.removeItem(layout)
            layout.deleteLater()
        self.quiz_vert_edge_layout = QHBoxLayout()
        
    def on_quiz_graph(self):
        for layout in [self.vert_layout, self.edge_layout]:
            self.clear_layout(layout)
            self.main_layout.removeItem(layout)
            layout.deleteLater()
        
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

    def on_degree_help(self):
        QMessageBox.information (self, "Degree Help", 
        "Ex:\n" 
        "Given the edges (1,2), (2,3), (1,4), (2,4)\n\n" 
        "count the number of times each vertex appears in an edge. \n\n" 
        "The degrees are then: 2,3,1,2\n\n"
        "Given a picture of the graph, simply count the number of edges attached to each vertex")

    def on_density_help(self):
        QMessageBox.information (self, "Density Help", 
        "Ex:\n" 
        "Given the edges (1,2), (2,3), (1,4), (2,4)\n\n" 
        "To find the density as a percentage, we use the formula:\n\n"
        "D=(2 x # of edges)/[(# of vertices)(# of vertices -1)]*100\n\n"
        "for our graph, we have: \n\n"
        "D=(2*4)/[(4)(3)]*100\n\n"
        "D=(8)/(12)*100\n\n"
        "D= 66.67%")

    def on_planar_help(self):
        QMessageBox.information (self, "Planar Help", "There are several methods that can be used to check planarity\n"
        "1: For smaller graphs, try drawing the graph with no edges crossing, then the graph IS planar\n\n"
        "2: if # of vertices is atleast 3 and # of edges < 3 x # of vertices - 6, then the graph is NOT planar\n\n"
        "3: if # of vertices is atleast 3 and there are no 3 cycles and # of edges < 2 x # of vertices - 4 then the graph is NOT planar\n\n")

    def on_eulerian_help(self):
        QMessageBox.information (self, "Eulerian",        
        "Ex:\n" 
        "Given the edges (1,2), (2,3), (1,4), (2,4)\n\n" 
        "To be eulerian, the degree of each vertex must be even\n\n"
        "Since the degree of both vertex 2 and 3 are odd, the graph is not eulerian")

    def on_hamilton_help(self):
        QMessageBox.information (self, "Hamiltonian Help",         
        "Ex:\n" 
        "Given the edges (1,2), (2,3), (1,4), (2,4)\n\n"
        "If each vertex does not have degree atleast (# of vertices)/2, then the graph is NOT hamiltonian\n\n"
        "Since vertex 3 has degree of 1, this graph is not hamiltonian")
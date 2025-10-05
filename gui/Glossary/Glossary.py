from PyQt5 import QtCore, QtGui, QtWidgets
from CollapsibleBox import *
from Content import *

# Ideally, we do not have to make one for each subject and can just import the appropriate content for the glossary (will talk to the group about this)
# Will need to look into finding a way to import math typesetting for this

if __name__ == "__main__":
    import sys
    
    # temporary headers, will import this into a list with their content
    #headers = [
    #    "Graph", "Neighbor Set", "Degree", "Complete Graph",
    #    "Handshaking Lemma", "Planar Graph", "Bipartite Graph",
    #    "Subgraph", "Kuratowski's Theorem", "Path", "Cycle"
    #]
    application = QtWidgets.QApplication(sys.argv)

    window = QtWidgets.QMainWindow()
    window.setCentralWidget(QtWidgets.QWidget())
    dock = QtWidgets.QDockWidget("Glossary")
    window.addDockWidget(QtCore.Qt.LeftDockWidgetArea, dock)

    scroll = QtWidgets.QScrollArea()
    dock.setWidget(scroll)
    content = QtWidgets.QWidget()
    scroll.setWidget(content)
    scroll.setWidgetResizable(True)

    vlay = QtWidgets.QVBoxLayout(content)

    searchBar = QtWidgets.QLineEdit()
    searchBar.setPlaceholderText("Search definitions or theorems...")
    searchBar.setClearButtonEnabled(True)
    vlay.addWidget(searchBar)

    listOfContentBoxes = []
    max_label_width = 0


    definitions=outputContent("GT")
    for word in definitions:
        contentBox = CollapsibleBox(word)
        vlay.addWidget(contentBox)
        lay = QtWidgets.QVBoxLayout()

        contentLabel = QtWidgets.QLabel(definitions[word])
        contentLabel.setStyleSheet("font-size: 10pt;")
        contentLabel.setMinimumSize(contentLabel.sizeHint())
        lay.addWidget(contentLabel)

        contentBox.setContentLayout(lay)
        listOfContentBoxes.append(contentBox)

        # track widest label (this scales the width the widget)
        max_label_width = max(max_label_width, contentLabel.sizeHint().width())

    # Apply the widest label width to all boxes
    for box in listOfContentBoxes:
        box.setMinimumWidth(max_label_width + 40)  # + some padding

    vlay.addStretch()

    window.resize(1024, 768)
    window.show()
    sys.exit(application.exec_())
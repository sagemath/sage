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
    app = QtWidgets.QApplication(sys.argv)

    w = QtWidgets.QMainWindow()
    w.setCentralWidget(QtWidgets.QWidget())
    dock = QtWidgets.QDockWidget("Glossary")
    w.addDockWidget(QtCore.Qt.LeftDockWidgetArea, dock)

    scroll = QtWidgets.QScrollArea()
    dock.setWidget(scroll)
    content = QtWidgets.QWidget()
    scroll.setWidget(content)
    scroll.setWidgetResizable(True)

    vlay = QtWidgets.QVBoxLayout(content)

    boxes = []
    max_label_width = 0


    definitions=outputContent("LA")
    for key in definitions:
        box = CollapsibleBox(key)
        vlay.addWidget(box)
        lay = QtWidgets.QVBoxLayout()

        label = QtWidgets.QLabel(definitions[key])
        label.setStyleSheet("font-size: 10pt;")
        label.setMinimumSize(label.sizeHint())
        lay.addWidget(label)

        box.setContentLayout(lay)
        boxes.append(box)

        # track widest label (this scales the width the widget)
        max_label_width = max(max_label_width, label.sizeHint().width())

    # Apply the widest label width to all boxes
    for box in boxes:
        box.setMinimumWidth(max_label_width + 40)  # + some padding

    vlay.addStretch()

    w.resize(1024, 768)
    w.show()
    sys.exit(app.exec_())
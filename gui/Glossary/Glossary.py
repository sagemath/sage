from PyQt5 import QtCore, QtGui, QtWidgets
from .CollapsibleBox import *
from .Content import *
import sys

# Ideally, we do not have to make one for each subject and can just import the appropriate content for the glossary (will talk to the group about this)
# Will need to look into finding a way to import math typesetting for this



    # Start an instance of the application
    

class GlossaryWidget(QtWidgets.QWidget):
    def __init__(self, subject, parent=None):
        super().__init__(parent)
        self.subject = subject
        self.listOfContentBoxes = []

        # Main layout
        mainLayout = QtWidgets.QVBoxLayout(self)

        # Scroll area
        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        mainLayout.addWidget(scroll)

        # Content widget inside scroll area
        content = QtWidgets.QWidget()
        scroll.setWidget(content)

        self.verticalLayout = QtWidgets.QVBoxLayout(content)

        # Search bar
        self.searchBar = QtWidgets.QLineEdit()
        self.searchBar.setPlaceholderText("Search definitions or theorems...")
        self.searchBar.setClearButtonEnabled(True)
        self.verticalLayout.addWidget(self.searchBar)

        # Add content boxes
        self.definitions = outputContent(self.subject)
        for word in self.definitions:
            contentBox = CollapsibleBox(word)
            self.verticalLayout.addWidget(contentBox)
            layout = QtWidgets.QVBoxLayout()

            contentLabel = QtWidgets.QLabel(self.definitions[word])
            contentLabel.setStyleSheet("font-size: 10pt;")
            contentLabel.setWordWrap(True)
            contentLabel.setSizePolicy(
                QtWidgets.QSizePolicy.Expanding,
                QtWidgets.QSizePolicy.Preferred)
            layout.addWidget(contentLabel)

            contentBox.setContentLayout(layout)
            self.listOfContentBoxes.append(contentBox)

        self.searchBar.textChanged.connect(self.filterContentBoxes)
        self.verticalLayout.addStretch()

    def filterContentBoxes(self, text):
        text = text.lower().strip()
        for box in self.listOfContentBoxes:
            title = box.toggleButton.text().lower()
            content = box.contentArea.widget().findChild(QtWidgets.QLabel).text().lower()
            if text in title or text in content or text == "":
                box.show()
            else:
                box.hide()


# Stuff here to use the wigdet for testing purposes
#application = QtWidgets.QApplication(sys.argv)
#glossaryWindow,listOfContentBoxes=glossary("GT")
#glossaryWindow.resize(1024, 768)
#glossaryWindow.show()

#sys.exit(application.exec_())
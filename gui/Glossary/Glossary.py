from PyQt5 import QtCore, QtGui, QtWidgets
from CollapsibleBox import *
from Content import *

# Ideally, we do not have to make one for each subject and can just import the appropriate content for the glossary (will talk to the group about this)
# Will need to look into finding a way to import math typesetting for this


if __name__ == "__main__":
    import sys
    
    # Start an instance of the application
    application = QtWidgets.QApplication(sys.argv)
    def glossary(subject):
        # Make window for the application with widgets
        window = QtWidgets.QMainWindow()
        window.setCentralWidget(QtWidgets.QWidget())

        # Make a dockable widget
        dockableWidget = QtWidgets.QDockWidget("Glossary")
        window.addDockWidget(QtCore.Qt.LeftDockWidgetArea, dockableWidget)

        # Scroll bar
        scroll = QtWidgets.QScrollArea()
        dockableWidget.setWidget(scroll)
        content = QtWidgets.QWidget()
        scroll.setWidget(content)
        scroll.setWidgetResizable(True)


        verticalLayout = QtWidgets.QVBoxLayout(content)

        searchBar = QtWidgets.QLineEdit()
        searchBar.setPlaceholderText("Search definitions or theorems...")
        searchBar.setClearButtonEnabled(True)
        verticalLayout.addWidget(searchBar)
        
        listOfContentBoxes = []

        # Flag for the glossary to differentiate bewteen subjects (this is temporary)
        definitions=outputContent(subject)

        # Adding content to the glossary
        for word in definitions:
            contentBox = CollapsibleBox(word)
            verticalLayout.addWidget(contentBox)
            layout = QtWidgets.QVBoxLayout()

            contentLabel = QtWidgets.QLabel(definitions[word])
            contentLabel.setStyleSheet("font-size: 10pt;")
            contentLabel.setWordWrap(True)
            contentLabel.setSizePolicy(
                QtWidgets.QSizePolicy.Expanding,
                QtWidgets.QSizePolicy.Preferred)
            
            layout.addWidget(contentLabel)

            contentBox.setContentLayout(layout)
            listOfContentBoxes.append(contentBox)
        
        searchBar.textChanged.connect(filterContentBoxes)

        verticalLayout.addStretch()

        return window,listOfContentBoxes

    # Function to filter through content boxes for the search bar
    def filterContentBoxes(text):
        text = text.lower().strip()
        
        for box in listOfContentBoxes:
            title = box.toggleButton.text().lower()
            content = box.contentArea.widget().findChild(QtWidgets.QLabel).text().lower()

            # Show the queried item if the search text is in the title or content
            if text in title or text in content or text == "":
                box.show()
            else:
                box.hide()

    
    # Stuff here to use the wigdet for testing purposes
    glossaryWindow,listOfContentBoxes=glossary("GT")
    glossaryWindow.resize(1024, 768)
    glossaryWindow.show()
    
    sys.exit(application.exec_())
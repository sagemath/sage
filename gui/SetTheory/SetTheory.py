import sys
import tempfile
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QPushButton, QMessageBox, QLabel, QWidget,
    QVBoxLayout, QHBoxLayout, QLineEdit, QSizePolicy, QComboBox
)
from PyQt5.QtGui import QPixmap, QIcon
from PyQt5.QtCore import QSize, Qt
from Glossary.Glossary import GlossaryWidget
# from sage.all import set, Set


class SetTheoryTab(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.glossaryWindow = None

        
        self.leftLabel = QLabel("Set A")
        self.leftLabel.setStyleSheet("font-weight: bold;")
        self.leftTextbox = QLineEdit()
        self.leftTextbox.setPlaceholderText("Left set")

        leftLayout = QVBoxLayout()
        leftLayout.addWidget(self.leftLabel)
        leftLayout.addWidget(self.leftTextbox)
        leftLayout.addStretch()

       
        self.rightLabel = QLabel("Proposed partition of A")
        self.rightLabel.setStyleSheet("font-weight: bold;")
        self.rightBoxesLayout = QVBoxLayout()
        self.rightBoxes = []
        self.maxBoxes = 5

        
        firstRow = QHBoxLayout()
        firstBox = QLineEdit()
        firstBox.setPlaceholderText("Right set 1")
        self.rightBoxes.append(firstBox)

        self.addButton = QPushButton("+")
        self.addButton.setFixedWidth(30)
        self.addButton.clicked.connect(self.addTextBox)

        firstRow.addWidget(firstBox)
        firstRow.addWidget(self.addButton)
        self.rightBoxesLayout.addLayout(firstRow)

        
        rightLayout = QVBoxLayout()
        rightLayout.addWidget(self.rightLabel)
        rightLayout.addLayout(self.rightBoxesLayout)
        rightLayout.addStretch()

        
        mainLayout = QHBoxLayout()
        mainLayout.addLayout(leftLayout)
        mainLayout.addSpacing(20)
        mainLayout.addLayout(rightLayout)
        mainLayout.addStretch()

        
        self.setABox = QLineEdit()
        self.setABox.setPlaceholderText("Set A")
        self.setABox.setMaximumWidth(200)

        self.operationDropdown = QComboBox()
        self.operationDropdown.addItems([
                "∪  (Union)",
                "∩  (Intersection)",
                "\  (Difference)",
                "Δ  (Symmetric Difference)",
                "⊆  (Subset)",
                "=  (Equality)"
            ])

        self.operationDropdown.setMaximumWidth(130)
        self.operationDropdown.currentIndexChanged.connect(self.updateEqualButton)

        self.setBBox = QLineEdit()
        self.setBBox.setPlaceholderText("Set B")
        self.setBBox.setMaximumWidth(200)

        self.equalButton = QPushButton("=")
        self.equalButton.setFixedWidth(40)
        self.equalButton.clicked.connect(self.computeOperation)

        self.resultLabel = QLabel("")
        self.resultLabel.setStyleSheet("font-weight: bold; padding-left: 8px;")

        operationRow = QHBoxLayout()
        operationRow.addWidget(self.setABox)
        operationRow.addWidget(self.operationDropdown)
        operationRow.addWidget(self.setBBox)
        operationRow.addWidget(self.equalButton)
        operationRow.addWidget(self.resultLabel)

        
        self.glossaryButton = QPushButton("Glossary")
        self.glossaryButton.clicked.connect(self.showGlossary)
        self.glossaryButton.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)

        
        outerLayout = QVBoxLayout()
        outerLayout.addLayout(mainLayout)
        outerLayout.addSpacing(20)
        outerLayout.addLayout(operationRow)
        outerLayout.addSpacing(20)

        outerLayout.addStretch()
        hbox = QHBoxLayout()
        hbox.addStretch()
        hbox.addWidget(self.glossaryButton)
        hbox.addStretch()
        outerLayout.addLayout(hbox)
        self.setLayout(outerLayout)

    def addTextBox(self):
        """Add another right text box, up to maxBoxes."""
        if len(self.rightBoxes) >= self.maxBoxes:
            return
        newBox = QLineEdit()
        newBox.setPlaceholderText(f"Right set {len(self.rightBoxes) + 1}")
        self.rightBoxes.append(newBox)
        self.rightBoxesLayout.addWidget(newBox)

    def updateEqualButton(self):
        """Update the equal button symbol depending on the operation."""
        op = self.operationDropdown.currentText()
        if "⊆" in op:
            self.equalButton.setText("≡")
        elif "=" in op:
            self.equalButton.setText("≡")
        else:
            self.equalButton.setText("=")

    def computeOperation(self):
        """Placeholder for future set operation logic."""
        # TODO: Implement set operation logic here
        self.resultLabel.setText("Result will appear here")

    def showGlossary(self):
        """Show glossary window."""
        if self.glossaryWindow is None:
            self.glossaryWindow = QMainWindow()
            self.glossaryWindow.setWindowTitle("Set Theory Glossary")
            self.glossaryWindow.setCentralWidget(GlossaryWidget("ST"))

        self.glossaryWindow.show()
        self.glossaryWindow.raise_()
        self.glossaryWindow.activateWindow()



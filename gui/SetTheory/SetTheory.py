from PyQt5.QtWidgets import (
    QWidget, QMainWindow, QPushButton, QMessageBox, QLabel,
    QVBoxLayout, QHBoxLayout, QLineEdit, QSizePolicy, QComboBox
)
from Glossary.Glossary import GlossaryWidget
from sage.all import Set  # ← restored Sage Set


class SetTheoryTab(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.glossaryWindow = None
        self.maxBoxes = 8
        self.rightBoxes = []  
        
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

        
        firstRow = QHBoxLayout()
        firstBox = QLineEdit()
        firstBox.setPlaceholderText("Right set 1")

        self.addButton = QPushButton("Add")
        self.addButton.setFixedWidth(70)
        self.addButton.clicked.connect(self.addTextBox)

        firstRow.addWidget(firstBox)
        firstRow.addWidget(self.addButton)

        self.rightBoxesLayout.addLayout(firstRow)
        self.rightBoxes.append((firstRow, firstBox, None))  

        rightLayout = QVBoxLayout()
        rightLayout.addWidget(self.rightLabel)
        rightLayout.addLayout(self.rightBoxesLayout)
        rightLayout.addStretch()

        #
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
            "\\  (Difference)",
            "Δ  (Symmetric Difference)"
        ])
        self.operationDropdown.setMaximumWidth(130)

        self.setBBox = QLineEdit()
        self.setBBox.setPlaceholderText("Set B")
        self.setBBox.setMaximumWidth(200)

        self.equalButton = QPushButton("=")
        self.equalButton.setFixedWidth(40)
        self.equalButton.clicked.connect(self.performSetOperation)

        self.resultLabel = QLabel("")
        self.resultLabel.setStyleSheet("font-weight: bold; padding-left: 8px;")

        operationRow = QHBoxLayout()
        operationRow.addWidget(self.setABox)
        operationRow.addWidget(self.operationDropdown)
        operationRow.addWidget(self.setBBox)
        operationRow.addWidget(self.equalButton)
        operationRow.addWidget(self.resultLabel)

        
        self.logicalDropdown = QComboBox()
        self.logicalDropdown.addItems([
            "⊆  (Subset)",
            "=  (Equality)"
        ])
        self.logicalDropdown.setMaximumWidth(130)

        self.logicalABox = QLineEdit()
        self.logicalABox.setPlaceholderText("Set A")
        self.logicalABox.setMaximumWidth(200)

        self.logicalBBox = QLineEdit()
        self.logicalBBox.setPlaceholderText("Set B")
        self.logicalBBox.setMaximumWidth(200)

        self.logicalEqualButton = QPushButton("≡")
        self.logicalEqualButton.setFixedWidth(40)
        self.logicalEqualButton.clicked.connect(self.performLogicalOperation)

        self.logicalResultLabel = QLabel("")
        self.logicalResultLabel.setStyleSheet("font-weight: bold; padding-left: 8px;")

        logicalRow = QHBoxLayout()
        logicalRow.addWidget(self.logicalABox)
        logicalRow.addWidget(self.logicalDropdown)
        logicalRow.addWidget(self.logicalBBox)
        logicalRow.addWidget(self.logicalEqualButton)
        logicalRow.addWidget(self.logicalResultLabel)

        
        self.glossaryButton = QPushButton("Glossary")
        self.glossaryButton.clicked.connect(self.showGlossary)
        self.glossaryButton.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)

        
        outerLayout = QVBoxLayout()
        outerLayout.addLayout(mainLayout)
        outerLayout.addSpacing(20)
        outerLayout.addLayout(operationRow)
        outerLayout.addLayout(logicalRow)
        outerLayout.addSpacing(20)
        outerLayout.addStretch()

        self.glossaryButton.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        outerLayout.addWidget(self.glossaryButton)

        self.setLayout(outerLayout)


       
        self.checkPartitionButton = QPushButton("Check Partition")
        self.checkPartitionButton.setFixedWidth(150)
        self.checkPartitionButton.clicked.connect(self.onCheckPartition)
        rightLayout.addWidget(self.checkPartitionButton)


    def addTextBox(self):
        """Add another right text box row, up to maxBoxes."""
        if len(self.rightBoxes) >= self.maxBoxes:
            QMessageBox.warning(self, "Limit Reached", f"Maximum of {self.maxBoxes} boxes allowed.")
            return

        row = QHBoxLayout()
        newBox = QLineEdit()
        newBox.setPlaceholderText(f"Right set {len(self.rightBoxes) + 1}")

        deleteButton = QPushButton("Remove")
        deleteButton.setFixedWidth(70)
        deleteButton.clicked.connect(lambda _, b=newBox: self.deleteTextBox(b))

        row.addWidget(newBox)
        row.addWidget(deleteButton)

        self.rightBoxesLayout.addLayout(row)
        self.rightBoxes.append((row, newBox, deleteButton))

    def deleteTextBox(self, boxToRemove):
        for i, (layout, box, deleteButton) in enumerate(self.rightBoxes):
            if box == boxToRemove:
                self.clearLayout(layout)
                self.rightBoxesLayout.removeItem(layout)
                self.rightBoxes.pop(i)
                break

        for idx, (_, box, _) in enumerate(self.rightBoxes):
            box.setPlaceholderText(f"Right set {idx + 1}")

    # Helper function to delete widgets in a layout
    def clearLayout(self, layout):
        while layout.count():
            item = layout.takeAt(0)
            widget = item.widget()
            if widget:
                widget.deleteLater()

    def checkIfPartition(self, setFamily, setX):
        # Use Sage Set operations. Note: Set.union returns a new Set.
        union_all = Set([])
        for part in setFamily:
            # Ensure part is a Sage Set (assumed elsewhere) and union it
            union_all = union_all.union(part)

        # Check that the union equals the target set
        if union_all != setX:
            return False

        # No empty parts (use Sage's cardinality)
        for part in setFamily:
            try:
                if part.cardinality() == 0:
                    return False
            except Exception:
                # Fallback: try Python len() if cardinality() is unavailable
                if len(part) == 0:
                    return False

        # Pairwise disjoint: intersection should have cardinality 0
        n = len(setFamily)
        for i in range(n):
            for j in range(i + 1, n):
                try:
                    if setFamily[i].intersection(setFamily[j]).cardinality() != 0:
                        return False
                except Exception:
                    # Fallback: convert to Python sets and test
                    if len(set(setFamily[i]).intersection(set(setFamily[j]))) != 0:
                        return False

        return True
    
    def onCheckPartition(self):
        setAText = self.leftTextbox.text().strip()
        if not setAText:
            QMessageBox.warning(self, "Error", "Please enter Set A.")
            return

        try:
            setA = self.parseSetInput(setAText)
        except Exception:
            QMessageBox.warning(self, "Error", "Invalid Set A format.")
            return

        # Collect non-empty right sets
        setFamily = []
        for _, box, _ in self.rightBoxes:
            text = box.text().strip()
            if text:
                try:
                    s = self.parseSetInput(text)
                    setFamily.append(s)
                except Exception:
                    QMessageBox.warning(self, "Error", f"Invalid set format: {text}")
                    return

        if not setFamily:
            QMessageBox.warning(self, "Error", "Please enter at least one subset.")
            return

        # Check partition
        isPartition = self.checkIfPartition(setFamily, setA)
        if isPartition:
            QMessageBox.information(self, "Result", "The sets form a valid partition of A.")
        else:
            QMessageBox.warning(self, "Result", "The sets do not form a partition of A.")


    # Solely to parse set input strings into Sage Sets
    def parseSetInput(self, inputStr):
        stringStripped = inputStr.strip()
        if stringStripped.startswith("{") and stringStripped.endswith("}"):
            stringStripped = stringStripped[1:-1]  # remove braces
        elements = [item.strip() for item in stringStripped.split(",") if item.strip()]
        return Set(elements)

    def performSetOperation(self):
        setA = self.parseSetInput(self.setABox.text())
        setB = self.parseSetInput(self.setBBox.text())
        operation = self.operationDropdown.currentText()

        if operation == "∪  (Union)":
            self.resultLabel.setText(str(setA.union(setB)))
        elif operation == "∩  (Intersection)":
            self.resultLabel.setText(str(setA.intersection(setB)))
        elif operation == "\\  (Difference)":
            self.resultLabel.setText(str(setA.difference(setB)))
        elif operation == "Δ  (Symmetric Difference)":
            self.resultLabel.setText(str(setA.difference(setB).union(setB.difference(setA))))

    def performLogicalOperation(self):
        setA = self.parseSetInput(self.logicalABox.text())
        setB = self.parseSetInput(self.logicalBBox.text())
        operation = self.logicalDropdown.currentText()

        if operation == "⊆  (Subset)":
            self.logicalResultLabel.setText(str(setA.issubset(setB)))
        elif operation == "=  (Equality)":
            self.logicalResultLabel.setText(str(setA == setB))

    def showGlossary(self):
        if self.glossaryWindow is None:
            self.glossaryWindow = QMainWindow()
            self.glossaryWindow.setWindowTitle("Set Theory Glossary")
            self.glossaryWindow.setCentralWidget(GlossaryWidget("ST"))
        self.glossaryWindow.show()
        self.glossaryWindow.raise_()
        self.glossaryWindow.activateWindow()

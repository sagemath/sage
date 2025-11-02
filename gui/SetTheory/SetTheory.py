from PyQt5.QtWidgets import (
    QWidget, QMainWindow, QPushButton, QMessageBox, QLabel,
    QVBoxLayout, QHBoxLayout, QLineEdit, QSizePolicy, QComboBox
)
from Glossary.Glossary import GlossaryWidget
#from sage.all import Set


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

        mainLayout = QHBoxLayout()
        mainLayout.addLayout(leftLayout)
        mainLayout.addSpacing(20)
        mainLayout.addLayout(rightLayout)
        mainLayout.addStretch()

        # Operation rows
        unionRow, self.unionABox, self.unionButton, self.unionBBox, self.unionResult = self.makeOperationRow("∪", "Union")
        interRow, self.interABox, self.intersectionButton, self.interBBox, self.interResult = self.makeOperationRow("∩", "Intersection")
        diffRow, self.diffABox, self.differenceButton, self.diffBBox, self.diffResult = self.makeOperationRow("\\", "Difference")
        symRow, self.symABox, self.symmetricDiffButton, self.symBBox, self.symResult = self.makeOperationRow("Δ", "Symmetric Difference")

        self.unionButton.clicked.connect(lambda: self.performSetOperation("union"))
        self.intersectionButton.clicked.connect(lambda: self.performSetOperation("intersection"))
        self.differenceButton.clicked.connect(lambda: self.performSetOperation("difference"))
        self.symmetricDiffButton.clicked.connect(lambda: self.performSetOperation("symmetric_difference"))

        operationLayout = QVBoxLayout()
        operationLayout.setSpacing(4)
        operationLayout.setContentsMargins(0, 0, 0, 0)
        operationLayout.addLayout(unionRow)
        operationLayout.addLayout(interRow)
        operationLayout.addLayout(diffRow)
        operationLayout.addLayout(symRow)

        self.resultLabel = QLabel("")
        self.resultLabel.setStyleSheet("font-weight: bold; margin-top: 10px;")
        operationLayout.addWidget(self.resultLabel)

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
        outerLayout.addLayout(operationLayout)
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


    def makeOperationRow(self, symbol, label):
        leftBox = QLineEdit()
        leftBox.setPlaceholderText("Set A")
        leftBox.setMaximumWidth(200)
        rightBox = QLineEdit()
        rightBox.setPlaceholderText("Set B")
        rightBox.setMaximumWidth(200)
        button = QPushButton(f"{symbol}  ({label})")
        eqLabel = QLabel("=")
        resultLabel = QLabel("")
        resultLabel.setStyleSheet("font-weight: bold; padding-left: 8px;")
        button.setFixedWidth(150)
        row = QHBoxLayout()
        row.setSpacing(8)
        row.setContentsMargins(0, 0, 0, 0)
        row.addWidget(leftBox)
        row.addWidget(button)   
        row.addWidget(rightBox)
        row.addWidget(eqLabel)
        row.addWidget(resultLabel)
        return row, leftBox, button, rightBox, resultLabel

    def performSetOperation(self, operation):
        try:
            if operation == "union":
                setA = self.parseSetInput(self.unionABox.text())
                setB = self.parseSetInput(self.unionBBox.text())
                self.unionResult.setText(str(setA.union(setB)))
            elif operation == "intersection":
                setA = self.parseSetInput(self.interABox.text())
                setB = self.parseSetInput(self.interBBox.text())
                self.interResult.setText(str(setA.intersection(setB)))
            elif operation == "difference":
                setA = self.parseSetInput(self.diffABox.text())
                setB = self.parseSetInput(self.diffBBox.text())
                self.diffResult.setText(str(setA.difference(setB)))
            elif operation == "symmetric_difference":
                setA = self.parseSetInput(self.symABox.text())
                setB = self.parseSetInput(self.symBBox.text())
                self.symResult.setText(str(setA.difference(setB).union(setB.difference(setA))))
        except Exception:
            QMessageBox.warning(self, "Error", "Invalid set format.")

    def addTextBox(self):
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

    def clearLayout(self, layout):
        while layout.count():
            item = layout.takeAt(0)
            widget = item.widget()
            if widget:
                widget.deleteLater()

    def checkIfPartition(self, setFamily, setX):
        union_all = Set([])
        for part in setFamily:
            union_all = union_all.union(part)
        if union_all != setX:
            return False
        for part in setFamily:
            if part.cardinality() == 0:
                return False
        n = len(setFamily)
        for i in range(n):
            for j in range(i + 1, n):
                if setFamily[i].intersection(setFamily[j]).cardinality() != 0:
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
        isPartition = self.checkIfPartition(setFamily, setA)
        if isPartition:
            QMessageBox.information(self, "Result", "The sets form a valid partition of A.")
        else:
            QMessageBox.warning(self, "Result", "The sets do not form a partition of A.")

    def parseSetInput(self, inputStr):
        stringStripped = inputStr.strip()
        if stringStripped.startswith("{") and stringStripped.endswith("}"):
            stringStripped = stringStripped[1:-1]
        elements = [item.strip() for item in stringStripped.split(",") if item.strip()]
        return Set(elements)

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

from PyQt5.QtWidgets import (
    QWidget, QMainWindow, QPushButton, QMessageBox, QLabel,
    QVBoxLayout, QHBoxLayout, QLineEdit, QSizePolicy, QComboBox,
    QGroupBox, QTextEdit, QDialog
)
from PyQt5.QtCore import QUrl
from PyQt5.QtWebEngineWidgets import QWebEngineView
from Glossary.Glossary import GlossaryWidget
from SetTheory.QuestionGenerator import generateQuestion
from sage.all import Set

class SetTheoryLearningTab(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        # Glossary initialization
        self.glossaryWindow = None

        # operation initalizations
        self.unionResult = ""
        self.interResult = ""
        self.subsetResult = ""
        self.symResult = ""
        self.partitionResult = ""
        self.diffResult = ""
        self.complimentResultOne = ""
        self.complimentResultTwo = ""

        # miss countere initialization
        self.missCounter = 0

        # Main layout
        main_layout = QVBoxLayout(self)

        # Top box: display area for a string (read-only, can show multi-line)
        self.display_box = QGroupBox("Display")
        display_layout = QVBoxLayout(self.display_box)
        self.display_view = QTextEdit()
        self.display_view.setReadOnly(True)
        self.display_view.setPlaceholderText("Displayed text will appear here...")
        self.display_view.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        display_layout.addWidget(self.display_view)
        self.display_box.setLayout(display_layout)

        # Bottom box: input baseline with text input and a submit button
        self.input_box = QGroupBox("Input")
        input_layout = QHBoxLayout(self.input_box)
        self.input_field = QLineEdit()
        self.input_field.setPlaceholderText("Type here and press Enter or Submit")
        self.input_field.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.submit_btn = QPushButton("Submit")
        self.submit_btn.setDefault(True)

        input_layout.addWidget(self.input_field)
        input_layout.addWidget(self.submit_btn)

        # Glossary button
        self.glossaryButton = QPushButton("Glossary")
        self.glossaryButton.clicked.connect(self.showGlossary)
        self.glossaryButton.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        input_layout.addWidget(self.glossaryButton)

        self.input_box.setLayout(input_layout)

        # Assemble
        main_layout.addWidget(self.display_box, stretch=3)
        main_layout.addWidget(self.input_box, stretch=0)
        self.setLayout(main_layout)

        # Connections
        self.submit_btn.clicked.connect(self._on_submit)
        self.input_field.returnPressed.connect(self._on_submit)

        # Initial example text
        self.setQuestion()

    def _on_submit(self):
        correct = "false"
        text = self.input_field.text().replace(" ", "").lower()
        correct = self.answerCheck(text)
        if not correct:
            if self.missCounter == 0:
                QMessageBox.warning(self, "Incorrect", "That answer is incorrect.")
                self.missCounter += 1
            else:
                self.incorrectBuzzer()
                self.missCounter += 1
        else:
            QMessageBox.information(self, "Correct!", "That answer is correct! Heres another question.")
            self.missCounter = 0
            self.cleanAnswers()
            self.setQuestion()

        

    def setQuestion(self):
        qNum, question, setsDict = generateQuestion()
        self.display_view.setPlainText(question)
        self.performSetOperation(qNum, setsDict)

    def showGlossary(self):
        if self.glossaryWindow is None:
            self.glossaryWindow = QMainWindow()
            self.glossaryWindow.setWindowTitle("Set Theory Glossary")
            self.glossaryWindow.setCentralWidget(GlossaryWidget("ST"))
        self.glossaryWindow.show()
        self.glossaryWindow.raise_()
        self.glossaryWindow.activateWindow()

    def performSetOperation(self, opNum, setsDict):
        try:
            if opNum == 1: # union
                setA = self.parseSetInput(setsDict[0])
                setB = self.parseSetInput(setsDict[1])
                self.unionResult = (str(setA.union(setB)))
            elif opNum == 2: # intersection
                setA = self.parseSetInput(setsDict[0])
                setB = self.parseSetInput(setsDict[1])
                self.interResult = (str(setA.intersection(setB)))
            elif opNum == 3: # subset
                self.subsetResult = ("true")
            elif opNum == 4: # symmetric difference
                setA = self.parseSetInput(setsDict[0])
                setB = self.parseSetInput(setsDict[1])
                self.symResult = (str(setA.difference(setB).union(setB.difference(setA))))
            elif opNum == 5: # partition
                setA = self.parseSetInput(setsDict[0])
                setB = self.parseSetInput(setsDict[1])
                setC = self.parseSetInput(setsDict[2])
                setD = self.parseSetInput(setsDict[3])
                self.partitionResult = (self.checkIfPartition([setB, setC, setD], setA))
            elif opNum == 6: # difference
                setA = self.parseSetInput(setsDict[0])
                setB = self.parseSetInput(setsDict[1])
                self.diffResult = (str(setA.difference(setB)))
            elif opNum == 7: # B compliment and C compliment
                setA = self.parseSetInput(setsDict[0])
                setB = self.parseSetInput(setsDict[1])
                setC = self.parseSetInput(setsDict[2])
                self.complimentResultOne = (str(setA.difference(setB)))
                self.complimentResultTwo = (str(setA.difference(setC)))
        except Exception:
            QMessageBox.warning(self, "Error", "Invalid set format.")

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
    
    def parseSetInput(self, inputStr):
        stringStripped = inputStr.strip()
        if stringStripped.startswith("{") and stringStripped.endswith("}"):
            stringStripped = stringStripped[1:-1]
        elements = [item.strip() for item in stringStripped.split(",") if item.strip()]
        return Set(elements)
    
    def answerCheck(self, text):
        if text == "true" or text == "false":
            return self.boolAnswerCheck(text)
        return self.setAnswerCheck(text)

    def boolAnswerCheck(self, answer):
        if self.subsetResult == "":
            return answer == self.partitionResult
        else:
            return answer == self.subsetResult

    def setAnswerCheck(self, text):
        if self.unionResult != "":
            answer = self.parseSetInput(text)
            return answer == self.unionResult
        elif self.interResult != "":
            answer = self.parseSetInput(text)
            return answer == self.interResult
        elif self.symResult != "":
            answer = self.parseSetInput(text)
            return answer == self.symResult
        elif self.diffResult != "":
            answer = self.parseSetInput(text)
            return answer == self.diffResult
        else:
            answerList = text.split(';')
            answerOne = self.parseSetInput(answerList[0])
            answerTwo = self.parseSetInput(answerList[1])
            return (answerOne == self.complimentResultOne and answerTwo == self.complimentResultTwo)

    def incorrectBuzzer(self):
        hint = "Please review the glossary definition of "
        if self.missCounter < 3:
            glossaryRec = self.fetchQuestion()
            hint = hint + glossaryRec + "."
            QMessageBox.warning(self, "Incorrect", hint)
        else:
            self.showYoutubeHint("https://www.youtube.com/watch?v=iTZfATpm0Yk")
    
    def fetchQuestion(self):
        if self.unionResult != "":
            return "Unions"
        elif self.interResult != "":
            return "Intersections"
        elif self.symResult != "":
            return "Symmetric Differences"
        elif self.diffResult != "":
            return "Difference"
        elif self.subsetResult != "":
            return "Subsets"
        elif self.partitionResult != "":
            return "Partitions"
        else:
            return "Compliments"
        
    def showYoutubeHint(self, youtube_url):
        #shows pop up of youtube video
        dialog = QDialog(self)
        dialog.setWindowTitle("INCORRECT\n Here's a video to help!")

        layout = QVBoxLayout(dialog)

        view = QWebEngineView()
        layout.addWidget(view)

        view.load(QUrl(youtube_url))

        dialog.resize(800, 450)
        dialog.exec_()
    
    def cleanAnswers(self):
        self.unionResult = ""
        self.interResult = ""
        self.subsetResult = ""
        self.symResult = ""
        self.partitionResult = ""
        self.diffResult = ""
        self.complimentResultOne = ""
        self.complimentResultTwo = ""
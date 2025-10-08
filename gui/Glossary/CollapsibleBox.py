from PyQt5 import QtCore, QtGui, QtWidgets

# Note: PyQt does NOT have a default widget for this, so I had help from https://stackoverflow.com/questions/52615115/how-to-create-collapsible-box-in-pyqt and used it as a base
# Goal: clicking a glossary button will have this window pop out
# Want to use this as a general class for the glossary

class CollapsibleBox(QtWidgets.QWidget):
    def __init__(self, title="", parent=None):
        super(CollapsibleBox, self).__init__(parent)

        self.labelContent = QtWidgets.QLabel(
            wordWrap=True,
            alignment=QtCore.Qt.AlignLeft|QtCore.Qt.AlignTop)

        self.toggleButton = QtWidgets.QToolButton(
            text=title, checkable=True, checked=False
        )
        self.toggleButton.setStyleSheet("QToolButton { border: none; }")
        self.toggleButton.setToolButtonStyle(
            QtCore.Qt.ToolButtonTextBesideIcon
        )
        self.toggleButton.setArrowType(QtCore.Qt.RightArrow)
        self.toggleButton.pressed.connect(self.pressButton)

        self.toggleAnimation = QtCore.QParallelAnimationGroup(self)


        
        self.contentArea = QtWidgets.QScrollArea(
            maximumHeight=0, minimumHeight=0
        )
        self.contentArea.setWidgetResizable(True)
        self.contentArea.setSizePolicy(
            QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed
        )
        self.contentArea.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.toggleButton.setStyleSheet("QToolButton { font-size: 10pt; }")

        layout = QtWidgets.QVBoxLayout(self)
        layout.setSpacing(0)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.toggleButton)
        layout.addWidget(self.contentArea)

        self.toggleAnimation.addAnimation(
            QtCore.QPropertyAnimation(self, b"minimumHeight")
        )
        self.toggleAnimation.addAnimation(
            QtCore.QPropertyAnimation(self, b"maximumHeight")
        )
        self.toggleAnimation.addAnimation(
            QtCore.QPropertyAnimation(self.contentArea, b"maximumHeight")
        )
        
    # Function for the buttons
    @QtCore.pyqtSlot()
    def pressButton(self):
        checked = self.toggleButton.isChecked()
        self.toggleButton.setArrowType(
            QtCore.Qt.DownArrow if not checked else QtCore.Qt.RightArrow
        )
        self.toggleAnimation.setDirection(
            QtCore.QAbstractAnimation.Forward
            if not checked
            else QtCore.QAbstractAnimation.Backward
        )
        self.toggleAnimation.start()

    # Formatting the content on the labels
    def setContentLayout(self, layout):
        contentWidget = QtWidgets.QWidget()
        contentWidget.setLayout(layout)
        self.contentArea.setWidget(contentWidget)
        self.contentArea.setWidgetResizable(True)
        collapsedHeight = (
            self.sizeHint().height() - self.contentArea.maximumHeight()
        )
        contentHeight = layout.sizeHint().height()
        for i in range(self.toggleAnimation.animationCount()):
            animation = self.toggleAnimation.animationAt(i)
            animation.setDuration(500)
            animation.setStartValue(collapsedHeight)
            animation.setEndValue(collapsedHeight + contentHeight)

        contentAnimation = self.toggleAnimation.animationAt(
            self.toggleAnimation.animationCount() - 1
        )
        contentAnimation.setDuration(500)
        contentAnimation.setStartValue(0)
        contentAnimation.setEndValue(contentHeight)
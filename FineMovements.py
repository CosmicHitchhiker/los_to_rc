from PySide6.QtCore import Slot, Signal
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QFormLayout,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QFileDialog,
    QDoubleSpinBox,
)
from PySide6.QtGui import QShortcut, QKeySequence


class FineMovementsDialog(QDialog):
    move_ra = Signal(float)
    move_dec = Signal(float)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.ra_step = QDoubleSpinBox()
        self.dec_step = QDoubleSpinBox()
        self.ra_step.setDecimals(2)
        self.dec_step.setDecimals(2)
        self.ra_step.setValue(0.01)
        self.dec_step.setValue(0.01)
        self.ra_step.setSingleStep(0.01)
        self.dec_step.setSingleStep(0.01)

        self.quit_button = QPushButton(text="Quit")
        self.move_button = QPushButton(text="Move")
        self.move_button.setCheckable(True)

        layout = QFormLayout()
        layout.addRow("RA Step (hourangle sec):", self.ra_step)
        layout.addRow("Dec Step (arcsec):", self.dec_step)
        layout.addRow(self.move_button, self.quit_button)
        self.setLayout(layout)

        # Keyboard Shortcuts (for fine movements)
        self.left_shortcut = QShortcut(QKeySequence('Left'), self)
        self.right_shortcut = QShortcut(QKeySequence('Right'), self)
        self.up_shortcut = QShortcut(QKeySequence('Up'), self)
        self.down_shortcut = QShortcut(QKeySequence('Down'), self)
        self.fine_shortcuts = [self.left_shortcut,
                               self.right_shortcut,
                               self.up_shortcut,
                               self.down_shortcut]

        self.move_button.toggled.connect(self.activate_movements)
        self.quit_button.clicked.connect(self.close)

    def activate_movements(self, activate: bool):
        if activate:
            self.left_shortcut.activated.connect(lambda: self.move_ra.emit(self.ra_step.value()))
            self.right_shortcut.activated.connect(lambda: self.move_ra.emit(-self.ra_step.value()))
            self.up_shortcut.activated.connect(lambda: self.move_dec.emit(self.dec_step.value()))
            self.down_shortcut.activated.connect(lambda: self.move_dec.emit(-self.dec_step.value()))
        else:
            self.left_shortcut.activated.disconnect()
            self.right_shortcut.activated.disconnect()
            self.up_shortcut.activated.disconnect()
            self.down_shortcut.activated.disconnect()
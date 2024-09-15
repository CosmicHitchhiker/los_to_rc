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
    QCheckBox,
)
from PySide6.QtGui import QShortcut, QKeySequence
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u


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

        self.pa = 320.5
        self.center = SkyCoord(0*u.deg, 0*u.deg)

        self.quit_button = QPushButton(text="Quit")
        self.slit_checkbox = QCheckBox()
        self.move_button = QPushButton(text="Move")
        self.move_button.setCheckable(True)

        layout = QFormLayout()
        layout.addRow("RA Step (hourangle sec):", self.ra_step)
        layout.addRow("Dec Step (arcsec):", self.dec_step)
        layout.addRow("Move along PA", self.slit_checkbox)
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
            self.left_shortcut.activated.connect(lambda: self.move_lr())
            self.right_shortcut.activated.connect(lambda: self.move_lr(-1))
            self.up_shortcut.activated.connect(lambda: self.move_ud())
            self.down_shortcut.activated.connect(lambda: self.move_ud(-1))
        else:
            self.left_shortcut.activated.disconnect()
            self.right_shortcut.activated.disconnect()
            self.up_shortcut.activated.disconnect()
            self.down_shortcut.activated.disconnect()

    def move_lr(self, k=1):
        if self.slit_checkbox.isChecked():
            new_coords = self.center.directional_offset_by((self.pa+90)*u.deg,
                                                           k * self.ra_step.value() * u.arcsec)
            dra = (new_coords.ra - self.center.ra).to(u.arcsec).value
            ddec = (new_coords.dec - self.center.dec).to(u.arcsec).value
            self.move_ra.emit(dra)
            self.move_dec.emit(ddec)
        else:
            self.move_ra.emit(k * self.ra_step.value())

    def move_ud(self, k=1):
        if self.slit_checkbox.isChecked():
            new_coords = self.center.directional_offset_by(self.pa*u.deg,
                                                           k * self.dec_step.value() * u.arcsec)
            dra = (new_coords.ra - self.center.ra).to(u.arcsec).value
            ddec = (new_coords.dec - self.center.dec).to(u.arcsec).value
            self.move_ra.emit(dra)
            self.move_dec.emit(ddec)
        else:
            self.move_dec.emit(k * self.dec_step.value())


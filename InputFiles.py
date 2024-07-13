from PySide6.QtCore import Slot, Signal
from PySide6.QtWidgets import (
    QWidget,
    QDialog,
    QDoubleSpinBox,
    QSpinBox,
    QVBoxLayout,
    QHBoxLayout,
    QPushButton,
    QLabel,
    QButtonGroup,
    QTableWidget,
    QTableWidgetItem,
    QFileDialog,
)
from PySide6.QtGui import QColor

from pathlib import Path
from matplotlib import colormaps
import numpy as np


class slitParams:
    def __init__(self, n=0, csv_path=None, is_used=True):
        colors = colormaps['tab20'](np.linspace(0, 1, 20))
        self.colors = colors[[2 * n, 2 * n + 1]]
        self.csv_path = csv_path
        self.is_used = is_used


class InputDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        # block any actions in the parent window
        # noinspection PyUnresolvedReferences
        self.setModal(QDialog.Accepted)

        self.table = QTableWidget()
        self.add_button = QPushButton(text='+')
        self.delete_button = QPushButton(text='Del')
        self.apply_button = QPushButton(text='Apply')

        self.dir = str(Path.home())
        self.data = []

        self.__configure_widgets()
        self.__setLayout()

        self.add_button.clicked.connect(self.add_file)
        self.delete_button.clicked.connect(self.delete_row)

    def __configure_widgets(self):
        self.table.setRowCount(0)
        self.table.setColumnCount(4)
        self.table.setHorizontalHeaderLabels(["Path", "Color 1", "Color 2", "Use"])

    def __setLayout(self):
        self.vlayout = QVBoxLayout()

        button_row = QHBoxLayout()
        button_row.addWidget(self.add_button)
        button_row.addWidget(self.delete_button)
        button_row.addWidget(self.apply_button)

        self.vlayout.addWidget(self.table)
        self.vlayout.addLayout(button_row)
        self.setLayout(self.vlayout)
        w = self.table.width()
        self.setFixedWidth(w)

    @Slot()
    def add_file(self):
        file_path = QFileDialog.getOpenFileName(self, "CSV", self.dir,
                                                "All (*)")[0]
        self.dir = "/".join(file_path.split('/')[:-1])
        n = len(self.data)
        self.data.append(slitParams(n, csv_path=file_path))
        self.table.insertRow(n)
        self.row_from_slitP(n, self.data[n])

    def delete_row(self):
        n = self.table.currentRow()
        self.table.removeRow(n)
        self.data.pop(n)

    def row_from_slitP(self, n: int, params: slitParams):
        path = QTableWidgetItem(params.csv_path)
        c1 = QTableWidgetItem()
        c1.setBackground(QColor(self.rgb_from_color(params.colors[0])))
        c2 = QTableWidgetItem()
        c2.setBackground(QColor(self.rgb_from_color(params.colors[1])))

        self.table.setItem(n, 0, path)
        self.table.setItem(n, 1, c1)
        self.table.setItem(n, 2, c2)

    @staticmethod
    def rgb_from_color(color):
        rgb_array = np.array(color)[:3] * 255
        rgb = '#'
        for c in rgb_array.astype(int):
            col = hex(c)[2:]
            rgb += '0' * (2-len(col)) + col
        return rgb


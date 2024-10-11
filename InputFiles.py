from pathlib import Path

import numpy as np
import pandas as pd
from PySide6.QtCore import Slot
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QFileDialog,
)
from slitParams import slitParams


def check_valid_path(path):
    try:
        pd.read_csv(path)
        return True
    except:
        print('INVALID CSV FILE OR PATH')
        return False


class InputDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        # block any actions in the parent window
        # noinspection PyUnresolvedReferences
        self.setModal(QDialog.Accepted)

        self.table = QTableWidget()
        self.add_button = QPushButton(text='+')
        self.delete_button = QPushButton(text='Del')
        self.ok_button = QPushButton(text='Ok')

        self.dir = str(Path.home())
        self.data = []

        self.__configure_widgets()
        self.__setLayout()

        self.add_button.clicked.connect(lambda: self.add_file(file_path=None))
        self.delete_button.clicked.connect(self.delete_row)
        self.ok_button.clicked.connect(self.accept)
        self.table.itemChanged.connect(self.change_info)

    def __configure_widgets(self):
        self.table.setRowCount(0)
        self.table.setColumnCount(4)
        self.table.setHorizontalHeaderLabels(["Path", "Color 1", "Color 2", "Label"])

    def __setLayout(self):
        self.vlayout = QVBoxLayout()

        button_row = QHBoxLayout()
        button_row.addWidget(self.add_button)
        button_row.addWidget(self.delete_button)
        button_row.addWidget(self.ok_button)

        self.vlayout.addWidget(self.table)
        self.vlayout.addLayout(button_row)
        self.setLayout(self.vlayout)
        w = self.table.width()
        self.setFixedWidth(w)

    @Slot()
    def add_file(self, file_path=None):
        if file_path is None:
            file_path = QFileDialog.getOpenFileName(self, "CSV", self.dir,
                                                    "All (*)")[0]
        if check_valid_path(file_path):
            self.dir = "/".join(file_path.split('/')[:-1])
            n = len(self.data)
            self.data.append(slitParams(csv_path=file_path))
            self.table.insertRow(n)
            self.row_from_slitP(n, self.data[n])

    def delete_row(self):
        n = self.table.currentRow()
        m = self.data[n].m
        self.table.removeRow(n)
        self.data.pop(n)
        slitParams.avaliable_n.insert(0, m)
        slitParams.avaliable_n.sort()
        # for num, dat in enumerate(self.data):
        #     dat.set_n(num)
        #     self.row_from_slitP(num, dat)

    def row_from_slitP(self, n: int, params: slitParams):
        path = QTableWidgetItem(params.csv_path)
        c1 = QTableWidgetItem()
        c1.setBackground(QColor(self.rgb_from_color(params.colors[0])))
        c2 = QTableWidgetItem()
        c2.setBackground(QColor(self.rgb_from_color(params.colors[1])))

        self.table.setItem(n, 0, path)
        self.table.setItem(n, 1, c1)
        self.table.setItem(n, 2, c2)

    @Slot()
    def apply_files(self):
        super().accept()

    @Slot()
    def change_info(self, item):
        if item.column() == 3:
            self.data[item.row()].label = item.text()

    @staticmethod
    def rgb_from_color(color):
        rgb_array = np.array(color)[:3] * 255
        rgb = '#'
        for c in rgb_array.astype(int):
            col = hex(c)[2:]
            rgb += '0' * (2-len(col)) + col
        return rgb


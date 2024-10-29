import requests
import pandas as pd
from io import StringIO

from PySide6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QFileDialog,
)

class HyperledaDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        # block any actions in the parent window
        # noinspection PyUnresolvedReferences
        self.setModal(QDialog.Accepted)

        self.table = QTableWidget()
        self.req_button = QPushButton(text='Request')
        self.select_button = QPushButton(text='Set chosen parameters')
        self.exit_button = QPushButton(text='Exit')

        self.__configure_widgets()
        self.__set_layout()

        self.exit_button.clicked.connect(self.close)

    def __configure_widgets(self):
        self.table.setRowCount(0)
        self.table.setColumnCount(7)
        self.table.setHorizontalHeaderLabels(["Name", "V mag","RA", "DEC", "PA", "i", "Vsys"])

    def __set_layout(self):
        self.vlayout = QVBoxLayout()

        button_row = QHBoxLayout()
        button_row.addWidget(self.req_button)
        button_row.addWidget(self.select_button)
        button_row.addWidget(self.exit_button)

        self.vlayout.addWidget(self.table)
        self.vlayout.addLayout(button_row)
        self.setLayout(self.vlayout)
        w = self.table.width()
        self.setFixedWidth(w)
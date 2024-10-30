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

        # Define the base URL
        self.base_url = "http://atlas.obs-hp.fr/hyperleda/fG.cgi"

        # Define the query parameters
        self.req_params = {
            'n': 'meandata',
            'c': 'o',
            'of': '1,leda,simbad',
            'nra': 'l',
            'nakd': '1',
            'd': 'objname,vt,bt,al2000,de2000,pa,incl,v',
            'sql': 'al2000>(18.2153528-0.05) and al2000<(18.2153528+0.05) and '
                   'de2000>(68.3634444-0.05) and de2000<(68.3634444+0.05)',
            'ob': '',
            'a': 'csv[,]'
        }

        # Changes when image or slit changes
        self.ra_lims = [0, 0]
        self.dec_lims = [0, 0]

        self.__configure_widgets()
        self.__set_layout()

        self.exit_button.clicked.connect(self.close)
        self.req_button.clicked.connect(self.hyperleda_request)
        self.table.itemSelectionChanged.connect(lambda: self.table.selectRow(self.table.currentRow()))

    def __configure_widgets(self):
        self.table.setRowCount(0)
        self.table.setColumnCount(7)
        header_labels = ["Name", "V mag", "B mag", "RA", "DEC", "PA", "i", "Vsys"]
        self.header_dict = {k:v for (v,k) in enumerate(header_labels)}
        self.table.setHorizontalHeaderLabels(header_labels)

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

    def hyperleda_request(self):
        self.table.setRowCount(0)
        response = requests.get(self.base_url, params=self.req_params)
        if response.status_code == 200:
            # Load the CSV data into a pandas DataFrame
            data = pd.read_csv(StringIO(response.text), comment='#')
            print(data.head())  # Display the first few rows
            #TODO: Сначала распарсить датафрейм (подсчитать расстояния до щели, сделать нормальный формат RA и DEC)
            for i in range(len(data)):
                self.table.insertRow(i)
                self.fill_row_from_pandas(i, data.loc[i])
        else:
            print("Failed to retrieve data:", response.status_code)



    def fill_row_from_pandas(self, n, series):
        #Name
        self.table.setItem(n, self.header_dict['Name'], QTableWidgetItem(series['objname']))
        #V mag
        #B mag

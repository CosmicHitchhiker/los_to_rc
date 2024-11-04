import requests
import pandas as pd
from io import StringIO

from PySide6.QtCore import Slot, Signal
from PySide6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QFileDialog,
)
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
from partd.utils import frame


class HyperledaDialog(QDialog):
    params_set = Signal()
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
        self.slit_pos = SkyCoord(0, 0, frame='icrs', unit='deg')

        self.__configure_widgets()
        self.__set_layout()

        self.exit_button.clicked.connect(self.close)
        self.req_button.clicked.connect(self.hyperleda_request)
        self.select_button.clicked.connect(self.params_set.emit)
        self.table.itemSelectionChanged.connect(lambda: self.table.selectRow(self.table.currentRow()))

    def __configure_widgets(self):
        self.table.setRowCount(0)
        self.header_labels = ["Name", "V mag", "B mag", "separation", "RA", "DEC", "PA", "i", "Vsys"]
        # Only labels to show in the table should correspond to header_labels
        self.pandas_labels = ["objname", "vt", "bt", "sep", "RA", "DEC", "pa", "incl", "v"]
        self.header_dict = {k:v for (v,k) in enumerate(self.header_labels)}
        self.table.setColumnCount(len(self.header_labels))
        self.table.setHorizontalHeaderLabels(self.header_labels)

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

    @Slot()
    def hyperleda_request(self):
        self.table.setRowCount(0)
        response = requests.get(self.base_url, params=self.req_params)
        if response.status_code == 200:
            # Load the CSV data into a pandas DataFrame
            try:
                data = pd.read_csv(StringIO(response.text), comment='#')
            except pd.errors.EmptyDataError:
                print("No data found.")
                return
            self.parse_data(data)
            for i in range(len(self.data)):
                self.table.insertRow(i)
                self.fill_row_from_pandas(i, self.data.loc[i])
        else:
            print("Failed to retrieve data:", response.status_code)

    def parse_data(self, data):
        self.data = data
        coords = SkyCoord(data['al2000'], data['de2000'], unit=(u.hourangle, u.deg), frame='icrs')
        self.data['RA'] = coords.ra.to_string(unit=u.hourangle, sep=":")
        self.data['DEC'] = coords.dec.to_string(unit=u.deg, sep=":")
        self.data['sep'] = np.round(self.slit_pos.separation(coords).to(u.arcmin).value, 2)
        self.data['coords'] = coords

    def fill_row_from_pandas(self, n, series):
        for i, label in enumerate(self.header_labels):
            print(series[self.pandas_labels[i]] is None)
            self.table.setItem(n, i, QTableWidgetItem(str(series[self.pandas_labels[i]])))

    def update_lims(self, lims=None):
        if lims is not None:
            self.ra_lims = lims[0]
            self.dec_lims = lims[1]
        self.req_params['sql'] = f'al2000>{self.ra_lims[0]} and ' + \
                                 f'al2000<{self.ra_lims[1]} and ' + \
                                 f'de2000>{self.dec_lims[0]} and ' + \
                                 f'de2000<{self.dec_lims[1]}'

    def update_slit(self, slit_pos):
        self.slit_pos = slit_pos

    def get_selected_data(self):
        n = self.table.currentRow()
        print(n)
        if n >= 0:
            return self.data.loc[n].to_dict()

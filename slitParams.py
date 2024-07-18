import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from matplotlib import colormaps


class slitParams:
    def __init__(self, n=0, csv_path=None, is_used=True):
        colors = colormaps['tab20'](np.linspace(0, 1, 20))
        self.n = n
        self.colors = colors[[2 * n, 2 * n + 1]]
        self.csv_path = csv_path
        self.is_used = is_used
        try:
            self.dataFrame = pd.read_csv(self.csv_path)
            self.slitpos = SkyCoord(self.dataFrame['RA'],
                                    self.dataFrame['DEC'],
                                    frame='icrs',
                                    unit=(u.hourangle, u.deg))
        except (UnicodeDecodeError, FileNotFoundError, KeyError):
            print('INVALID CSV FILE OR PATH')

    def set_n(self, new_n):
        self.n = new_n
        colors = colormaps['tab20'](np.linspace(0, 1, 20))
        self.colors = colors[[2 * self.n, 2 * self.n + 1]]

    def set_csv_path(self, new_path):
        self.csv_path = new_path
        try:
            self.dataFrame = pd.read_csv(self.csv_path)
        except (UnicodeDecodeError, FileNotFoundError):
            print('INVALID CSV FILE OR PATH')


#! /usr/bin/env python3

# Copyright (C) 2022 The Qt Company Ltd.
# SPDX-License-Identifier: LicenseRef-Qt-Commercial OR BSD-3-Clause

import sys

import numpy as np
import scipy.spatial
import matplotlib
from astropy.visualization import simple_norm
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from astropy.coordinates import SkyCoord
from matplotlib.backend_bases import MouseButton
from matplotlib.legend_handler import HandlerTuple
import pandas as pd

from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from PySide6.QtCore import Slot, Signal
from PySide6.QtWidgets import (
    QApplication,
    QWidget,
    QDoubleSpinBox,
    QHBoxLayout,
    QFormLayout,
    QGridLayout,
    QPushButton,
    QFileDialog,
    QCheckBox,
)

from OpenFile import OpenFile
from radecSpinBox import radecSpinBox
from astroPatches import rotatedEllipse
from InputFiles import InputDialog
from FineMovements import FineMovementsDialog
from slitParams import slitParams
matplotlib.use('QtAgg')


class galParams:
    def __init__(self, i=0., pa=0., vel=0., dist=0., center=None):
        self.i = i
        self.pa = pa
        self.vel = vel
        self.dist = dist
        if center is None:
            self.center = SkyCoord(0, 0, unit=(u.hourangle, u.deg),
                                   frame='icrs')
        else:
            self.center = center
        self.frame = self.center.skyoffset_frame(rotation=self.pa)

    def update_frame(self):
        self.frame = self.center.skyoffset_frame(rotation=self.pa)


class galaxyImage:
    def __init__(self, figure):
        self.figure = figure
        self.wcs: WCS | None = None
        self.image: np.ndarray | None = None
        self.axes_gal: matplotlib.axes.Axes | None = None
        self.norm_im: matplotlib.colors.Normalize | None = None
        self.ellipses = []
        # just default value, replaced with the real galaxy coordinate frame later
        self.gal_frame = SkyCoord(0, 0, unit=(u.deg, u.deg),
                                  frame='icrs').skyoffset_frame()

    def add_image(self, image):
        self.figure.clear()
        self.wcs = WCS(image.header)
        self.axes_gal = self.figure.subplots(
            subplot_kw={'projection': self.wcs})
        self.image = image
        self.norm_im = simple_norm(image.data, 'linear', percent=99.3)
        self.axes_gal.imshow(self.image.data, cmap='bone', norm=self.norm_im)
        # self.plot_galaxy()

    def plot_galaxy(self, gal_p=None):
        if self.axes_gal is None:
            return
        while self.axes_gal.patches:
            list(self.axes_gal.patches)[0].remove()

        if gal_p is not None:
            self.gal_frame = gal_p.frame
            self.plot_ellipses(gal_p.dist, gal_p.i)

    def plot_slit(self, data):
        if self.axes_gal is None:
            return
        [x.remove() for x in self.axes_gal.lines]
        for dat in data:
            mask1, mask2 = dat.dataFrame['mask1'], dat.dataFrame['mask2']
            if len(mask1[mask1]) > 0:
                self.axes_gal.plot(
                    dat.slitpos.ra[mask1],
                    dat.slitpos.dec[mask1],
                    marker='.',
                    linestyle='',
                    transform=self.axes_gal.get_transform('icrs'),
                    color=dat.colors[0])
            if len(mask2[mask2]) > 0:
                self.axes_gal.plot(
                    dat.slitpos.ra[mask2],
                    dat.slitpos.dec[mask2],
                    marker='.',
                    linestyle='',
                    transform=self.axes_gal.get_transform('icrs'),
                    color=dat.colors[1])
        self.figure.canvas.draw()

    def get_image_center(self):
        center = np.array(np.shape(self.image)) * 0.5
        cent_coord = self.wcs.pixel_to_world(*center)
        return cent_coord

    def plot_ellipses(self, d, i):
        r = np.array([5, 10, 15, 20]) * 1e-3
        ang_r = u.radian * r / d
        for a in ang_r:
            self.plot_ellips(a, i)

        # axis
        rmax = u.radian * 20e-3 / d
        self.plot_ellips(rmax * np.cos(np.radians(i)), 90, theta=0 * u.deg, col='tab:olive', lw=0.5)
        self.plot_ellips(rmax, 90, col='tab:olive', lw=0.5)

        # center marker
        self.plot_ellips(u.arcsec * 3, 90, theta=0 * u.deg, col='red', lw=1.5)
        self.plot_ellips(u.arcsec * 3, 90, col='red', lw=1.5)

    def plot_ellips(self, r, i, theta=90 * u.deg, col='tab:olive', lw=0.5):
        el = rotatedEllipse([0 * u.deg, 0 * u.deg],
                            r, r*np.cos(np.radians(i)),
                            theta=theta,
                            edgecolor=col, facecolor='none', lw=lw,
                            transform=self.axes_gal.get_transform(self.gal_frame))
        self.axes_gal.add_patch(el)
        pass


class csvPlot(QWidget):
    del_point = Signal(list)

    def __init__(self, figure):
        super().__init__()
        self.figure = figure
        self.data: list[slitParams] = []
        self.axes_plot: matplotlib.axes.Axes | None = None
        self.figure.canvas.mpl_connect('button_press_event', self.on_click)
        # object to find the closest point in a set
        self.ckdtree = None
        self.all_points = None
        self.all_points_n_line = None
        self.point_chosen = False
        self.point = None
        self.index_chosen = None
        self.last_gal_p = None
        self.scale_x = 1.
        self.scale_y = 1.

    def add_data(self, data):
        self.figure.clear()
        self.data = data
        self.axes_plot = self.figure.subplots()

    def on_click(self, event):
        if event.inaxes != self.axes_plot: return
        if not self.point_chosen and event.button is MouseButton.RIGHT:
            self.index_chosen = self.ckdtree.query([event.xdata / self.scale_x,
                                                    event.ydata / self.scale_y])[1]
            x, y = self.all_points[self.index_chosen]
            self.point = self.axes_plot.plot(x, y, 'ro')[0]
            self.figure.canvas.draw()
            self.point_chosen = True
            return
        elif self.point_chosen:
            if event.button is MouseButton.RIGHT:
                n_line = self.all_points_n_line[self.index_chosen]
                min_i_n_line = np.arange(len(self.all_points_n_line))[(self.all_points_n_line == n_line)].min()
                self.data[n_line].del_element(self.index_chosen - min_i_n_line)
                self.calc_rc(self.last_gal_p)
                self.del_point.emit(self.data)
                self.point_chosen = False
            if event.button is MouseButton.LEFT:
                self.point_chosen = False
                self.point.remove()
                self.figure.canvas.draw()

    def calc_rc(self, gal_p):
        # if it is None, that means there are no data
        if self.axes_plot is None:
            return
        self.last_gal_p = gal_p
        self.axes_plot.clear()
        self.all_points = np.array([]).reshape((0, 2))
        self.all_points_n_line = np.array([], dtype=int)
        for i, dat in enumerate(self.data):
            dat.los_to_rc(gal_p)
            new_points = np.array([dat.dataFrame['R_pc'],
                                   dat.dataFrame['Circular_v']]).T
            new_mask = dat.dataFrame['mask1'] | dat.dataFrame['mask2']
            self.all_points = np.concatenate((self.all_points,
                                              new_points[new_mask]))
            new_index = np.ones(len(new_points[new_mask]), dtype=int) * i
            self.all_points_n_line = np.concatenate((self.all_points_n_line,
                                                     new_index))
        ckd_data = self.all_points.copy()
        self.scale_x = (ckd_data[:, 0].max() - ckd_data[:, 0].min())
        self.scale_y = (ckd_data[:, 1].max() - ckd_data[:, 1].min())
        ckd_data[:, 0] = ckd_data[:, 0] / self.scale_x
        ckd_data[:, 1] = ckd_data[:, 1] / self.scale_y
        self.ckdtree = scipy.spatial.cKDTree(ckd_data)
        self.plot_rc()

    def plot_rc(self):
        self.axes_plot.set_ylabel('Circular Velocity, km/s')
        self.axes_plot.set_xlabel('R, parsec')
        lines = []
        labels = []
        for dat_sp in self.data:
            line = []
            dat = dat_sp.dataFrame
            verr = dat['Circular_v_err'].to_numpy()
            mask1, mask2 = dat['mask1'], dat['mask2']
            if len(mask1[mask1]) > 0:
                line.append(self.axes_plot.errorbar(
                    dat['R_pc'][mask1],
                    dat['Circular_v'][mask1],
                    yerr=verr[mask1],
                    linestyle='-',
                    marker='.',
                    color=dat_sp.colors[0])[0])
            if len(mask2[mask2]) > 0:
                line.append(self.axes_plot.errorbar(
                    dat['R_pc'][mask2],
                    dat['Circular_v'][mask2],
                    yerr=verr[mask2],
                    linestyle='-',
                    marker='.',
                    color=dat_sp.colors[1])[0])
            if line and dat_sp.label:
                lines.append(tuple(line))
                labels.append(dat_sp.label)

        self.axes_plot.legend(lines, labels,
                       handler_map={tuple: HandlerTuple(ndivide=None)})
        self.axes_plot.axhline(y=0, color='black', linestyle='--', lw=0.5)
        self.axes_plot.axvline(x=0, color='black', linestyle='--', lw=0.5)
        self.figure.canvas.draw()

    def return_rc(self):
        res = []
        for n, dat in enumerate(self.data):
            preresult = dat.dataFrame[['RA', 'DEC', 'Circular_v', 'Circular_v_err',
                                       'R_pc', 'R_arcsec', 'mask1', 'mask2']]
            preresult['source_index'] = n
            res.append(preresult)
        res = pd.concat(res, ignore_index=True)
        return res


class PlotWidget(QWidget):
    def __init__(self, parent=None, csv=None, frame=None, refcenter=None, pa=0.,
                 inclination=0., velocity=0.):
        super().__init__(parent)

        self.fine_active = False

        # create widgets
        ################
        # Plots (rotation cuves and galaxy image)
        self.plot_fig = FigureCanvas(Figure(figsize=(5, 3)))
        self.gal_fig = FigureCanvas(Figure(figsize=(5, 3)))
        self.toolbar_plot = NavigationToolbar2QT(self.plot_fig, self)
        self.toolbar_gal = NavigationToolbar2QT(self.gal_fig, self)
        # String input for the path to the galaxy image field
        self.image_field = OpenFile(text='image', mode='o')
        # String input for the pathes to the csvs with LOS velocities
        # self.csv_field = OpenFile(text='csv', mode='n')
        # Inclination input
        self.i_input = QDoubleSpinBox()
        # PA input
        self.PA_input = QDoubleSpinBox()
        # Galaxy center coordinates input
        self.ra_input = radecSpinBox(radec='ra')
        self.dec_input = radecSpinBox(radec='dec')
        # Galaxy center velocity input
        self.vel_input = QDoubleSpinBox()
        # Distance to the galaxy input
        self.dist_input = QDoubleSpinBox()
        self.dist_checkbox = QCheckBox('Calculate from velocity')
        # Buttons
        self.redraw_button = QPushButton(text='Redraw')
        self.saveres_button = QPushButton(text='Save Results')
        self.fine_button = QPushButton(text='Fine movements')
        self.manage_csv_button = QPushButton(text='Manage CSV files')
        # Dialogs
        self.manage_csv = InputDialog(self)
        self.fine_dialog = FineMovementsDialog(self)

        # Configure all widgets
        self.configureElements(frame, csv, inclination, pa, refcenter, velocity)
        self.configureLayout()

        self.galIm = galaxyImage(self.gal_fig.figure)
        self.csvGraph = csvPlot(self.plot_fig.figure)
        self.gal_p = galParams()

        self.redraw_button.clicked.connect(self.redraw)
        self.PA_input.valueChanged.connect(self.galFrameChanged)
        self.ra_input.valueChanged.connect(self.galFrameChanged)
        self.dec_input.valueChanged.connect(self.galFrameChanged)
        self.vel_input.valueChanged.connect(self.galFrameChanged)
        self.i_input.valueChanged.connect(self.galFrameChanged)
        self.dist_input.valueChanged.connect(self.galFrameChanged)
        self.saveres_button.clicked.connect(self.save_rc)
        self.dist_checkbox.stateChanged.connect(self.galFrameChanged)
        self.fine_button.clicked.connect(lambda: self.fine_dialog.show())
        self.manage_csv_button.clicked.connect(lambda: self.manage_csv.show())
        self.fine_dialog.move_ra.connect(self.ra_input.stepByAngle)
        self.fine_dialog.move_dec.connect(self.dec_input.stepByAngle)
        self.csvGraph.del_point.connect(self.galIm.plot_slit)

    def configureElements(self, frame, csv, inclination, pa, refcenter,
                          velocity):
        if frame is not None:
            self.image_field.fill_string(frame)

        if csv is not None:
            for csv_path in csv:
                self.manage_csv.add_file(csv_path)
            # csv = ', '.join(csv)
            # self.csv_field.fill_string(csv)

        self.i_input.setKeyboardTracking(False)
        self.i_input.setMinimum(1)
        self.i_input.setValue(inclination)

        self.PA_input.setKeyboardTracking(False)
        self.PA_input.setMaximum(360.0)
        self.PA_input.setValue(pa)

        self.ra_input.setKeyboardTracking(False)
        self.dec_input.setKeyboardTracking(False)
        if refcenter is not None:
            self.ra_input.setValue(refcenter[0])
            self.dec_input.setValue(refcenter[1])

        self.vel_input.setKeyboardTracking(False)
        self.vel_input.setSuffix('km/s')
        self.vel_input.setMaximum(500000)
        self.vel_input.setValue(velocity)

        self.dist_input.setKeyboardTracking(False)
        self.dist_input.setSuffix('Mpc')
        self.dist_input.setMinimum(0.01)
        self.dist_input.setValue(velocity / 70)
        self.dist_input.setSingleStep(0.1)
        self.dist_input.setDisabled(True)

        self.dist_checkbox.setToolTip('Assuming H0=70km/s/Mpc')
        self.dist_checkbox.setChecked(True)

    def configureLayout(self):
        # Layout
        r_button_layout = QHBoxLayout()
        r_button_layout.addWidget(self.redraw_button)
        r_button_layout.addWidget(self.saveres_button)

        l_button_layout = QHBoxLayout()
        # l_button_layout.addWidget(self.manage_csv_button)
        l_button_layout.addWidget(self.fine_button)

        left_layout = QFormLayout()
        # left_layout.addRow(self.csv_field)
        left_layout.addRow(self.manage_csv_button)
        left_layout.addRow('i', self.i_input)
        left_layout.addRow('RA', self.ra_input)
        left_layout.addRow('system velocity', self.vel_input)
        left_layout.addRow(l_button_layout)
        right_layout = QFormLayout()
        right_layout.addRow(self.image_field)
        right_layout.addRow('PA', self.PA_input)
        right_layout.addRow('DEC', self.dec_input)
        right_layout.addRow(self.dist_checkbox, self.dist_input)
        right_layout.addRow(r_button_layout)

        glayout = QGridLayout()
        glayout.addWidget(self.toolbar_plot, 0, 0)
        glayout.addWidget(self.toolbar_gal, 0, 1)
        glayout.addWidget(self.plot_fig, 1, 0)
        glayout.addWidget(self.gal_fig, 1, 1)
        glayout.addLayout(left_layout, 2, 0)
        glayout.addLayout(right_layout, 2, 1)
        glayout.setRowStretch(0, 0)
        glayout.setRowStretch(1, 1)
        glayout.setRowStretch(2, 0)
        self.setLayout(glayout)

    @Slot()
    def calc_dist(self):
        if self.dist_checkbox.isChecked():
            self.dist_input.setValue(self.vel_input.value() / 70.)
            self.dist_input.setDisabled(True)
        else:
            self.dist_input.setDisabled(False)

    @Slot()
    def galFrameChanged(self):
        self.updateValues()
        self.galIm.plot_galaxy(self.gal_p)
        self.csvGraph.calc_rc(self.gal_p)
        self.galIm.plot_slit(self.csvGraph.data)

    @Slot()
    def redraw(self):
        """ Update the plot with the current input values """
        self.updateValues()

        image = fits.open(self.image_field.files)[0]
        self.galIm.add_image(image)
        im_center = self.galIm.get_image_center()
        if im_center.separation(self.gal_p.center) > 1 * u.deg:
            self.ra_input.setValue(im_center.ra)
            self.dec_input.setValue(im_center.dec)
        self.galIm.plot_galaxy(self.gal_p)

        self.csvGraph.add_data(self.manage_csv.data)
        self.csvGraph.calc_rc(self.gal_p)
        self.galIm.plot_slit(self.csvGraph.data)

    @Slot()
    def save_rc(self):
        # filenames = self.csv_field.return_filenames()
        fname = self.manage_csv.data[0].csv_path
        result = self.csvGraph.return_rc()
        regexps = "CSV (*.csv)"
        fname_temp = '/'.join(fname.split('/')[:-1]) + '/rot_curve.csv'
        file_path = QFileDialog.getSaveFileName(self,
                                                "Save rotation curve",
                                                fname_temp,
                                                regexps)[0]
        result.to_csv(file_path)
        print('Saving ', file_path)

    def updateValues(self):
        self.gal_p.i = self.i_input.value() * u.deg
        self.gal_p.pa = self.PA_input.value() * u.deg
        self.gal_p.center = SkyCoord(self.ra_input.getAngle(),
                                     self.dec_input.getAngle(),
                                     frame='icrs')
        self.gal_p.vel = self.vel_input.value()
        self.calc_dist()
        self.gal_p.dist = self.dist_input.value()
        self.gal_p.update_frame()

        self.fine_dialog.center = SkyCoord(self.ra_input.getAngle(),
                                           self.dec_input.getAngle(),
                                           frame='icrs')
        self.fine_dialog.pa = self.PA_input.value()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--csv', nargs='+', default=None,
                        help='''csv or similar file with positions and
                        velocities''')
    parser.add_argument('-r', '--refcenter', nargs=2, default=None,
                        help='''coordinates of center of galaxy''')
    parser.add_argument('-v', '--velocity', type=float, default=0.0,
                        help='system velocity')
    parser.add_argument('-p', '--PA', type=float, default=0.0,
                        help='galaxy PA')
    parser.add_argument('-i', '--inclination', type=float, default=0.0,
                        help='inclination of galaxy')
    parser.add_argument('-f', '--frame', default=None,
                        help='frame with image')
    pargs = parser.parse_args(sys.argv[1:])

    app = QApplication(sys.argv)
    w = PlotWidget(None, pargs.csv, pargs.frame, pargs.refcenter, pargs.PA,
                   pargs.inclination, pargs.velocity)
    w.show()
    sys.exit(app.exec())

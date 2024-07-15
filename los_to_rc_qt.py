#! /usr/bin/env python3

# Copyright (C) 2022 The Qt Company Ltd.
# SPDX-License-Identifier: LicenseRef-Qt-Commercial OR BSD-3-Clause

import sys

import numpy as np
# from scipy.stats import norm
import matplotlib
from astropy.visualization import simple_norm
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from astropy.coordinates import SkyCoord, ICRS
import pandas as pd
# from itertools import zip_longest, chain
from matplotlib import colormaps

# from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from PySide6.QtCore import Slot, Qt
from PySide6.QtGui import QShortcut, QKeySequence
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
from InputFiles import InputDialog, slitParams
matplotlib.use('QtAgg')


def plot_slit_points(ax, rel_slit, masks=None, gal_frame=None):
    if masks is None:
        masks = [np.ones(len(rel_slit)).astype(bool)]
    for mask in masks:
        if len(mask[mask]) > 0:
            ax.plot(rel_slit.ra[mask], rel_slit.dec[mask], marker='.',
                    linestyle='', transform=ax.get_transform(gal_frame))


def los_to_rc(data, slit, gal_frame, inclination, sys_vel, dist,
              verr_lim=200):
    """
    data - pd.DataFrame
    slit - SkyCoord (inerable - SkyCoord of list)
    gal_frame - coordinates frame
    inclination - float * u.deg
    sys_vel - float
    obj_name - str
    verr_lim - float
    """
    # H0 = 70 / (1e+6 * u.parsec)
    if 'position' in data:
        slit_pos = data['position']
    else:
        slit_pos = slit.separation(slit[0]).to(u.arcsec)

    gal_center = SkyCoord(0 * u.deg, 0 * u.deg, frame=gal_frame)
    rel_slit = slit.transform_to(gal_frame)

    dist = dist * 1000000 * u.parsec
    # dist = sys_vel / H0
    # Исправляем за наклон галактики
    rel_slit_corr = SkyCoord(rel_slit.lon / np.cos(inclination), rel_slit.lat,
                             frame=gal_frame)
    # Угловое расстояние точек щели до центра галактики
    # (с поправкой на наклон галактики)
    separation = rel_slit_corr.separation(gal_center)
    # Физическое расстояние
    r_slit = dist * np.sin(separation)

    # Угол направления на центр галактики
    gal_frame_center = gal_center.transform_to(gal_frame)
    slit_gal_pa = gal_frame_center.position_angle(rel_slit_corr)

    vel_lon = (data['velocity'].to_numpy() - sys_vel) / np.sin(inclination)
    if 'velocity_err' in data:
        vel_lon_err = np.abs(data['velocity_err'].to_numpy() / np.sin(inclination))
    else:
        vel_lon_err = data['velocity'].to_numpy() * 0

    vel_r = vel_lon / np.cos(slit_gal_pa)
    vel_r_err = np.abs(vel_lon_err / np.cos(slit_gal_pa))

    mask = (vel_r_err < verr_lim)
    # mask_center = (separation.to(u.arcsec) > 5 * u.arcsec)
    # cone_angle = 20 * u.deg
    # mask_cone_1 = (slit_gal_pa > 90 * u.deg - cone_angle) & \
    #               (slit_gal_pa < 90 * u.deg + cone_angle)
    # mask_cone_2 = (slit_gal_pa > 270 * u.deg - cone_angle) & \
    #               (slit_gal_pa < 270 * u.deg + cone_angle)
    # mask_cone = ~(mask_cone_1 | mask_cone_2)
    # mask = mask & mask_center
    # mask = mask & mask_cone

    # lat = np.array(rel_slit_corr.lat.to(u.arcsec)/u.arcsec)
    # minor_ax = np.argmin(np.abs(lat))

    closest_point = np.argmin(np.abs(r_slit))

    first_side = (slit_pos >= slit_pos[closest_point])
    second_side = (slit_pos < slit_pos[closest_point])
    first_side_mask = (first_side & mask)
    second_side_mask = (second_side & mask)

    data['Circular_v'] = -vel_r
    data['Circular_v_err'] = vel_r_err
    data['R_pc'] = r_slit / u.parsec
    data['R_arcsec'] = separation.to(u.arcsec) / u.arcsec
    data['mask1'] = np.array(first_side_mask, dtype=bool)
    data['mask2'] = np.array(second_side_mask, dtype=bool)

    return data


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
    def __init__(self, figure, image):
        self.colors = colormaps['tab20'](np.linspace(0, 1, 20))
        self.wcs = WCS(image.header)
        self.figure = figure
        self.axes_gal = figure.subplots(
            subplot_kw={'projection': self.wcs})
        self.image = image
        self.norm_im = simple_norm(image.data, 'linear', percent=99.3)
        self.slits = None
        self.masks = None
        self.slit_draws = None
        self.ellipses = []
        # just default value, replaced with the real galaxy coordinate frame later
        self.gal_frame = SkyCoord(0, 0, unit=(u.deg, u.deg),
                                  frame='icrs').skyoffset_frame()
        self.overlay = self.axes_gal.get_coords_overlay(self.gal_frame)
        self.plot_galaxy()

    def plot_galaxy(self, gal_p=None):
        self.axes_gal.clear()
        self.axes_gal = self.figure.subplots(
            subplot_kw={'projection': self.wcs})
        self.axes_gal.imshow(self.image.data, cmap='bone', norm=self.norm_im)
        # "Стираем" чёрточки по краям картинки
        self.axes_gal.coords['ra'].set_ticks(color='white')
        self.axes_gal.coords['dec'].set_ticks(color='white')

        if gal_p is not None:
            self.gal_frame = gal_p.frame
            self.overlay = self.axes_gal.get_coords_overlay(self.gal_frame)
            # "Стираем" чёрточки по краям картинки
            self.overlay['lon'].set_ticks(color='white')
            self.overlay['lat'].set_ticks(color='white')
            self.overlay['lon'].set_ticklabel(alpha=0)
            self.overlay['lat'].set_ticklabel(alpha=0)
            self.overlay.grid(color='white', linestyle='solid', alpha=0.5)
            # print set center of the galaxy
            self.axes_gal.plot(0, 0, 'r+', ms=10,
                               transform=self.axes_gal.get_transform(self.gal_frame))
            self.plot_ellipses(gal_p.dist, gal_p.i)

        if self.slits is not None:
            self.plot_slit(self.slits, self.masks)

        # plot_galaxy(self.axes_gal, self.image, self.gal_frame)

    def plot_slit(self, slits, masks):
        self.slits = slits
        self.masks = masks
        for slit, mask in zip(slits, masks):
            plot_slit_points(self.axes_gal, slit, mask,
                             'icrs')

        for line in self.axes_gal.lines:
            list(self.axes_gal.lines).remove(line)

        for slit, mask, i in zip(slits, masks, range(0, 20, 2)):
            mask1, mask2 = mask
            if len(mask1[mask1]) > 0:
                self.axes_gal.plot(
                    slit.ra[mask1],
                    slit.dec[mask1],
                    marker='.',
                    linestyle='',
                    transform=self.axes_gal.get_transform('icrs'),
                    color=self.colors[i])
            if len(mask2[mask2]) > 0:
                self.axes_gal.plot(
                    slit.ra[mask2],
                    slit.dec[mask2],
                    marker='.',
                    linestyle='',
                    transform=self.axes_gal.get_transform('icrs'),
                    color=self.colors[i + 1])

    def get_center(self):
        center = np.array(np.shape(self.image)) * 0.5
        print(center)
        cent_coord = self.wcs.pixel_to_world(*center)
        return cent_coord

    def plot_ellipses(self, d, i):
        r = np.array([5, 10, 15, 20]) * 1e-3
        ang_r = u.radian * r / d
        for a in ang_r:
            self.plot_ellips(a, i)

    def plot_ellips(self, r, i):
        el = rotatedEllipse([0 * u.deg, 0 * u.deg],
                            r, r*np.cos(np.radians(i)),
                            theta=90 * u.deg,
                            edgecolor='tab:olive', facecolor='none', lw=0.5,
                            transform=self.axes_gal.get_transform(self.gal_frame))
        self.axes_gal.add_patch(el)
        pass


class csvPlot:
    def __init__(self, data, figure):
        self.colors = colormaps['tab20'](np.linspace(0, 1, 20))
        # data - list of pd.DataFrame
        self.data = data
        self.slits = []
        self.masks = []
        for dat in self.data:
            slit_ra = dat['RA']
            slit_dec = dat['DEC']
            self.slits.append(SkyCoord(slit_ra, slit_dec, frame='icrs',
                                       unit=(u.hourangle, u.deg)))
        self.axes_plot = figure.subplots()

    # def calc_rc(self, gal_frame, inclination, sys_vel, dist=None):
    def calc_rc(self, gal_p):
        # if dist is None:
        #     dist = sys_vel / 70.
        self.axes_plot.clear()
        self.masks = []
        for dat, slit in zip(self.data, self.slits):
            dat = los_to_rc(dat, slit, gal_p.frame, gal_p.i, gal_p.vel, gal_p.dist)
            self.masks.append([dat['mask1'].to_numpy(),
                               dat['mask2'].to_numpy()])
        self.plot_rc()
        return self.slits, self.masks

    def plot_rc(self):
        self.axes_plot.set_ylabel('Circular Velocity, km/s')
        self.axes_plot.set_xlabel('R, parsec')
        for dat, mask, i in zip(self.data, self.masks, range(0, 20, 2)):
            verr = dat['Circular_v_err'].to_numpy()
            mask1, mask2 = mask
            if len(mask1[mask1]) > 0:
                self.axes_plot.errorbar(
                    dat['R_pc'][mask1],
                    dat['Circular_v'][mask1],
                    yerr=verr[mask1],
                    linestyle='',
                    marker='.',
                    color=self.colors[i])
            if len(mask2[mask2]) > 0:
                self.axes_plot.errorbar(
                    dat['R_pc'][mask2],
                    dat['Circular_v'][mask2],
                    yerr=verr[mask2],
                    linestyle='',
                    marker='.',
                    color=self.colors[i + 1])

    def return_rc(self):
        return self.data


class PlotWidget(QWidget):
    def __init__(self, parent=None, csv=None, frame=None, refcenter=None, pa=0.,
                 inclination=0., velocity=0.):
        super().__init__(parent)

        self.gal_changed = False
        self.csv_changed = False
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
        # Keyboard Shortcuts (for fine movements)
        self.left_shortcut = QShortcut(QKeySequence('Left'), self)
        self.right_shortcut = QShortcut(QKeySequence('Right'), self)
        self.up_shortcut = QShortcut(QKeySequence('Up'), self)
        self.down_shortcut = QShortcut(QKeySequence('Down'), self)
        self.fine_shortcuts = [self.left_shortcut,
                               self.right_shortcut,
                               self.up_shortcut,
                               self.down_shortcut]
        # Dialogs
        self.manage_csv = InputDialog(self)

        # Configure all widgets
        self.configureElements(frame, csv, inclination, pa, refcenter, velocity)
        self.configureLayout()

        self.galIm = None
        self.csvGraph = None
        self.gal_p = galParams()

        self.redraw_button.clicked.connect(self.redraw)
        # self.csv_field.changed_path.connect(self.csvChanged)
        self.image_field.changed_path.connect(self.galChanged)
        self.PA_input.valueChanged.connect(self.galFrameChanged)
        self.ra_input.valueChanged.connect(self.galFrameChanged)
        self.dec_input.valueChanged.connect(self.galFrameChanged)
        self.vel_input.valueChanged.connect(self.galFrameChanged)
        self.i_input.valueChanged.connect(self.galFrameChanged)
        self.dist_input.valueChanged.connect(self.galFrameChanged)
        self.saveres_button.clicked.connect(self.save_rc)
        self.dist_checkbox.stateChanged.connect(self.galFrameChanged)
        self.fine_button.clicked.connect(self.fineMovements)
        self.left_shortcut.activated.connect(lambda: self.ra_input.stepBy(0.1))
        self.right_shortcut.activated.connect(lambda: self.ra_input.stepBy(-0.1))
        self.up_shortcut.activated.connect(lambda: self.dec_input.stepBy(0.1))
        self.down_shortcut.activated.connect(lambda: self.dec_input.stepBy(-0.1))
        self.manage_csv_button.clicked.connect(lambda: self.manage_csv.show())
        self.manage_csv.accepted.connect(self.csvChanged)

    def configureElements(self, frame, csv, inclination, pa, refcenter,
                          velocity):
        if frame is not None:
            self.image_field.fill_string(frame)
            self.gal_changed = True

        if csv is not None:
            for csv_path in csv:
                self.manage_csv.add_file(csv_path)
            # csv = ', '.join(csv)
            # self.csv_field.fill_string(csv)
            self.csv_changed = True

        self.i_input.setKeyboardTracking(False)
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
        self.dist_input.setValue(velocity / 70)
        self.dist_input.setSingleStep(0.1)
        self.dist_input.setDisabled(False)

        self.dist_checkbox.setToolTip('Assuming H0=70km/s/Mpc')

        self.fine_button.setCheckable(True)
        self.fineMovements()

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
    def galChanged(self):
        self.gal_changed = True

    @Slot()
    def csvChanged(self):
        self.csv_changed = True

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
        slits, masks = self.csvGraph.calc_rc(self.gal_p)
        self.galIm.plot_slit(slits, masks)
        self.gal_fig.draw()
        self.plot_fig.draw()

    @Slot()
    def kinematicsChanged(self):
        self.updateValues()
        self.csvGraph.calc_rc(self.gal_p)
        self.plot_fig.draw()

    @Slot()
    def redraw(self):
        """ Update the plot with the current input values """
        self.updateValues()

        if self.gal_changed:
            self.gal_fig.figure.clear()
            image = fits.open(self.image_field.files)[0]
            self.galIm = galaxyImage(self.gal_fig.figure, image)
            im_center = self.galIm.get_center()
            if im_center.separation(self.gal_p.center) > 1 * u.deg:
                self.ra_input.setValue(im_center.ra)
                self.dec_input.setValue(im_center.dec)
            self.galIm.plot_galaxy(self.gal_p)
            self.gal_changed = False

        if self.csv_changed:
            self.plot_fig.figure.clear()
            data = [x.dataFrame for x in self.manage_csv.data]
            self.csvGraph = csvPlot(data, self.plot_fig.figure)
            slits, masks = self.csvGraph.calc_rc(self.gal_p)
            if self.galIm is not None:
                self.galIm.plot_galaxy(self.gal_p)
                self.galIm.plot_slit(slits, masks)
            self.csv_changed = False

        self.gal_fig.draw()
        self.plot_fig.draw()

    @Slot()
    def save_rc(self):
        # filenames = self.csv_field.return_filenames()
        filenames = [x.csv_path for x in self.manage_csv.data]
        dataframes = self.csvGraph.return_rc()
        regexps = "CSV (*.csv)"
        for fname, dat in zip(filenames, dataframes):
            fname_temp = '.'.join(fname.split('.')[:-1]) + '_rc.csv'
            file_path = QFileDialog.getSaveFileName(self,
                                                    "Save rotation curve",
                                                    fname_temp,
                                                    regexps)[0]
            print('Saving ', file_path)
            dat[['RA', 'DEC', 'Circular_v', 'Circular_v_err',
                 'R_pc', 'R_arcsec', 'mask1', 'mask2']].to_csv(file_path)

    def updateValues(self):
        print('Fine button is pressed: ', self.fine_button.isChecked())
        self.gal_p.i = self.i_input.value() * u.deg
        self.gal_p.pa = self.PA_input.value() * u.deg
        self.gal_p.center = SkyCoord(self.ra_input.getAngle(),
                                     self.dec_input.getAngle(),
                                     frame='icrs')
        self.gal_p.vel = self.vel_input.value()
        self.calc_dist()
        self.gal_p.dist = self.dist_input.value()
        self.gal_p.update_frame()
        self.calc_dist()

    def fineMovements(self):
        # QApplication.setOverrideCursor(Qt.CrossCursor)
        self.fine_active = self.fine_button.isChecked()
        print(self.fine_active)
        if self.fine_active:
            for s in self.fine_shortcuts:
                s.blockSignals(False)
        else:
            for s in self.fine_shortcuts:
                s.blockSignals(True)


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

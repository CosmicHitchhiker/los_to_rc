import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from matplotlib import colormaps


class slitParams:
    palette = colormaps['tab20'](np.linspace(0, 1, 20))
    avaliable_n = np.arange(20).tolist()

    def __init__(self, csv_path=None, is_used=True):
        m = slitParams.avaliable_n[0]
        slitParams.avaliable_n.pop(0)
        self.colors = self.palette[[2 * m, 2 * m + 1]]
        self.csv_path = csv_path
        self.is_used = is_used
        try:
            self.dataFrame = pd.read_csv(self.csv_path)
            self.dataFrame.dropna(subset=['velocity', 'velocity_err'],
                                  inplace=True)
            self.slitpos = SkyCoord(self.dataFrame['RA'],
                                    self.dataFrame['DEC'],
                                    frame='icrs',
                                    unit=(u.hourangle, u.deg))
        except (UnicodeDecodeError, FileNotFoundError, KeyError):
            print('INVALID CSV FILE OR PATH')

    def set_csv_path(self, new_path):
        self.csv_path = new_path
        try:
            self.dataFrame = pd.read_csv(self.csv_path)
        except (UnicodeDecodeError, FileNotFoundError):
            print('INVALID CSV FILE OR PATH')

    def los_to_rc(self, gal_p, verr_lim=200.):
        # los_to_rc(data, slit, gal_frame, inclination, sys_vel, dist,
        #           verr_lim=200):
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
        if 'position' in self.dataFrame:
            slit_pos = self.dataFrame['position'].to_numpy()
        else:
            slit_pos = self.slitpos.separation(self.slitpos[0]).to(u.arcsec)

        gal_center = SkyCoord(0 * u.deg, 0 * u.deg, frame=gal_p.frame)
        rel_slit = self.slitpos.transform_to(gal_p.frame)

        dist = gal_p.dist * 1e+6 * u.parsec
        # dist = sys_vel / H0
        # Исправляем за наклон галактики
        rel_slit_corr = SkyCoord(rel_slit.lon / np.cos(gal_p.i), rel_slit.lat,
                                 frame=gal_p.frame)
        # Угловое расстояние точек щели до центра галактики
        # (с поправкой на наклон галактики)
        separation = rel_slit_corr.separation(gal_center)
        # Физическое расстояние
        r_slit = dist * np.sin(separation)

        # Угол направления на центр галактики
        gal_frame_center = gal_center.transform_to(gal_p.frame)
        slit_gal_pa = gal_frame_center.position_angle(rel_slit_corr)

        vel_lon = (self.dataFrame['velocity'].to_numpy() - gal_p.vel) / np.sin(gal_p.i)
        if 'velocity_err' in self.dataFrame:
            vel_lon_err = np.abs(self.dataFrame['velocity_err'].to_numpy() / np.sin(gal_p.i))
        else:
            vel_lon_err = self.dataFrame['velocity'].to_numpy() * 0

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

        # the closest point to the center of the galaxy
        closest_point = np.argmin(np.abs(r_slit))

        first_side = (slit_pos >= slit_pos[closest_point])
        second_side = (slit_pos < slit_pos[closest_point])
        first_side_mask = (first_side & mask)
        second_side_mask = (second_side & mask)

        self.dataFrame['Circular_v'] = -vel_r
        self.dataFrame['Circular_v_err'] = vel_r_err
        self.dataFrame['R_pc'] = r_slit / u.parsec
        self.dataFrame['R_arcsec'] = separation.to(u.arcsec) / u.arcsec
        self.dataFrame['mask1'] = np.array(first_side_mask, dtype=bool)
        self.dataFrame['mask2'] = np.array(second_side_mask, dtype=bool)

        return self.dataFrame

    def del_element(self, index):
        mask = (self.dataFrame['mask1'] | self.dataFrame['mask2'])
        i = self.dataFrame[mask].iloc[[index]].index.values[0]
        self.dataFrame.drop(i, inplace=True)
        self.slitpos = SkyCoord(self.dataFrame['RA'],
                                self.dataFrame['DEC'],
                                frame='icrs',
                                unit=(u.hourangle, u.deg))

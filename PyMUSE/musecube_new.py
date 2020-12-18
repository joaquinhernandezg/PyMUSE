import copy
import os

import numpy as np
import astropy.units as u
import astropy

from astropy.io import fits
from mpdaf.obj import WCS, WaveCoord

from .image import Image
from .base import Base

__all__ = ('MuseCube')


def remove_dims_from_header(header, dim="3"):
    # no es buena idea hacer esto si se van a hacer espectros
    header = copy.deepcopy(header)
    l = []

    for i in header.keys():
        if dim in i:
            l.append(i)
    for i in l:
        del header[i]

    return header


class MuseCube(Base):
    """
    Class to handle VLT/MUSE data

    """
    n_dims = 3

    def __init__(self, filename_cube=None, filename_white=None, pixelsize=0.2 * u.arcsec,
                 flux_units=1E-20 * u.erg / u.s / u.cm ** 2 / u.angstrom, input_wave_cal='air', data=None, stat=None,
                 header_0=None, header_1=None):
        """
        Parameters
        ----------
        filename_cube: string
            Name of the MUSE datacube .fits file
        filename_white: string
            Name of the MUSE white image .fits file
        pixel_size : float or Quantity, optional
            Pixel size of the datacube, if float it assumes arcsecs.
            Default is 0.2 arcsec
        n_fig : int, optional
            Figure to display the canvas
        flux_units : Quantity
        vmin, vmax: Display paramters
        input_wave_cal: Wavelength calibration of the data cube, can be either air or vac
        output_wave_cal: Wavelength calibration of the spectra extractied with this package

        """
        super(MuseCube, self).__init__(filename_cube, data, stat, header_0, header_1)

        if input_wave_cal not in ['air', 'vac']:
            raise Warning('input_wave_cal and output_wave_cal should be either "air" or "vac"')
        # agregar mas verificacion de parametros

        self.wave_cal = input_wave_cal
        self.pixel_size = pixelsize
        self.flux_units = flux_units

        self.filename_white = filename_white

        print("MuseCube: Ready!")

    def _Base__set_coordinates_from_header(self):
        """
        Sets de WCS for x, y spatial coordinates and another for the z spectral coordinate
        :return: None
        """
        self.wcs = WCS(hdr=self.header_1)
        self.wave = WaveCoord(hdr=self.header_1)


    def create_white(self, new_white_fitsname='white_from_colapse.fits', stat=False, save=True):
        """
        Function that collapses all wavelengths available to produce a new white image
        :param new_white_fitsname: Name of the new image
        :param stat. Default = False. If True, the new white image will be done using the stat extension of the cube.
        :param Save. Default = True. If False, the new returned image will not be saved in the hard disk
        :return: white_image.
        """
        white = self.sum(stat=stat)
        if new_white_fitsname:
            white.write_to_fits(new_white_fitsname)
        return white

    def __iter__(self):
        for i in range(self.flux.shape[0]):
            yield self[i, :, :]

    def __getitem__(self, item):

        """
        modify for including wavelenght ranges in astropy units
        modify header when cutting
        """
        """ 
        if isinstance(item, tuple):
            print(item[0].start, item[0].stop)
            if isinstance(item[0].start, astropy.units.Quantity) and isinstance(item[0].stop, astropy.units.Quantity):
                print(item[0].start, item[0].stop)"""  # implement unit detection

        flux = self.flux.__getitem__(item)
        stat = self.stat.__getitem__(item)

        if isinstance(flux, np.ndarray):
            if flux.ndim == 3:
                # should modify headers
                return self._to_obj(MuseCube, flux=flux, stat=stat)
            elif flux.ndim == 2:
                # should modify headers
                return self._to_obj(Image, flux=flux, stat=stat)
            elif flux.ndim == 1:
                return Spectrum
        else:
            return flux

    def copy(self):
        """
        Makes full copy of the cube
        :return: MuseCube
        """
        return MuseCube(flux_units=self.flux_units, pixelsize=self.pixel_size, data=self.flux.copy(),
                        stat=self.stat.copy(),
                        header_1=self.header_1.copy(), header_0=self.header_0.copy(),
                        input_wave_cal=self.wave_cal)  # modify copy

    def sum(self, stat=False):
        """
        Collapse the cube over the z axis
        :param stat: True if you want to conserve the statistics in the summed image, False if not
        :return: Image
        """
        flux = self.flux.sum(axis=0)
        if stat:
            stat = self.stat.sum(axis=0)
        else:
            stat = None
        return Image(flux_units=self.flux_units, pixelsize=self.pixel_size, data=flux, stat=stat,
                     header_1=self.header_1.copy(), header_0=self.header_0.copy())  # modify copy

    def median(self, stat=None):
        """
        Calculates the median over the z axis
        :param stat: True if you want to conserve the statistics in the summed image, False if not
        :return: Image
        """
        flux = self.flux.median(axis=0)
        if stat:
            stat = self.stat.median(axis=0)
        else:
            stat = None
        return Image(flux_units=self.flux_units, pixelsize=self.pixel_size, data=flux, stat=stat,
                     header_1=self.header_1.copy(), header_0=self.header_0.copy())  # modify copy

    @property
    def wavelength(self):
        """
        Creates the wavelength array for the spectrum. The values of dw, and limits will depend
        of the data and should be revised.
        :return: w: array[]
                 array which contain an evenly sampled wavelength range
        """
        return self.wave.coord()

    def _to_obj(self, obj, flux, stat=None):
        """
        Given flux and stat. To create another object conserving the information of the cube. It is intended for
        private use only, but it can be used.
        :param obj:
        :param flux:
        :param stat:
        :return:
        """

        header_0 = self.header_0
        header_1 = self.header_1

        if obj == Image:
            # change header
            return Image(data=flux, stat=stat, header_0=header_0, header_1=header_1)

        elif obj == MuseCube:

            return MuseCube(data=flux, stat=stat, header_0=header_0, header_1=header_1)

        elif obj == Spectrum:
            return Spectrum()

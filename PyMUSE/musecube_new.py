import copy
import os

import numpy as np
import astropy.units as u
import astropy

from astropy.io import fits
from mpdaf.obj import WCS, WaveCoord

from .image import Image

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


class MuseCube:
    """
    Class to handle VLT/MUSE data

    """

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

        # init
        if input_wave_cal not in ['air', 'vac']:
            raise Warning('input_wave_cal and output_wave_cal should be either "air" or "vac"')
        # agregar mas verificacion de parametros

        self.wave_cal = input_wave_cal
        self.flux_units = flux_units
        self.pixel_size = pixelsize

        self.header_0 = header_0
        self.header_1 = header_1
        self.filename_white = filename_white

        self._flux = None
        self._stat = None

        if filename_cube:
            self.__load_from_fits_file(filename_cube)
        else:
            self.flux = data
            self.stat = stat

        if self.header_1:
            self.wcs = WCS(hdr=self.header_1)
            self.wave = WaveCoord(hdr=self.header_1)
            # read pixelsize
            # read flux_units BUNIT

        print("MuseCube: Ready!")

    @property
    def flux(self):
        return self._flux

    @flux.setter
    def flux(self, value):
        """
        :param value:
        :return Mone:
        """
        # modify when value is None
        if not isinstance(value, np.ndarray) or value.ndim != 3:
            raise ValueError(f"Invalid flux dimensions, got {value.ndim}, expected 3")
        if self.flux is None:
            self._flux = value
        elif value.shape != self.stat.shape:
            raise ValueError(f"Stat and flux can not have different dimensions, try creating a copy instead")
        else:
            self._flux = value

    @property
    def stat(self):
        return self._stat

    @stat.setter
    def stat(self, value):
        """

        :param value:
        :return:
        """
        #mofify when value is None
        if not isinstance(value, np.ndarray) or value.ndim != 3:
            raise ValueError(f"Invalid flux dimensions, got {value.ndim}, expected 3")
        if self.stat is None:
            self._stat = value
        elif value.shape != self.flux.shape:
            raise ValueError(f"Stat and flux can not have different dimensions, try creating a cube instead")
        else:
            self._flux = value

    def __load_from_fits_file(self, filename):
        """

        :param filename: string
        :return:
        """
        hdulist = fits.open(filename)
        self.header_0 = hdulist[0].header
        self.header_1 = hdulist[1].header

        self.flux = hdulist[1].data.astype(np.float64)
        self.stat = hdulist[2].data.astype(np.float64)
        hdulist.close()
        return

    def create_white(self, new_white_fitsname='white_from_colapse.fits', stat=False, save=True):
        """
        Function that collapses all wavelengths available to produce a new white image
        :param new_white_fitsname: Name of the new image
        :param stat. Default = False. If True, the new white image will be done using the stat extension of the cube.
        :param Save. Default = True. If False, the new returned image will not be saved in the hard disk
        :return: white_image.
        """
        white = self.sum(stat=stat)
        return white

    def __getitem__(self, item):
        """
        modify for including wavelenght ranges in astropy units
        """
        """ 
        if isinstance(item, tuple):
            print(item[0].start, item[0].stop)
            if isinstance(item[0].start, astropy.units.Quantity) and isinstance(item[0].stop, astropy.units.Quantity):
                print(item[0].start, item[0].stop)""" #implement unit detection

        flux = self.flux.__getitem__(item)
        stat = self.stat.__getitem__(item)

        if isinstance(flux, np.ndarray):
            if flux.ndim == 3:
                # should modify headers
                return self.__to_obj(MuseCube, flux=flux, stat=stat)
            elif flux.ndim == 2:
                # should modify headers
                return self.__to_obj(Image, flux=flux, stat=stat)
            elif flux.ndim == 1:
                return Spectrum
        else:
            return flux

    def __iter__(self):
        for i in range(self.flux.shape[0]):
            yield self[i, :, :]

    def __to_obj(self, obj, flux, stat=None):
        pixelsize = self.pixel_size
        flux_units = self.flux_units
        header_0 = self.header_0
        header_1 = self.header_1

        if obj == Image:
            header_0 = remove_dims_from_header(header_0)
            header_1 = remove_dims_from_header(header_1)
            return Image(flux_units=flux_units, pixelsize=pixelsize,
                         data=flux, stat=stat, header_0=header_0, header_1=header_1)

        elif obj == MuseCube:
            # modify header for wcs
            return MuseCube(flux_units=flux_units, pixelsize=pixelsize,
                            data=flux, stat=stat, header_0=header_0, header_1=header_1)

        elif obj == Spectrum:
            return Spectrum()

    def copy(self):
        return MuseCube(flux_units=self.flux_units, pixelsize=self.pixel_size, data=self.flux.copy(),
                        stat=self.stat.copy(),
                        header_1=self.header_1.copy(), header_0=self.header_0.copy(),
                        input_wave_cal=self.wave_cal)  # modify copy

    def sum(self, stat=False):
        flux = self.flux.sum(axis=0)
        if stat:
            stat = self.stat.sum(axis=0)
        else:
            stat = None
        return Image(flux_units=self.flux_units, pixelsize=self.pixel_size, data=flux, stat=stat,
                     header_1=self.header_1.copy(), header_0=self.header_0.copy())  # modify copy

    def median(self, stat=None):
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

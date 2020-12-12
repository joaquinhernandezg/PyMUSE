import astropy.units as u
from astropy.io import fits
import numpy as np
from mpdaf.obj import WCS
import copy


__all__ = ('Image')

class Image:
    def __init__(self, filename=None, pixelsize=0.2*u.arcsec,  flux_units=1E-20 * u.erg / u.s / u.cm ** 2 / u.angstrom, data=None, stat=None,
                 header_0=None, header_1=None):

        # agregar mas verificacion de parametros
        # verificar archivo existe
        # verificar headers son astropy headers
        # verificar pixel size and flux units are apropiated astropy units
        # verificar data y stat son np array
        self.flux_units = flux_units #shoud be read_from_file
        self.pixel_size = pixelsize #should be read from file

        self.header_0 = header_0
        self.header_1 = header_1

        if filename:
            self.__load_from_fits_file(filename)
        else:
            self.flux = data
            self.stat = stat

        if self.header_1:
            self.wcs = WCS(hdr=self.header_1)
            #read pixelsize
            #read flux_units

    def __load_from_fits_file(self, filename):
        hdulist = fits.open(filename)
        self.header_0 = hdulist[0]
        self.header_1 = hdulist[1]

        self.flux = hdulist[1].data.astype(np.float64)
        self.stat = hdulist[2].data.astype(np.float64)
        pass

    def __getitem__(self, item):
        flux = self.flux.__getitem__(item)
        stat = self.stat.__getitem__(item)

        if isinstance(flux, np.ndarray):

            if flux.ndim == 2:
                return self.__to_obj(Image, flux, stat)
            elif flux.ndim == 1:
                return Spectrum
            else:
                raise ValueError(f"Invalid Dimensions for index, {item}")
        else:
            return flux

    def __add__(self, other):
        if (np.isreal(other) and not isinstance(other, np.ndarray) and not isinstance(other, Image)): #make it more pythonic
            new_flux = self.flux + other
            new_stat = self.stat + other
        elif self.__is_operable(other):
            new_flux = self.flux + other.flux
            new_stat = self.stat + other.stat
        else:
            raise ValueError(f"Unknow type {type(other)}")

        return self.__change_flux_and_stat(new_flux, new_stat)

    def __sub__(self, other):
        if (np.isreal(other) and not isinstance(other, np.ndarray) and not isinstance(other, Image)):
            new_flux = self.flux - other
            new_stat = self.stat - other
        elif self.__is_operable(other):
            new_flux = self.flux - other.flux
            new_stat = self.stat - other.stat
        else:
            raise ValueError(f"Unknow type {type(other)}")

        return self.__change_flux_and_stat(new_flux, new_stat)

    def __mul__(self, other):
        if (np.isreal(other) and not isinstance(other, np.ndarray) and not isinstance(other, Image)):
            new_flux = self.flux - other
            new_stat = self.stat - other
        elif self.__is_operable(other):
            new_flux = self.flux*other.flux
            new_stat = self.stat*other.stat
        else:
            raise ValueError(f"Unknow type {type(other)}")

        return self.__change_flux_and_stat(new_flux, new_stat)

    def __is_operable(self, other):
        if self.flux.shape != other.flux.shape:
            raise ValueError
        if self.flux_units != other.flux_units:
            raise ValueError
        return True

    def __to_obj(self, obj, flux, stat=None):
        pixelsize = self.pixel_size
        flux_units = self.flux_units
        header_0 = self.header_0
        header_1 = self.header_1

        if obj == Image:
            return Image(flux_units=flux_units, pixelsize=pixelsize,
                         data=flux, stat=stat, header_0=header_0, header_1=header_1)

    def copy(self):
        return Image(flux_units=self.flux_units, pixelsize=self.pixel_size, data=self.flux.copy(), stat=self.stat.copy(),
                     header_1=self.header_1.copy(), header_0=self.header_0.copy()) #modify copy

    def __change_flux_and_stat(self, new_flux, new_stat):
        new_image = self.copy()
        new_image.flux = new_flux
        new_image.stat = new_stat
        return new_image







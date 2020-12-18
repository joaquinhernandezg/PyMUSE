import astropy.units as u
from astropy.io import fits
import numpy as np
from mpdaf.obj import WCS
import copy
import operator
from .base import Base


__all__ = ('Image')

class Image(Base):
    n_dims = 2
    def __init__(self, filename=None, pixelsize=0.2*u.arcsec,  flux_units=1E-20 * u.erg / u.s / u.cm ** 2 / u.angstrom, data=None, stat=None,
                 header_0=None, header_1=None):

        # agregar mas verificacion de parametros
        # verificar archivo existe
        # verificar headers son astropy headers
        # verificar pixel size and flux units are apropiated astropy units
        # verificar data y stat son np array
        
        super(Image, self).__init__(filename, data, stat, header_0, header_1)
        self.flux_units = flux_units #shoud be read_from_file
        self.pixel_size = pixelsize #should be read from file

    def _Base__set_coordinates_from_header(self):
        self.wcs = WCS(hdr=self.header_1)


    def __getitem__(self, item):
        flux = self.flux.__getitem__(item)
        stat = self.stat.__getitem__(item)

        if isinstance(flux, np.ndarray):

            if flux.ndim == 2:
                return self._to_obj(Image, flux, stat)
            elif flux.ndim == 1:
                return Spectrum
            else:
                raise ValueError(f"Invalid Dimensions for index, {item}")
        else:
            return flux

    def __add__(self, other):
        return self.__operate(other, operator.add)

    def __sub__(self, other):
        return self.__operate(other, operator.sub)

    def __mul__(self, other):
        return self.__operate(other, operator.mul)

    def __truediv__(self, other):
        return self.__operate(other, operator.truediv)

    def __operate(self, other, func):
        if (np.isreal(other) and not isinstance(other, np.ndarray) and not isinstance(other, Image)):
            new_flux = func(self.flux, other)
            new_stat = func(self.stat, other)
        elif self.__is_operable(other):
            new_flux = func(self.flux, other.flux)
            new_stat = func(self.stat, other.stat)
        else:
            raise ValueError(f"Unknow type {type(other)}")

        return self.__change_flux_and_stat(new_flux, new_stat)

    def __is_operable(self, other):
        if self.flux.shape != other.flux.shape:
            raise ValueError
        if self.flux_units != other.flux_units:
            raise ValueError
        return True

    def _to_obj(self, obj, flux, stat=None):
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
        new_image = self._to_obj(Image, new_flux, new_stat)
        return new_image








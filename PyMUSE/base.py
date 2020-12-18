import numpy as np
import astropy.units as u
import os
import shutil
from abc import ABC, abstractmethod
from astropy.io import fits


class Base(ABC):
    n_dims = None

    def __init__(self, filename=None, data=None, stat=None, flux_units=1E-20 * u.erg / u.s / u.cm ** 2 / u.angstrom,
                 header_0=None, header_1=None):

        self._flux = None
        self._stat = None

        if filename:
            self.__load_from_fits_file(filename)

        else:
            self.header_0 = header_0
            self.header_1 = header_1
            self.flux = data
            self.stat = stat

        if self.header_1:
            self.__set_coordinates_from_header()

        pass

    @property
    def shape(self):
        """
        :param self:
        :return:
        """
        return self.flux.shape

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
        if not isinstance(value, np.ndarray):
            raise ValueError(f"Flux must be a numpy array")
        if value.ndim != self.n_dims:
            raise ValueError(f"Invalid flux dimensions, got {value.ndim}, expected {self.n_dims}")
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
        # mofify when value is None
        if value is None:
            self._stat = None
        elif not isinstance(value, np.ndarray):
            raise ValueError(f"Flux must be a numpy array")
        elif value.ndim != self.n_dims:
            raise ValueError(f"Invalid flux dimensions, got {value.ndim}, expected {self.n_dims}")
        elif value.shape != self.flux.shape:
            raise ValueError(f"Stat and flux can not have different dimensions, try creating a copy instead")
        else:
            self._stat = value

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

    @abstractmethod
    def __set_coordinates_from_header(self):
        pass

    #@abstractmethod
    #def __verify_headers(self, header_0, header_1):
    #    #verify dimensions, etc
    #    pass

    @abstractmethod
    def __getitem__(self, item):
        pass


    @abstractmethod
    def copy(self):
        pass

    def write_to_fits(self, filename, overwrite=False):
        hdr_0 = fits.PrimaryHDU(header=self.header_0)
        hdr_1 = fits.ImageHDU(header=self.header_1, data=self.flux, name="DATA")
        hdr_2 = fits.ImageHDU(header=self.header_1, data=self.stat, name="STAT")

        hdulist = fits.HDUList([hdr_0, hdr_1, hdr_2])
        hdulist.writeto(filename, overwrite=overwrite)

    @abstractmethod
    def _to_obj(self, obj, flux, stat=None):
        pass

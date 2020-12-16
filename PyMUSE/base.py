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
        self.header_0 = header_0
        self.header_1 = header_1

        self._flux = None
        self._stat = None

        if filename:
            self.__load_from_fits_file(filename)

        else:
            self.flux = data
            self.stat = stat

        if self.header_1:
            self.__set_coordinates_from_header()

        pass

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

    @abstractmethod
    def write_to_fits(self):
        pass



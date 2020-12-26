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

        if filename is not None:
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

    def mask(self, mask_array, crop_nans=True):
        if not mask_array.dtype == 'bool':
            raise ValueError("Invalid Mask type, should be boolean")
        masked_flux = self.flux*mask_array
        masked_flux[masked_flux == 0.0] = np.nan

        masked_stat = self.stat * mask_array
        masked_stat[masked_stat == 0.0] = np.nan

        if crop_nans:
            #this is for crop and reduce the size of the masked cube, obtained from
            #https://stackoverflow.com/questions/25831023/numpy-crop-2d-array-to-non-nan-values
            nans = np.isnan(masked_flux[0, :, :])
            nancols = np.all(nans, axis=0)  # 10 booleans, True where col is all NAN
            nanrows = np.all(nans, axis=1)  # 15 booleans

            firstcol = nancols.argmin()  # 5, the first index where not NAN
            firstrow = nanrows.argmin()  # 7

            lastcol = len(nancols) - nancols[::-1].argmin()  # 8, last index where not NAN
            lastrow = len(nanrows) - nanrows[::-1].argmin()  #

            masked_flux = masked_flux[:, firstrow:lastrow, firstcol:lastcol]
            masked_stat = masked_stat[:, firstrow:lastrow, firstcol:lastcol]

        object_copy = self._to_obj(type(self), flux=masked_flux, stat=masked_stat)
        return object_copy


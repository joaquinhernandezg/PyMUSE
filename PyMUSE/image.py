import astropy.units as u
from astropy.io import fits
import numpy as np
from mpdaf.obj import WCS


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
        self.pixelsize = pixelsize #should be read from file

        self.header_0 = header_0
        self.header_1 = header_1

        if filename:
            self.__load_from_fits_file(filename)
        else:
            self.flux = data
            self.stat = stat

        if self.header_1:
            self.wcs = WCS(hdr=self.header_1.header)
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
        data = self.data.__getitem__(item)
        var = self.var.__getitem__(item)

        if isinstance(data, np.ndarray):

            if data.ndim == 2:
                return Imagewwwwwwwwwwwwwwwwwww
            elif data.ndim == 1:
                return Spectrum
            else:
                raise ValueError(f"Invalid Dimensions for index, {item}")
        else:
            return data







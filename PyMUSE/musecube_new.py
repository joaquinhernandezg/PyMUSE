import copy
import os

import numpy as np
import astropy.units as u

from astropy.io import fits
from mpdaf.obj import WCS

from .image import Image



__all__ = ('MuseCube')

def remove_dims_from_header(header, dim="3"):
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
                 flux_units=1E-20 * u.erg / u.s / u.cm ** 2 / u.angstrom, input_wave_cal='air', data=None, stat=None, header_0=None, header_1=None):
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
        #agregar mas verificacion de parametros

        self.wave_cal = input_wave_cal
        self.flux_units = flux_units
        self.pixel_size = pixelsize

        self.header_0 = header_0
        self.header_1 = header_1
        self.filename_white = filename_white

        if filename_cube:
            self.__load_from_fits_file(filename_cube)
        else:
            self.flux = data
            self.stat = stat

        if self.header_1:
            self.wcs = WCS(hdr=self.header_1)
            #read pixelsize
            #read flux_units

        return
        self.white_data = fits.open(self.filename_white)[1].data
        self.hdulist_white = fits.open(self.filename_white)
        self.hdulist_white_temp = fits.open(self.filename_white)
        self.white_data = np.where(self.white_data < 0, 0, self.white_data)
        self.white_data_orig = fits.open(self.filename_white)[1].data

        print("MuseCube: Ready!")

    def __load_from_fits_file(self, filename):
        hdulist = fits.open(filename)
        self.header_0 = hdulist[0].header
        self.header_1 = hdulist[1].header

        self.flux = hdulist[1].data.astype(np.float64)
        self.stat = hdulist[2].data.astype(np.float64)
        return

    def load_data(self):

        if self.filename_white is None:
            print("MuseCube: No white image given, creating one.")

            w_data = self.create_white(save=False)
            w_data = np.where(w_data == 0,np.nan, w_data)

            w_header_0 = copy.deepcopy(self.header_0)
            w_header_1 = copy.deepcopy(self.header_1)

            # These loops remove the third dimension from the header's keywords. This is neccesary in order to
            # create the white image and preserve the cube astrometry
            w_header_0 = remove_dims_from_header(w_header_0)
            w_header_1 = remove_dims_from_header(w_header_1)
            
            w_header_1['WCSAXES'] = 2
            # prepare the header
            hdu = fits.HDUList()
            hdu_0 = fits.PrimaryHDU(header=w_header_0)
            hdu_1 = fits.ImageHDU(data=w_data, header=w_header_1)
            hdu.append(hdu_0)
            hdu.append(hdu_1)
            hdu.writeto('new_white.fits', clobber=True)
            self.filename_white = 'new_white.fits'
            print("MuseCube: `new_white.fits` image saved to disk.")

    def create_white(self, new_white_fitsname='white_from_colapse.fits', stat=False, save=True):
        """
        Function that collapses all wavelengths available to produce a new white image
        :param new_white_fitsname: Name of the new image
        :param stat. Default = False. If True, the new white image will be done using the stat extension of the cube.
        :param Save. Default = True. If False, the new returned image will not be saved in the hard disk
        :return: white_image.
        """
        wave = self.wavelength
        n = len(wave)
        wv_input = [[wave[0], wave[n - 1]]]
        white_image = self.get_image(wv_input, fitsname=new_white_fitsname, stat=stat, save=save)
        return white_image

    def __getitem__(self, item):
        """
        modify for including wavelenght ranges in astropy units
        """

        flux = self.flux.__getitem__(item)
        stat = self.stat.__getitem__(item)

        if isinstance(flux, np.ndarray):
            print(flux.ndim)
            if flux.ndim == 3:
                #should modify headers
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
            #modify header for wcs
            return MuseCube(flux_units=flux_units, pixelsize=pixelsize,
                         data=flux, stat=stat, header_0=header_0, header_1=header_1)


    def get_image(self, wv_input, fitsname='new_collapsed_cube.fits', type='sum', n_figure=2, save=False, stat=False,
                  maskfile=None, inverse_mask=True):
        """
        Function used to colapse a determined wavelength range in a sum or a median type
        :param wv_input: tuple or list
                         can be a list of wavelengths or a tuple that will represent a  range
        :param fitsname: str
                         The name of the fits that will contain the new image
        :param type: str, possible values: 'sum' or 'median'
                     The type of combination that will be done.
        :param n_figure: int
                         Figure to display the new image if it is saved
        :param maskfile: str, Default = None
                         Region file containing 1 more region to mask-in or mask-out
        :param inverse_mask. Boolean, Default = True, If False, the mask will be reversed.
        :return:
        """
        if maskfile:
            r = pyregion.open(maskfile)
            n = len(r)
            masks = []
            for i in range(n):
                masks.append(self.region_2dmask(pyregion.ShapeList([r[i]])))

            mask_final = masks[0]
            for i in range(n):
                mask_final = np.logical_and(mask_final, masks[i])
            if inverse_mask:
                mask_final = np.where(~mask_final, True, False)

        sub_cube = self.sub_cube(wv_input, stat=stat)
        if type == 'sum':
            matrix_flat = np.nansum(sub_cube, axis=0)
        elif type == 'median':
            matrix_flat = np.nanmedian(sub_cube, axis=0)
        else:
            raise ValueError('Unknown type, please chose sum or median')
        if maskfile:
            matrix_flat = np.where(mask_final == 1, matrix_flat, np.nan)
            if save:
                self.__save2fits(fitsname, matrix_flat, type='white', n_figure=n_figure)

        else:
            if save:
                self.__save2fits(fitsname, matrix_flat, type='white', n_figure=n_figure)

        return matrix_flat

    @property
    def wavelength(self):
        """
        Creates the wavelength array for the spectrum. The values of dw, and limits will depend
        of the data and should be revised.
        :return: w: array[]
                 array which contain an evenly sampled wavelength range
        """
        dw = self.header_1['CD3_3']
        w_ini = self.header_1['CRVAL3']
        N = self.header_1['NAXIS3']
        w_fin = w_ini + (N - 1) * dw
        # w_aux = w_ini + dw*np.arange(0, N) #todo: check whether w_aux and w are the same
        w = np.linspace(w_ini, w_fin, N)
        # print('wavelength in range ' + str(w[0]) + ' to ' + str(w[len(w) - 1]) + ' and dw = ' + str(dw))
        return w


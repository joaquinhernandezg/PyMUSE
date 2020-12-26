import pyregion
from astropy.utils import isiterable

from PyMUSE.musecube_new import MuseCube
from PyMUSE.image import Image
from PyMUSE.spectrum import Spectrum
import numpy as np

from abc import ABC, abstractmethod

class RegionMasker(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def from_params(self, *args, **kwargs):
        pass

    @staticmethod
    def from_region_string(object, region_string):
        mask = Ds9Mask.get_2d_mask_from_region_string(object, region_string)
        return mask


class MaskBox(RegionMasker):
    def __init__(self):
        super(MaskBox, self).__init__()

    def from_params(self):
        pass



class MaskEllipse(RegionMasker):
    def __init__(self):
        super(MaskEllipse, self).__init__()

    def from_params(self, object, x_c, y_c, params, coord_system='pix', origin=0):
        if not isinstance(params, (int, float, tuple, list, np.ndarray)):
            raise ValueError('Not ready for this `radius` type.')

        if origin == 0 and coord_system == 'pix':
            x_c += 1
            y_c += 1
        if coord_system == 'wcs':
            y_c, x_c = object.wcs.sky2pix([y_c, x_c], nearest=True)

        if isinstance(params, (int, float)):
            a = params
            b = params
            theta = 0
        elif isiterable(params) and (len(params) == 3):
            a = max(params[:2])
            b = min(params[:2])
            theta = params[2]
        else:
            raise ValueError('If iterable, the length of params must be == 3; otherwise try float.')

        region_string = Ds9RegionString.ellipse_params_to_region_string(x_c, y_c, a, b, theta, color="green")
        mask = self.from_region_string(object, region_string)
        return object.mask(mask)



class MaskPolygon(RegionMasker):
    def __init__(self):
        super(MaskPolygon, self).__init__()

    def from_params(self):
        pass


class Ds9Mask:
    @staticmethod
    def get_2d_mask_from_region_string(object, region_string, inverse=True):
        from pyregion.region_to_filter import as_region_filter

        if isinstance(object, MuseCube):
            im_aux = np.ones_like(object[0].flux)
        elif isinstance(object, Image):
            im_aux = np.ones_like(object.flux)
        else:
            raise ValueError(f"Invalid type {type(object)}")

        header = object.header_1
        r = pyregion.parse(region_string).as_imagecoord(header)
        shape = im_aux.shape
        region_filter = as_region_filter(r, origin=1)
        mask_new = region_filter.mask(shape)
        if inverse:
            mask_new = np.where(~mask_new, True, False)
        return mask_new

    def from_regfile(self):
        pass


class Ds9RegionString(object):
    @staticmethod
    def ellipse_params_to_region_string(x_c, y_c, a, b, theta, color="green"):
        #coord_system only in pixels

        x_center, y_center, radius = x_c, y_c, [a, b, theta]
        region_string = 'physical;ellipse({},{},{},{},{}) # color = {}'.format(x_center, y_center, radius[0],
                                                                               radius[1],
                                                                               radius[2], color)
        return region_string



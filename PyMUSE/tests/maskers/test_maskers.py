from PyMUSE.region_masker import Ds9RegionString, RegionMasker, MaskEllipse
from PyMUSE.musecube_new import  MuseCube
from PyMUSE.image import Image

import os

import unittest

import numpy as np

cube_filename = os.path.join(os.path.dirname(__file__), "../minicube.fits")



class TestMaskEllipse(unittest.TestCase):
    def setUp(self):
        self.cube = MuseCube(cube_filename)
        shape = self.cube.shape

        self.x_c = int(shape[2] / 2)
        self.y_c = int(shape[1] / 2)

        a = np.random.randint(0, int(shape[2] / 2) - 1, 1)[0]
        b = np.random.randint(0, int(shape[1] / 2) - 1, 1)[0]
        theta = 0
        self.params = [a, b, theta]

    def test_mask_ellipse_from_params_with_appropiated_params_returns_cube(self):

        masked_cube = MaskEllipse().from_params(self.cube, x_c =self.x_c, y_c=self.y_c, params=self.params)
        self.assertIsInstance(masked_cube, MuseCube)

    def test_mask_ellipse_from_params_with_appropiated_params_returns_image(self):
        image = self.cube[0]
        masked_image = MaskEllipse().from_params(image, x_c=self.x_c, y_c=self.y_c, params=self.params)
        self.assertIsInstance(masked_image, Image)

    def test_mask_ellipse_from_params_with_cube_changes_flux_in_cube(self):
        masked_cube = MaskEllipse().from_params(self.cube, x_c =self.x_c, y_c=self.y_c, params=self.params)
        self.assertFalse((masked_cube.flux == self.cube.flux).all())

    def test_mask_ellipse_from_params_with_cube_changes_flux_in_image(self):
        image = self.cube[0]
        masked_image = MaskEllipse().from_params(image, x_c=self.x_c, y_c=self.y_c, params=self.params)
        self.assertFalse((masked_image.flux == image.flux).all())

    def test_mask_ellipse_from_params_changes_stat(self):
        shape = self.cube.shape
        x_c = int(shape[2] / 2)
        y_c = int(shape[1] / 2)

        a = np.random.randint(0, int(shape[2] / 2) - 1, 1)[0]
        b = np.random.randint(0, int(shape[1] / 2) - 1, 1)[0]
        theta = 0

        masked_cube = MaskEllipse().from_params(self.cube, x_c=x_c, y_c=y_c, params=[a, b, theta])
        self.assertFalse((masked_cube.stat == MuseCube.stat).all())


if __name__ == '__main__':
    unittest.main()
import unittest
from PyMUSE.musecube_new import MuseCube
from PyMUSE.image import Image
import numpy as np


class TestCubeFromFile(unittest.TestCase):
    def setUp(self):
        self.filename = "minicube.fits"
        self.cube = MuseCube(self.filename)
        self.image = self.cube[0]

    def testImageHasflux(self):
        self.assertTrue(hasattr(self.image, "flux"))

    def testImageHasStat(self):
        self.assertTrue(hasattr(self.image, "stat"))

    def testImageHasheader_0(self):
        self.assertTrue(hasattr(self.image, "header_0"))

    def testImageHasheader_1(self):
        self.assertTrue(hasattr(self.image, "header_1"))

    def testImageHasfluxhasFlux_units(self):
        self.assertTrue(hasattr(self.image, "flux_units"))

    def testImageHasfluxhasPixel_size(self):
        self.assertTrue(hasattr(self.image, "pixel_size"))

class TestGetItem(unittest.TestCase):
    def setUp(self) -> None:
        self.filename = "minicube.fits"
        self.cube = MuseCube(self.filename)
        self.image = self.cube[0]

    def test_get_item_slice_return_a_value(self):
        value = self.image[0, 0]
        self.assertIsInstance(value, np.float)

    def test_slice_return_an_image(self):
        # corregir en la clase
        subimage = self.image[1:3, 1:3]
        self.assertIsInstance(subimage, Image)


class TestIter(unittest.TestCase):
    def setUp(self):
        self.filename = "minicube.fits"
        self.cube = MuseCube(self.filename)


if __name__ == '__main__':
    unittest.main()

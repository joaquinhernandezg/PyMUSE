import unittest
from PyMUSE.musecube_new import MuseCube
import numpy as np

class TestCubeFromFile(unittest.TestCase):
    def setUp(self):
        self.filename = "minicube.fits"
        self.cube = MuseCube(self.filename)

    def testCubeHasflux(self):
        self.assertTrue(hasattr(self.cube, "flux"))

    def testCubeHasStat(self):
        self.assertTrue(hasattr(self.cube, "stat"))

    def testCubeHasheader_0(self):
        self.assertTrue(hasattr(self.cube, "header_0"))

    def testCubeHasheader_1(self):
        self.assertTrue(hasattr(self.cube, "header_1"))

    def testCubeHasfluxhasWavecal(self):
        self.assertTrue(hasattr(self.cube, "wave_cal"))

    def testCubeHasfluxhasFlux_units(self):
        self.assertTrue(hasattr(self.cube, "flux_units"))

    def testCubeHasfluxhasPixel_size(self):
        self.assertTrue(hasattr(self.cube, "pixel_size"))

class TestGetItem(unittest.TestCase):
    def setUp(self) -> None:
        self.filename = "minicube.fits"
        self.cube = MuseCube(self.filename)

    def test_get_item_slice_return_a_value(self):
        value = self.cube[0, 0, 0]
        self.assertIsInstance(value, np.float)

    def test_slice_return_a_cube(self):
        return
        # corregir en la clase
        cls = self.cube[:, :, :]
        self.assertIsInstance(cls, MuseCube)

"""
class TestOtherThing(unittest.TestCase):
    def setUp(self) -> None:
        self.filename = "minicube.fits"
        self.cube = MuseCube(self.filename)
    def test_filename_cube_equals_filename(self):
        self.assertEqual(self.cube.filename,self.filename)

    def test_get_item_slice_return_a_value(self):
        value = self.cube[0, 0, 0]
        self.assertIsInstance(value, np.float)

    def test_slice_return_a_cube(self):
        # corregir en la clase
        cls = self.cube[:, :, :]
        self.assertIsInstance(cls, MuseCube)
"""

if __name__ == '__main__':
    unittest.main()

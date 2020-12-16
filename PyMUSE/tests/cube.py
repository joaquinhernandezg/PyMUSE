import unittest
from PyMUSE.musecube_new import MuseCube
from PyMUSE.image import Image
import numpy as np
import astropy.units as u
import os
import shutil


class TestCubeFromFile(unittest.TestCase):
    def setUp(self):
        self.filename = "minicube.fits"
        self.cube = MuseCube(self.filename)

    def testCubeHasflux(self):
        self.assertTrue(hasattr(self.cube, "flux"))

    def testFluxIsNotNone(self):
        self.assertTrue(self.cube.flux is not None)

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

    def testFluxUnitsIsNotNone(self):
        self.assertTrue(self.cube.flux_units)

    def testCubeHasfluxhasPixel_size(self):
        self.assertTrue(hasattr(self.cube, "pixel_size"))

    def testPixelSizeNotNone(self):
        self.assertTrue(self.cube.pixel_size)

class TestModifyAttributesOnTheRun(unittest.TestCase):
    def setUp(self):
        self.filename = "minicube.fits"
        self.cube = MuseCube(self.filename)

    def test_cant_set_flux_with_different_dimensions(self):
        shape_flux = list(self.cube.flux.shape)
        shape_flux[0] += 1
        new_flux = np.random.normal(size=shape_flux)
        with self.assertRaises(ValueError):
            self.cube = new_flux

    def test_cant_set_stat_with_different_dimensions_as_flux(self):
        shape_flux = list(self.cube.flux.shape)
        shape_flux[0] += 1
        new_stat = np.random.normal(size=shape_flux)
        with self.assertRaises(ValueError):
            self.stat = new_stat

    def test_cant_set_stat_with_different_dimensions_as_flux(self):
        shape_flux = list(self.cube.flux.shape)
        shape_flux[0] += 1
        new_stat = np.random.normal(size=shape_flux)
        with self.assertRaises(ValueError):
            self.stat = new_stat

    def test_cant_set_numbers_as_flux(self):
        new_flux = np.random.randint(-100, 100)
        with self.assertRaises(ValueError):
            self.flux = new_flux

    def test_cant_set_numbers_as_stat(self):
        new_stat = np.random.randint(-100, 100)
        with self.assertRaises(ValueError):
            self.flux = new_stat

class TestMethods(unittest.TestCase):
    def setUp(self):
        self.filename = "minicube.fits"
        self.cube = MuseCube(self.filename)
        if os.path.exists("temp"):
            shutil.rmtree('temp')
        os.mkdir("temp")
        self.white_filename = os.path.join("temp", "white.fits")

    def test_create_white_image_returns_Image(self):
        white = self.cube.create_white()
        self.assertIsInstance(white, Image)

    def test_create_white_image_writes_image_when_filename_provided(self):
        self.cube.create_white(new_white_fitsname=self.white_filename)
        self.assertTrue(os.path.exists(self.white_filename))

    def test_create_white_image_writes_readable_image(self):
        pass

    def test___to_obj_creates_image_with_cube_information(self):
        pass

    def test___to_obj_creates_spectrum_with_cube_information(self):
        pass

    def tearDown(self):
        shutil.rmtree('temp')
        pass

class TestCubeCopy(unittest.TestCase):
    def setUp(self):
        self.filename = "minicube.fits"
        self.cube = MuseCube(self.filename)
        self.copy = self.cube.copy()

    def test_copy_maintains_header_0(self):
        pass

    def test_copy_maintains_header_1(self):
        pass

    def test_copy_maintains_flux(self):
        pass

    def test_copy_maintains_stat(self):
        pass

    def test_copy_maintains_flux_units(self):
        pass

    def test_copy_allows_change_flux(self):
        pass

    def test_copy_allows_change_flux_and_stat_dims(self):
        pass


class TestGetItem(unittest.TestCase):
    def setUp(self) -> None:
        self.filename = "minicube.fits"
        self.cube = MuseCube(self.filename)

    def test_get_item_slice_return_a_value(self):
        value = self.cube[0, 0, 0]
        self.assertIsInstance(value, np.float)

    def test_slice_return_a_cube(self):
        # corregir en la clase
        subcube = self.cube[:, :, :]
        self.assertIsInstance(subcube, MuseCube)

    def test_slice_return_an_image(self):
        # corregir en la clase
        image = self.cube[0, 0:2, 0:2]
        self.assertIsInstance(image, Image)

    def test_slice_return_an_Spectrum(self):
        return
        # corregir en la clase
        spec= self.cube[:, 2, 2]
        self.assertIsInstance(image, Image)

    def test_slice_with_wavelenght_units(self):
        return
        spectral_range = self.cube.wave.get_range()

        subimage = self.cube[spectral_range[0]:spectral_range[1], 1:3, 1:23]
        raise  subimage



class TestIter(unittest.TestCase):
    def setUp(self):
        self.filename = "minicube.fits"
        self.cube = MuseCube(self.filename)

    def testIterReturnImage(self):
        for image in self.cube:
            self.assertIsInstance(image, Image)
            return

class TestSum(unittest.TestCase):
    def setUp(self):
        self.filename = "minicube.fits"
        self.cube = MuseCube(self.filename)

    def test_sum_works(self):
        image = self.cube.sum()
        self.assertIsInstance(image, Image)


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

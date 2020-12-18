import unittest
from PyMUSE.musecube_new import MuseCube
from PyMUSE.image import Image
import numpy as np
import astropy.units as u
import os
import shutil


class TestCubeFromFile(unittest.TestCase):
    """
    Tests presence of attributes when reading from file
    """
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

    def testCubeHasfluxhasWCS(self):
        self.assertTrue(hasattr(self.cube, "wcs"))

    def testCubeHasfluxhasWAVE(self):
        self.assertTrue(hasattr(self.cube, "wave"))

    def testFluxUnitsIsNotNone(self):
        self.assertTrue(self.cube.flux_units)

    def testCubeHasfluxhasPixel_size(self):
        self.assertTrue(hasattr(self.cube, "pixel_size"))

    def testPixelSizeNotNone(self):
        self.assertTrue(self.cube.pixel_size)

class TestModifyAttributesOnTheRun(unittest.TestCase):
    """
    Test changesa of attributes on the run
    """
    def setUp(self):
        self.filename = "minicube.fits"
        self.cube = MuseCube(self.filename)

    def test_cant_set_flux_with_different_dimensions(self):
        shape_flux = list(self.cube.flux.shape)
        shape_flux[0] += 1
        new_flux = np.random.normal(size=shape_flux)
        with self.assertRaises(ValueError):
            self.cube.flux = new_flux

    def test_cant_set_stat_with_different_dimensions_as_flux(self):
        shape_flux = list(self.cube.flux.shape)
        shape_flux[0] += 1
        new_stat = np.random.normal(size=shape_flux)
        with self.assertRaises(ValueError):
            self.cube.stat = new_stat

    def test_cant_set_stat_with_different_shape_as_flux(self):
        shape_flux = list(self.cube.flux.shape)
        shape_flux[0] += 4
        new_stat = np.random.normal(size=tuple(shape_flux))
        with self.assertRaises(ValueError):
            self.cube.stat = new_stat

    def test_cant_set_numbers_as_flux(self):
        new_flux = np.random.randint(-100, 100, 1)[0]
        with self.assertRaises(ValueError):
            self.cube.flux = new_flux

    def test_cant_set_numbers_as_stat(self):
        new_stat = np.random.randint(-100, 100, 1)[0]
        with self.assertRaises(ValueError):
            self.cube.flux = new_stat

class TestMethods(unittest.TestCase):
    """
    Test generic methods
    """
    def setUp(self):
        self.filename = "minicube.fits"
        self.cube = MuseCube(self.filename)

        self.temp_path = os.path.join(os.path.dirname(__file__), "temp")

        os.makedirs(self.temp_path)

        self.white_filename = os.path.join(self.temp_path, "white.fits")

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

    def test_write_fits(self):
        new_filename = os.path.join(self.temp_path, "cube_2.fits")
        self.cube.write_to_fits(new_filename, overwrite=True)
        self.assertTrue(os.path.exists(new_filename))

    def tearDown(self):
        shutil.rmtree(self.temp_path)
        pass

class TestCubeCopy(unittest.TestCase):
    """
    Test copy method
    """
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
        """
        cube[int, int, int] returns a number, the value of the flux in the corresponding voxel
        """
        x = np.random.randint(0, self.cube.shape[2], 1)[0]
        y = np.random.randint(0, self.cube.shape[1], 1)[0]
        z = np.random.randint(0, self.cube.shape[0], 1)[0]
        value = self.cube[z, y, x]
        self.assertIsInstance(value, np.float)

    def test_slice_return_a_cube(self):
        """
        cube[:, :, :] returns a MuseCube object
        """

        x_min, x_max = np.sort(np.random.randint(0, self.cube.shape[2], 2))
        y_min, y_max = np.sort(np.random.randint(0, self.cube.shape[1], 2))
        z_min, z_max = np.sort(np.random.randint(0, self.cube.shape[0], 2))

        subcube = self.cube[z_min:z_max, y_min:y_max, x_min:x_max]
        self.assertIsInstance(subcube, MuseCube)

    def test_slice_return_an_image(self):
        """
        cube[int, int:int, int:int] with appropiated values returns an Image object
        """
        # corregir en la clase

        x_min, x_max = np.sort(np.random.randint(0, self.cube.shape[2], 2))
        y_min, y_max = np.sort(np.random.randint(0, self.cube.shape[1], 2))
        z = np.random.randint(0, self.cube.shape[0], 1)[0]
        image = self.cube[z, y_min:y_max, x_min:x_max]
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

    def test_raises_error_when_values_outside_limits(self):
        x_min, x_max = np.sort(np.random.randint(self.cube.shape[2], self.cube.shape[2]+100, 2))
        y_min, y_max = np.sort(np.random.randint(self.cube.shape[1], self.cube.shape[1]+100, 2))
        z_min, z_max = np.sort(np.random.randint(self.cube.shape[0], self.cube.shape[0]+100, 2))

        with self.assertRaises(Exception):
            print(self.cube[z_min:z_max, y_min:y_max, x_min:x_max])



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



if __name__ == '__main__':
    unittest.main()
